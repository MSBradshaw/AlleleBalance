import re
import pandas as pd
import sqlite3
import numpy as np
import argparse
import glob

class DatabaseAccess:
    def __init__(self, db_path):
        self.connection = sqlite3.connect(db_path)
        create_chromosome_table_template = """ CREATE TABLE IF NOT EXISTS {table}
        (
        position INTEGER PRIMARY KEY,
        genotype TEXT,
        alleles TEXT,
        allele_counts TEXT,
        allele_balances TEXT,
        first_alleles TEXT,
        reference_allele TEXT
        );
        """
        # create a table for each chromosome. the names for tables with be chromosome_1, chromosome_x etc...
        chromosomes = [str(x) for x in range(1, 24)]
        chromosomes.append('x')
        chromosomes.append('y')
        curs = self.connection.cursor()
        for c in chromosomes:
            statement = create_chromosome_table_template.format(table='chromosome_' + c)
            print(statement)
            log_file.write(statement+'\n')
            a = curs.execute(statement)

    def print_tables(self):
        curs = self.connection.cursor()
        curs.execute("SELECT name FROM sqlite_master WHERE type='table';")
        print(curs.fetchall())

    def update_db_with_tsv(self, path_to_tsv: str):
        """
        Given a tsv, add its information into the database
        :param path_to_tsv: path to a tsv (typically the output of analyze_sample.py) with these columns:
            chromosome	position	genotype	diff	allele_balance	geno_string	number_of_alleles
        :return: None
        """
        df = pd.read_csv(path_to_tsv, sep='\t', skiprows=1)
        get_template = """SELECT allele_balances, allele_counts, alleles, first_alleles 
        FROM chromosome_{chrom} WHERE position = {pos} AND genotype = "{genotype}"
        """
        update_template = """ UPDATE chromosome_{chrom}
        SET 
        allele_balances = "{ab}",
        allele_counts = "{ac}",
        alleles = "{alleles}",
        first_alleles = "{ma}"
        WHERE
        position = {pos}
        AND
        genotype = "{genotype}"
        """
        insert_template = """INSERT INTO chromosome_{chrom} 
        (position, allele_balances, allele_counts, genotype, first_alleles, alleles)
        VALUES ({pos},"{ab}","{ac}","{genotype}","{ma}","{alleles}")
        """
        # for each chromosome
        for c in set(df['chromosome']):
            sub = df[df['chromosome'] == c]
            sub.index = sub['position']
            # for each position
            curs = self.connection.cursor()
            for p in sub['position']:
                allele_balance = sub.loc[p]['allele_balance']
                first_alleles = sub.loc[p]['first_allele']
                allele_count = sub.loc[p]['allele_counts']
                genotype = sub.loc[p]['genotype']
                alleles = sub.loc[p]['alleles']
                curs.execute(get_template.format(chrom=c, pos=p, genotype=genotype))
                res = curs.fetchall()
                if len(res) == 0:
                    # there is no information at this position and genotype
                    # use the insert template
                    statement = insert_template.format(chrom=c, pos=p, ab=str(allele_balance),
                                                       ac=str(allele_count), genotype=genotype, ma=first_alleles,
                                                       alleles=alleles)
                    curs.execute(statement)
                else:
                    # use the update template
                    statement = update_template.format(chrom=c,
                                                       pos=p,
                                                       ab=res[0][0] + ';' + str(allele_balance),
                                                       ac=res[0][1] + ';' + str(allele_count),
                                                       genotype=genotype,
                                                       ma=res[0][2] + ';' + first_alleles,
                                                       alleles=res[0][3] + ';' + alleles)
                    curs.execute(statement)
            # commit changes after each chromosome
            self.connection.commit()

    def get_mean_and_std_dev(self, chromosome: str, position: int, genotype: str) -> [float, float]:
        """
        Get the mean and standard deviation of allele balance at a certain position
        :param chromosome: string, single character chromosome '1','2' ... 'x', 'y'
        :param position: integer, numerical position on the chromosome
        :param genotype: string, genotype of interest
        :return: list, [allele balance mean, allele balance standard deviation] or
                 [-1,-1] if there is no information about the desired location
        """
        get_template = """SELECT allele_balances, allele_counts, alleles, first_alleles 
        FROM chromosome_{chrom} WHERE position = {pos} AND genotype = "{genotype}"
        """
        curs = self.connection.cursor()
        statement = get_template.format(chrom=chromosome, pos=position, genotype=genotype)
        curs.execute(statement)
        res = curs.fetchall()
        # if there are no results, return [-1, -1]
        if len(res) == 0:
            return [-1, -1]
        # allele balance is the first item returned from the query
        # here I assume there is only 1 entry that matches the query
        ab = res[0][0]
        # split up the string into individual allele balances (';' denotes separation between samples)
        allele_balances = [float(x) for x in ab.split(';')]
        mean = sum(allele_balances) / len(allele_balances)
        std = np.std(allele_balances)
        return [mean, std]


def collect_and_output_genotypes(input_file, output_file, db_path, coverage_threshold, geno_dict):
    ab_info = {'chromosome': [],
               'position': [],
               'genotype': [],
               'alleles': [],
               'allele_counts': [],
               'allele_balance': [],
               'num_alleles': [],
               'first_allele': [],
               'reference_allele': [],
               'z_score': [],
               'strand_bias': []}
    da = DatabaseAccess(db_path)
    count = 0
    geno_count = 0
    geno_count_unknown = 0
    not_enough_count = 0
    for line in open(input_file, 'r'):
        count += 1
        row = line.split('\t')

        if row[0] != '1':
            break

        ref_allele = row[2]
        seq = row[4]
        # if there if not enough coverage, skip it
        if len(seq) < coverage_threshold:
            continue

        # remove indel information about bases to follow
        seq = re.sub('\\+[0-9]+[ACGTNacgtn]', '', seq)
        seq = re.sub('-[0-9]+[ACGTNacgtn]', '', seq)

        # remove > and < as these represent skips to the reference genome
        seq = seq.replace('>', '').replace('<', '')
        # remove * which denote gap in the read
        seq = seq.replace('*', '')
        # remove $ which denote end of a read
        seq = seq.replace('$', '')
        # remove beginning of read marker and the quality character that follows
        seq = re.sub('\\^.', '', seq)

        # what portion of reads come from the forward stand?
        strand_bias = sum(x == '.' or x.isupper() for x in seq) / len(seq)

        seq = seq.upper()
        seq = seq.replace(',', ref_allele)
        seq = seq.replace('.', ref_allele)

        # count of the occurrence of each allele
        seq_list = [char for char in seq]
        alleles = list(set(seq_list))
        allele_counts = [seq_list.count(x) for x in alleles]

        # sort alleles alphabetically
        alleles, allele_counts = zip(*sorted(zip(alleles, allele_counts)))

        # get the genotype
        genotype = None
        if str(row[1]) in geno_dict[str(row[0])]:
            genotype = geno_dict[str(row[0])][str(row[1])]
            geno_count += 1
        else:
            genotype = row[2] + '/' + row[2]
            geno_count_unknown += 1

        # calc AB as mode count / all count
        first_allele = genotype.split('/')[0]
        try:
            index = alleles.index(first_allele)
        except ValueError:
            log_file.write('Error: genotype not found in alleles\n')
            log_file.write(str(first_allele)+'\n')
            log_file.write(str(row) + '\n\n')
            log_file.write('Error: genotype not found in alleles\n')

        # calculate allele balance
        ab = allele_counts[index] / sum(allele_counts)

        # calc score for each site
        mean_std_arr = da.get_mean_and_std_dev(row[0], int(row[1]), genotype)
        z_score = -1
        if -1 in mean_std_arr:
            not_enough_count += 1
            log_file.write('Not enough info in DB to calculate Z-score\n')
        else:
            # caluclate the z-score, z = (x - mean) / standard deviation
            z_score = (ab - mean_std_arr[0]) / mean_std_arr[1]

        ab_info['chromosome'].append(row[0])
        ab_info['position'].append(row[1])
        ab_info['reference_allele'].append(row[3])
        ab_info['allele_balance'].append(ab)
        ab_info['num_alleles'].append(len(alleles))
        ab_info['first_allele'].append(alleles[index])
        # genotype will be formatted like this: 'A/C' or 'T/T'
        ab_info['genotype'].append(genotype)
        ab_info['alleles'].append('/'.join(alleles))
        ab_info['allele_counts'].append(','.join([str(x) for x in allele_counts]))
        ab_info['z_score'].append(z_score)
        ab_info['strand_bias'].append(strand_bias)

    df = pd.DataFrame(ab_info)
    df.to_csv(output_file, sep='\t')

    out_of_range_count = sum(df['z_score'] >= 3) + sum(df['z_score'] <= -3)
    with open(output_file, 'w') as file:
        file.write('# number of samples out of range:\t' + str(out_of_range_count))
    df.to_csv(output_file, mode='a', header=True)

    return df


def get_genotype_dict_from_vcf(vcf: str) -> dict:
    """
    Create a dictionary of the genotype in the given vcf file. Dictionary is formated as:
    a dictionary of dictionaries with chromosome (string) as the first key and position as the second
    The genotype is sorted alphabetically for example 'A/T' or 'C/G' and NEVER 'T/A' or 'G/C'
    :param vcf: path to vcf file
    :return: dictionary of genotypes
    """
    # dictionary of dictionaries with chromosome as the first key and position as the second
    geno_dict = {}
    for line in open(vcf, 'r'):
        # skip headers
        if line[:1] == '#':
            continue
        row = line.split('\t')
        ref = row[3]
        alt = row[4]
        chrom = row[0]
        pos = row[1]
        numerical_genotype = row[9].split(':')[0]
        numerical_genotype = numerical_genotype.replace('|', '/')
        genotype = ''
        if numerical_genotype == '1/1':
            genotype = alt[0] + '/' + alt[0]
        elif numerical_genotype == '0/1':
            genotype = ref[0] + '/' + alt[0]
        elif numerical_genotype == '1/2':
            alleles = row[4].split(',')
            alleles = [x[0] for x in alleles]
            # reverse soring so they get joined in the right order
            alleles.sort(reverse=True)
            genotype = '/'.join(alleles)
        else:
            log_file.write('Unknown genotype in get_genotype_dict_from_vcf:\t')
            log_file.write(str(numerical_genotype) + '\n')
            print('ERROR')
            print(numerical_genotype)
            genotype = 'ERROR'
        # if a the chromosome has not been seen before, create a new entry for it
        if chrom not in geno_dict:
            geno_dict[chrom] = {}
        # add the genotype in for the chromosome and position
        geno_dict[chrom][pos] = genotype
    return geno_dict


def run_sample(db_path):
    """
       Required Files
       1. sample.pileup
       2. sample.vcf
       3. Path to DB, as sys.argv[1]
       Output
       1. ab_report.txt
       """

    pileup_file = 'sample.pileup'
    vcf_file = 'sample.vcf'
    output_file = 'output.tsv'

    gd = get_genotype_dict_from_vcf(vcf_file)
    df = collect_and_output_genotypes(pileup_file, output_file, db_path,40, gd)

    out_of_range_count = sum(df['z_score'] >= 3) + sum(df['z_score'] <= -3)
    with open(output_file, 'w') as file:
        file.write('# number of samples out of range:\t' + str(out_of_range_count) + '\n')
    df.to_csv(output_file, mode='a', header=True, sep='\t')


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--db",
                        dest="db_path",
                        required=True,
                        help="path to exome or genome sqlite db")

    parser.add_argument("--pileup",
                        dest="pileup",
                        required=True,
                        help="path pileup file")

    parser.add_argument("--vcf",
                        dest="vcf",
                        required=True,
                        help="path vcf file")

    parser.add_argument("--type",
                        dest="run_type",
                        required=True,
                        help="function to be run: analyze sample `--type analyze` or update db `--type update`")

    args = parser.parse_args()

    return args


def check_update(files):
    files_over_thresh = 0
    for file in files:
        # the first line of each file has the number of sites out of range in each file
        for line in open(file, 'r'):
            out_of_range_count = line.split('\t')[1]
            if out_of_range_count > 0:
                files_over_thresh += 1
            break
    return files_over_thresh == len(files)


with open('allele_balance_log.txt', 'w') as log_file:
    if __name__ == "__main__":
        """
        Required Files
        1. sample.pileup
        2. sample.vcf
        3. Path to DB, as sys.argv[1]
        Output
        1. ab_report.txt
        2. allele_balance_log.txt
        """
        args = vars(get_args())
        if args['run_type'] == 'update':
            da = DatabaseAccess(args['db_path'])
            files = glob.glob('work/*/*/*.tsv')
            if check_update(files):
                for file in files:
                    da.update_db_with_tsv(file)
                    log_file.write('Updated database with file:' + file + ' \n')
            else:
                log_file.write('Not updating database, too many samples with an abnormal number of imbalanced '
                               'alleles. Bad batch suspected\n')
        elif args['run_type'] == 'analyze':
            run_sample(args['db_path'])
        else:
            print('Invalid `--type`. To analyze a sample use `--type analyze` or to update the db `--type update`')
            log_file.write('Invalid `--type`. To analyze a sample use `--type analyze` or to update the db `--type update`\n')
            quit()
