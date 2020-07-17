import sqlite3
import pandas as pd
import numpy as np


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
        df = pd.read_csv(path_to_tsv, sep='\t')
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
                # TODO when genotype calling is added this will no longer be the same information
                alleles = sub.loc[p]['genotype']
                curs.execute(get_template.format(chrom=c, pos=p, genotype=genotype))
                res = curs.fetchall()
                if len(res) == 0:
                    print('Nada')
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



# da = DatabaseAccess('ab.db')
# for i in range(4):
#     print()
#     print()
#     print()
#     print()
#     da.update_db_with_tsv('delete1.tsv')
# print(da.get_mean_and_std_dev('1', 14747, 'C/C'))
# curs = da.connection.cursor()
# curs.execute('SELECT * FROM chromosome_1')

# da.get_allele_balance_population_mean_and_std_dev('1', 0)
