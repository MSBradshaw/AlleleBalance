import re
import pandas as pd
import matplotlib.pyplot as plt
from sys import argv

"""
Parameters
1. path to pileup file
2. path to vcf file
3. path to output dir (needs to already exist and should be unique for the sample)
"""


def collect_and_output_genotypes(input_file, output_file, coverage_threshold, geno_dict):
    ab_info = {'chromosome': [],
               'position': [],
               'genotype': [],
               'alleles': [],
               'allele_counts': [],
               'allele_balance': [],
               'num_alleles': [],
               'first_allele': [],
               'reference_allele': []}
    count = 0
    geno_count = 0
    geno_count_unknown = 0
    for line in open(input_file, 'r'):
        if count % 10000 == 0:
            pass
            # print(count)
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

        seq = seq.upper()
        seq = seq.replace(',', ref_allele)
        seq = seq.replace('.', ref_allele)

        # remove > and < as these represent skips to the reference genome
        seq = seq.replace('>', '').replace('<', '')
        # remove * which denote gap in the read
        seq = seq.replace('*', '')
        # remove $ which denote end of a read
        seq = seq.replace('$', '')
        # remove beginning of read marker and the quality character that follows
        seq = re.sub('\\^.', '', seq)

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
        print()
        print(alleles)
        print(first_allele)
        index = alleles.index(first_allele)

        ab = allele_counts[index] / sum(allele_counts)
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

    pd.DataFrame(ab_info).to_csv(output_file, sep='\t')


# collect_and_output_genotypes(argv[1], argv[3], 40)

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
            genotype = alt + '/' + alt
        elif numerical_genotype == '0/1':
            genotype = ref + '/' + alt
        elif numerical_genotype == '1/2':
            alleles = row[4].split(',')
            # reverse soring so they get joined in the right order
            alleles.sort(reverse=True)
            genotype = '/'.join(alleles)
        else:
            print(numerical_genotype)
            print('ERROR')
            genotype = 'ERROR'
        # if a the chromosome has not been seen before, create a new entry for it
        if chrom not in geno_dict:
            geno_dict[chrom] = {}
        # add the genotype in for the chromosome and position
        geno_dict[chrom][pos] = genotype
    return geno_dict


gd = get_genotype_dict_from_vcf(argv[2])
collect_and_output_genotypes(argv[1], argv[3], 40, gd)
# TODO use first allele in genotype for allele balance calculations
# TODO calc z score
# TODO add strand bias to this script
# TODO strand bias DB creation
# TODO make a summary report
