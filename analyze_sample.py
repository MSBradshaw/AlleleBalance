import vcf
import re
import pandas as pd
import matplotlib.pyplot as plt
from sys import argv

"""
Parameters
1. path to g.vcf file
2. path to output dir (needs to already exist and should be unique for the sample)
"""
# argv = ['', '/Users/michael/TESTBAMs/HaplotypeCaller_GIABO1_S53.g.vcf', './']
# make sure the output dir ends with a '/'
output_dir = argv[2]
if output_dir[-1] != '/':
    output_dir += '/'

vcf_reader = vcf.Reader(open(argv[1], 'r'))

summary_stats = {'chromosome': [],
                 'position': [],
                 'genotype': [],
                 'diff': [],
                 'allele_balance': [],
                 'geno_string': [],
                 'number_of_alleles': []}
for i, record in enumerate(vcf_reader):
    if i % 1000 == 0:
        print(i)
    # check if the record has allelic depth information we need this for stats
    if 'AD' in record.samples[0].data._asdict():
        # get the genotype
        gt = record.samples[0]['GT']
        # if all allele depths are 0, skip this location
        if sum(record.samples[0]['AD']) == 0:
            continue
        # calc allele balance of reference allele / total # of alleles
        ab = record.samples[0]['AD'][0] / sum(record.samples[0]['AD'])
        # variable to store how much this site differs from the expected
        diff = -1
        # variable for storing the simple genotype (hetero or homo)
        geno_string = 'homo'
        # check if it is homozygous
        if re.search(r'(\d)\/\1', gt):
            # is it homozygous reference or not? if it is homozygous, expected allele balance is 1.0
            if gt[0] == '0':
                # is reference
                diff = 1.0 - ab
            else:
                # is alternate
                diff = 0.0 - ab
        else:
            geno_string = 'hetero'
            # if it is heterozygous, expected allele balance is 0.5
            diff = abs(0.5 - ab)
        summary_stats['chromosome'].append(record.CHROM)
        summary_stats['position'].append(record.POS)
        summary_stats['genotype'].append(gt)
        summary_stats['allele_balance'].append(ab)
        summary_stats['diff'].append(diff)
        summary_stats['geno_string'].append(geno_string)
        # the number of alleles is the number of non-zero allele counts
        summary_stats['number_of_alleles'].append(len(record.samples[0]['AD']) - record.samples[0]['AD'].count(0))
    else:
        # there is no allelic depth information so skip it
        pass

df = pd.DataFrame(summary_stats)

extreme = sum(df['diff'] >= 0.2)
extreme_hetero = sum(df[df['geno_string'] == 'hetero']['diff'] >= 0.2)
extreme_homo = sum(df[df['geno_string'] == 'homo']['diff'] >= 0.2)
moderate = sum(df['diff'] >= 0.05)
moderate_hetero = sum(df[df['geno_string'] == 'hetero']['diff'] >= 0.05)
moderate_homo = sum(df[df['geno_string'] == 'homo']['diff'] >= 0.05)

with open(output_dir + 'allelic_balance_counts.txt', 'w') as file:
    file.write('Percent heterozygous: ' + str(sum(df['geno_string'] == 'hetero') / df.shape[0]) + '\n')
    file.write('Heterozygous count: ' + str(sum(df['geno_string'] == 'hetero')) + '\n')
    file.write('Homozygous count: ' + str(sum(df['geno_string'] == 'homo')) + '\n')
    file.write('Number of extreme differences: ' + str(extreme) + '\n')
    file.write('Number of extreme heterozygous differences: ' + str(extreme_hetero) + '\n')
    file.write('Number of extreme homo differences: ' + str(extreme_homo) + '\n')
    file.write('Number of moderate differences: ' + str(moderate) + '\n')
    file.write('Number of moderate heterozygous differences: ' + str(moderate_hetero) + '\n')
    file.write('Number of moderate homo differences: ' + str(moderate_homo) + '\n')

# plot allelic depth
plt.hist(
    [df[df['geno_string'] != 'hetero']["allele_balance"], df[df['geno_string'] == 'hetero']["allele_balance"]],
    label=['Homozygous', 'Heterozygous'], bins=[x / 100 for x in range(0, 101, 2)])
plt.legend(loc='upper right')
plt.xlabel('Allele balance')
plt.ylabel('Count')
plt.savefig(output_dir + 'allele_balance.png')
plt.clf()

# plot difference of expected vs observed
plt.hist(
    [df[df['geno_string'] != 'hetero']["diff"], df[df['geno_string'] == 'hetero']["diff"]],
    label=['Homozygous', 'Heterozygous'], bins=[x / 100 for x in range(0, 101, 2)])
plt.legend(loc='upper right')
plt.xlabel('Difference from expected')
plt.ylabel('Count')
plt.savefig(output_dir + 'difference.png')
plt.clf()

# save summary stats df to a file
df.to_csv(output_dir + 'summary_stats.tsv', sep='\t')
