# Allele Balance Quality Control

This program is a QC metric for [Layer Lab's exome and genome analysis pipeline](https://github.com/javaidm/layer_lab_vc).

There are two main uses for this program.
1. Analyze samples for allele balance (AB)
2. Update the sqlite database for tracking AB across many samples

## Command line parameters
`--type` (required) function to be run: analyze sample `--type analyze` or update db `--type update`

`--db` (required) path to the sqlite database to be used. There should be a separate database for exomes and genomes.

`--vcf` (required if using `--type analyze`) path vcf file to be used

`--pileup` (required if using `--type analyze`) path pileup.gz (note that it is expected to be compressed) file to be used

## Analyze Sample
This function assumes that several files exist in the working directory
1. `sample.mpileup`
2. `sample.vcf`

Using the information found in `sample.mpileup` the allele balance will be calculated and compared to all previous samples in the db. 
For each base in `sample.mpileup` a [z-score](https://en.wikipedia.org/wiki/Standard_score) is calculated comparing it to all previous samples of the same genotype at that base.
z-scores cannot be calculated if there are lessons than 3 previous samples in the db at that base with the same genotype.
`sample.vcf` used to get the genotype called for each base and for nothing else.

Two files are produced by this function and saved in the working directory
1. `ab_report.txt` contains 1 header line with information about the number of sites in the sample with a z-score above 3 or below -3 below which is formated like a tsv with the following columns: '  chromosome   position  genotype alleles allele_counts   allele_balance    num_alleles    first_allele   reference_allele  z_score  strand_bias'
2. `allele_balance_log.txt` log file with errors, warnings and debugging output - probably only useful for the guy writing this documentation

## Update DB

This function takes the file given as the `--input` parameter and adds it's information to the database given in the `--db` parameter.

The tables in teh database are formatted with the following statement:

    CREATE TABLE IF NOT EXISTS chromosome_1
        (
        position INTEGER PRIMARY KEY,
        genotype TEXT,
        alleles TEXT,
        allele_counts TEXT,
        allele_balances TEXT,
        first_alleles TEXT,
        reference_allele TEXT
        );
There is one table for each chromosome (1-22, x and y)  named like `chromosome_1` or `chromosome_x`
 
## Usage Example

Analyze a sample in the current directory

`python --type analyze --db path/to/exomes.db --vcf path/to/sample.vcf --pileup path/to/sample.pileup.gz`

Update the database

`python --type update --db path/to/exomes.db --input ab_report.txt`