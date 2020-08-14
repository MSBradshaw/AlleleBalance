# Allele Balance Quality Control

This program is a QC metric for [Layer Lab's exome and genome analysis pipeline](https://github.com/javaidm/layer_lab_vc).

There are two main uses for this program.
1. Analyze samples for allele balance (AB)
2. Update the sqlite database for tracking AB across many samples

## Command line parameters
`--type` (required) function to be run: analyze sample `--type analyze` or update db `--type update`

`--db` (required) path to the sqlite database to be used. There should be a separate database for exomes and genomes.

`--vcf` (required if using `--type analyze`) path vcf.gz file to be used (note that it is expected to be compressed)

`--pileup` (required if using `--type analyze`) path pileup.gz (note that it is expected to be compressed) file to be used

`--inputs` (required if using `--type update`) comma separated list of file paths for `output.tsv` files generated by `--type analyze` to add to the data base. This option is only used when using `--type update`

## Analyze Sample
To use this function you must specific `--type analyze`

Required parametes:
1. `--db`
2. `--pileup`
3. `--vcf`

Using the information found in the file specified with the `--pileup` option (must be a .pileup.gz), the allele balance will be calculated and compared to all previous samples in the db. 
For each base in the pileup a [z-score](https://en.wikipedia.org/wiki/Standard_score) is calculated comparing it to all previous samples of the same genotype at that base.
z-scores cannot be calculated if there are lessons than 3 previous samples in the db at that base with the same genotype.
The file specified by `--vcf` (must be a .vcf.gz) is used to get the genotype called for each base and for nothing else.

Two files are produced by this function and saved in the working directory
1. `output.tsv` contains 1 header line with information about the number of sites in the sample with a z-score above 3 or below -3 below which is formatted like a tsv with the following columns: '  chromosome   position  genotype alleles allele_counts   allele_balance    num_alleles    first_allele   reference_allele  z_score  strand_bias'
2. `allele_balance_log.txt` log file with errors, warnings and debugging output - probably only useful for the guy writing this documentation

## Update DB
To use this function you must specific `--type update`

Required parametes:
1. `--db`
2. `--inputs`

This function takes the files specified in the `--inputs` parameter, checks if as a whole they meet the criteria for being added to the db and if the files do it updates the db. 

The current criteria for being added to the db are:
1. None of the files are over threshold
2. Threshold for an individual files is if more than 0 sites are out of range based on z-scores (so the z-score must be within -3 and 3)

The tables in the database are formatted with the following statement:

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
 
## Usage Examples

Analyze a sample in the current directory

`python --type analyze --db path/to/exomes.db --vcf path/to/sample.vcf --pileup path/to/sample.pileup.gz`

Update the database

`python --type update --db path/to/exomes.db --input results/sample_1/output.tsv,results/sample_2/output.tsv`