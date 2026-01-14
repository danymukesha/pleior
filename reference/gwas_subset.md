# Example GWAS subset data

A subset of GWAS summary statistics containing SNPs and their
associations with multiple traits, used for demonstrating pleiotropy
analysis.

## Usage

``` r
gwas_subset
```

## Format

A data.table with 10 rows and 5 columns:

- SNPS:

  Character. SNP identifiers (rs numbers)

- MAPPED_TRAIT:

  Character. Associated trait or phenotype

- PVALUE_MLOG:

  Numeric. -log10 transformed p-values

- CHR_ID:

  Character. Chromosome identifier

- CHR_POS:

  Character. Chromosomal position

## Source

Simulated data based on real GWAS catalog associations

## Examples

``` r
data(gwas_subset)
head(gwas_subset)
#>         SNPS          MAPPED_TRAIT PVALUE_MLOG CHR_ID   CHR_POS
#>       <char>                <char>       <num> <char>    <char>
#> 1:  rs814573     Alzheimer disease   672.69900     19  44908684
#> 2:  rs814573 myocardial infarction    15.00000     19  44908684
#> 3:    rs7412     Alzheimer disease   122.39790     19  44919689
#> 4:    rs7412       LDL cholesterol  9629.00000     19  44919689
#> 5: rs2817462    memory performance     6.30103      6 156588550
#> 6: rs2817462  tyrosine measurement     9.69897      6 156588550
```
