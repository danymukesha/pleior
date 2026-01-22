# GWAS Summary Statistics Subset

A curated subset of GWAS catalog summary statistics prepared for
pleiotropy analysis.

## Usage

``` r
gwas_subset
```

## Format

A data.table with 1000-5000 rows and 5 columns:

- SNPS:

  Character. SNP identifiers (rs numbers)

- MAPPED_TRAIT:

  Character. Standardized trait or phenotype name

- PVALUE_MLOG:

  Numeric. -log10 transformed p-values for significance

- CHR_ID:

  Character. Chromosome identifier (1-22, X, Y)

- CHR_POS:

  Numeric. Chromosomal position in base pairs

## Source

GWAS Catalog associations file (2024 release), downloaded from:
<https://www.ebi.ac.uk/gwas/api/search/downloads/alternative>

The subset was created by: 1. Loading the complete GWAS catalog
associations file 2. Filtering for genome-wide significant associations
(P \< 1e-5) 3. Selecting top traits by number of associations 4.
Creating a balanced subset of ~500 associations 5. Ensuring
representation across multiple chromosomes

## Details

This dataset contains genome-wide association study (GWAS) summary
statistics from the GWAS Catalog (https://www.ebi.ac.uk/gwas/), focusing
on diverse traits to demonstrate pleiotropy detection and visualization.

This dataset includes associations with diverse traits such as:

- Body mass index

- Alzheimer's disease

- Type 2 diabetes

- Coronary artery disease

- Blood pressure traits

- Lipid levels (LDL, HDL cholesterol)

The data is suitable for:

- Demonstrating pleiotropy detection (SNPs associated with multiple
  traits)

- Testing visualization functions

- Teaching GWAS analysis workflows

- Benchmarking analysis methods

## References

GWAS Catalog Team. (2024). GWAS Catalog: a knowledgebase for genome-wide
association studies. \*Nucleic Acids Research\*, 51(D1), D944-D949.
<https://doi.org/10.1093/nar/gkab1009>

## Examples

``` r
# Load the data
data(gwas_subset)
head(gwas_subset)
#>      DISEASE/TRAIT CHR_ID   CHR_POS REPORTED GENE(S)       SNPS P-VALUE
#>             <char> <char>    <char>           <char>     <char>  <char>
#> 1: Body mass index      1 201815159               NR  rs2820292   8E-11
#> 2: Body mass index      2    630323               NR  rs6725549   1E-74
#> 3: Body mass index      2  24927427               NR rs10182181   2E-30
#> 4: Body mass index      2  26732753               NR  rs3739081    2E-9
#> 5: Body mass index      2  58630284               NR rs13011109   1E-14
#> 6: Body mass index      2  59078490               NR  rs1016287   4E-13
#>    PVALUE_MLOG P-VALUE (TEXT)    MAPPED_TRAIT
#>          <num>         <char>          <char>
#> 1:    10.09691           <NA> body mass index
#> 2:    74.00000           <NA> body mass index
#> 3:    29.69897           <NA> body mass index
#> 4:     8.69897           <NA> body mass index
#> 5:    14.00000           <NA> body mass index
#> 6:    12.39794           <NA> body mass index

# Basic statistics
n_unique_snps <- length(unique(gwas_subset$SNPS))
n_unique_traits <- length(unique(gwas_subset$MAPPED_TRAIT))
cat("Unique SNPs:", n_unique_snps, "\n")
#> Unique SNPs: 452164 
cat("Unique traits:", n_unique_traits, "\n")
#> Unique traits: 16452 

# Preprocess and detect pleiotropy
gwas_clean <- preprocess_gwas(gwas_subset, pvalue_threshold = 1e-5)
pleio_results <- detect_pleiotropy(gwas_clean)
head(pleio_results)
#> # A tibble: 6 × 7
#>   SNPS       N_TRAITS TRAITS             MAPPED_TRAIT PVALUE_MLOG CHR_ID CHR_POS
#>   <chr>         <int> <chr>              <chr>              <dbl> <chr>  <chr>  
#> 1 esv3585367        2 rheumatoid arthri… rheumatoid …        27   1      173496…
#> 2 esv3585367        2 rheumatoid arthri… rheumatoid …        32.4 1      173496…
#> 3 esv3585367        2 rheumatoid arthri… rheumatoid …        29.7 1      173496…
#> 4 esv3585367        2 rheumatoid arthri… rheumatoid …        24.7 1      173496…
#> 5 esv3585367        2 rheumatoid arthri… rheumatoid …        10.7 1      173496…
#> 6 esv3585367        2 rheumatoid arthri… rheumatoid …         8   1      173496…
```
