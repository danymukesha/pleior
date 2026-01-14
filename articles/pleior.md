# Introduction

*pleior* is an R package for analyzing pleiotropy in genome-wide
association studies (GWAS).

Here, we provide a basic workflow for using *pleior* and how to use its
core functions to load, preprocess, and analyze GWAS data for
pleiotropic effects.

## Setup

First, load the *pleior* package and the example dataset.

``` r
library(pleior)
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

## Workflow

### Step 1: Load GWAS data

Although we use the example dataset here, you can load your own GWAS
summary statistics.

``` r
gwas_data <- gwas_subset
```

### Step 2: Preprocess GWAS data

Filter the data to retain significant associations.

``` r
gwas_clean <- preprocess_gwas(gwas_data, pvalue_threshold = 5e-8)
head(gwas_clean)
#>          SNPS          MAPPED_TRAIT PVALUE_MLOG CHR_ID   CHR_POS
#>        <char>                <char>       <num> <char>    <char>
#> 1:   rs814573     Alzheimer disease   672.69900     19  44908684
#> 2:   rs814573 myocardial infarction    15.00000     19  44908684
#> 3:     rs7412     Alzheimer disease   122.39790     19  44919689
#> 4:     rs7412       LDL cholesterol  9629.00000     19  44919689
#> 5:  rs2817462  tyrosine measurement     9.69897      6 156588550
#> 6: rs10182181       body mass index    29.69897      2  24927427
```

### Step 3: Detect pleiotropic SNPs

Identify SNPs associated with multiple traits.

``` r
pleio_results <- detect_pleiotropy(gwas_clean, traits = c("Alzheimer disease", "myocardial infarction", "LDL cholesterol"))
head(pleio_results)
#>        SNPS N_TRAITS                                  TRAITS
#>      <char>    <int>                                  <char>
#> 1:   rs7412        2       Alzheimer disease;LDL cholesterol
#> 2:   rs7412        2       Alzheimer disease;LDL cholesterol
#> 3: rs814573        2 Alzheimer disease;myocardial infarction
#> 4: rs814573        2 Alzheimer disease;myocardial infarction
#>             MAPPED_TRAIT PVALUE_MLOG CHR_ID  CHR_POS
#>                   <char>       <num> <char>   <char>
#> 1:     Alzheimer disease    122.3979     19 44919689
#> 2:       LDL cholesterol   9629.0000     19 44919689
#> 3:     Alzheimer disease    672.6990     19 44908684
#> 4: myocardial infarction     15.0000     19 44908684
```

### Step 4: Visualize results

Create a Manhattan plot to visualize pleiotropic SNPs.

``` r
plot_pleiotropy_manhattan(pleio_results, highlight_snp = "rs814573")
```

![](pleior_files/figure-html/plot-1.png)

For more advanced features, such as functional annotation or interactive
visualization, refer to the package documentation.
