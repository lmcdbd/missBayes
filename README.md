## OVERVIEW
missBayes is an R package for differential expression analysis of proteomics data. It applies a Bayesian framework to estimate posterior distributions and posterior probabilities of protein log2 fold changes between conditions.

A key feature of missBayes is its explicit modeling of missing values, allowing the method to account for missingness mechanisms commonly observed in bottom-up proteomics data and thereby improve inference robustness.

Details of the methods and evaluation used in this package are described in the accompanying bioRxiv preprint:
https://www.biorxiv.org/content/10.64898/2026.01.15.699650v1

## INSTALLATION
```r
install.packages("devtools")
devtools::install_github("lmcdbd/missBayes")
```
The model specification files are installed with the package; nothing further
needs to be downloaded.

## RUNNING THE MODEL

```r
library(limma)
library(parallel)
library(missBayes)
# log2all.df: a numeric matrix or data.frame of log2 intensities. Rows correspond to proteins (with protein IDs as rownames), and columns correspond to samples. Each row must contain at least one non-NA value

# metadata: a data.frame containing sample annotations. The column 'Condition' specifies experiment groups

group <- as.factor(metadata$Condition)

# Define comparison using limma::makeContrasts

comparison <- makeContrasts("Treatment - Control", levels = levels(group))

# Fit the Bayesian missing model
# output: a list of data.frames containing model results for each contrast specified in comparisons
set.seed(456)  # Set seed for reproducibility
output <- BayesMissingModel(log2all.df, group, comparison)

# Plot posterior distribution for a given protein
plotPost(log2all.df, "protein ID", group, "Treatment - Control")
