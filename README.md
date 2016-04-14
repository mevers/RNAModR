# RNAModR

**RNAModR** provides functions to map single-nucleotide RNA modifications to a 
reference mRNA transcriptome, and perform exploratory functional analyses of
sites across the transcriptome. All RNAModR functions are available through
an R package.

### Installing the **RNAModR** R package from Github (requires the [devtools](https://github.com/hadley/devtools) package)
```{r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mevers/RNAModR", build_vignettes = TRUE)
```
### About RNAModR

**RNAModR** provides plotting and statistical analysis routines to analyse
single-nucleotide RNA modifications across the transcriptome.

**RNAModR** allows the analysis and display of spatial and sequence-specific 
localisation of RNA modification sites. The statistical significance of
site enrichment/depletion is evaluated using multiple Fisher's exact tests.
Multiple testing corrections are applied using the method of Benjamini and
Hochberg.   

### Contributors

Please contact [maurits.evers](mailto:maurits.evers@anu.edu.au "Email Maurits Evers") in case of questions/suggestions.
In case of bugs/feature requests please open an issue on github.

### Licensing

The **RNAModR** R package is open source licensed under the 
GNU Public License, version 3 (GPLv3).
