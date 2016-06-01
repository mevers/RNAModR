# RNAModR

**RNAModR** provides functions to map lists of genomic loci of RNA modifications 
to a reference mRNA transcriptome, and perform exploratory functional analyses of
sites across the transcriptome trough visualisation and statistical analysis of
the distribution of sites across transcriptome sections (5'UTR, CDS, 3'UTR).

**RNAModR** performs enrichment analyses to assess the statistical significance of
the spatial and sequence-specific localisation of sites relative to null sites.
Null sites can be generated automatically or may be supplied manually in the form
of a list of genomic loci.

Note that enrichment analyses results may depend critically on the choice and validity 
of null sites. For example, in establishing a list of null sites for 5-methylcytidine 
modifications, the non-uniformity in the actual distribution of cytosines across 
transcript sections makes a simple position permutation approach inappropriate for 
generating null sites.

## Installing RNAModR 

### The github way (requires the [devtools](https://github.com/hadley/devtools) package)

1. Open an R terminal.

2. Within R, type

   ```{r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mevers/RNAModR", build_vignettes = FALSE)
   ```

### Additional dependencies
**RNAmodR** makes use of functions provided by the following R/Bioconductor packages:
 
 * AnnotationDbi
 * beanplot
 * Biostrings
 * GenomeInfoDb
 * GenomicFeatures
 * GenomicRanges
 * gplots
 * RSQLite
 * rtracklayer

Additionally, **RNAModR** requires _organism-specific_ R packages to contruct a custom transcriptome. 
For example, if genomic positions of RNA modification are based on the human GRCh38/hg38 reference genome,
RNAModR requires the following packages:

* BSgenome.Hsapiens.UCSC.hg38
* org.Hs.eg.db

Currently, **RNAModR** supports human and mouse data, based on the following reference genome versions

* Human: hg38, hg19
* Mouse: mm10, mm9.

Note that additional packages can be installed following the usual Bioconductor procedure; for example, if data 
are based on the hg38 reference genome, the following R commands will install the necessary dependencies:

```{r}
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi", "beanplot", "Biostrings", "BSgenome.Hsapiens.UCSC.hg38", "GenomeInfoDb", "GenomicFeatures", "GenomicRanges", "gplots", "org.Hs.eg.db", "RSQLite", "rtracklayer"))
```

## Getting started
The following lines of R code will load the **RNAModR** library, and plot the distribution of m6A sites [[Linder et al., Nature Methods 12, 767 (2015)](http://www.nature.com/nmeth/journal/v12/n8/abs/nmeth.3453.html)] across the 5'UTR, CDS and 3'UTR of the human hg38-based transcriptome.

```{r}
# Load the library.
library(RNAModR);

# Build reference transcriptome.
# This might take a few minutes.
BuildTx("hg38");

# Load and map m6A sites to reference transcriptome.
posSites <- ReadBED(system.file("extdata", "miCLIP_m6A_Linder2015_hg38.bed", package = "RNAModR"));
posSites <- SmartMap(posSites, "m6A_Linder");

# Keep sites located in the 5'UTR, CDS and 3'UTR
posSites <- FilterTxLoc(posSites, filter = c("5'UTR", "CDS", "3'UTR"));

# Plot distribution across transcript sections
PlotSectionDistribution(posSites);
```

## Documentation

The RNAModR manual can be downloaded [here](doc/RNAModR-manual.pdf).


## Contributors

Please contact [Maurits Evers](mailto:maurits.evers@anu.edu.au "Email Maurits Evers") in case of questions/suggestions.
In case of bugs/feature requests please open an issue on github.

## Licensing

The **RNAModR** R package is open source licensed under the 
GNU Public License, version 3 (GPLv3).
