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

---

**NOTE**

This is production code! Use it at your own risk. That means that documentation may be incomplete and functions may return unexpected errors.

<strike>[Update September 2019]

Due to some substantial changes in functions/methods from some R/Bioconductor packages that RNAModR depends on, RNAModR's functionality is currently limited. Specifically,

1. `BuildTx` (building a custom reference transcriptome) is broken (but you can still use the pre-generated transcriptomes from the links below); this seems to be related to major [changes in `GenomicRanges` version 1.32.0](https://github.com/Bioconductor/GenomicRanges/blob/master/NEWS) and changes in `RMariaDB` that completely [break functionality of `GenomicFeatures::makeTxDbFromUCSC`](https://github.com/r-dbi/RMariaDB/issues/135) (verified on MacOS Sierra). I am unsure how badly other OS are affected. A workaround/fix for MacOS Sierra (that should also work for other flavours) has been posted as part of [Issue #5](https://github.com/mevers/RNAModR/issues/5).
2. `plotRelDistDistribution` and `plotRelDistEnrichment()` are broken; this is related to afore-mentioned changes in `GenomicRanges`, which introduced a `CompressedGRangesList` class as a replacement for the `GRangesList` class.

I appreciate any and all testing; for issues, please open an official [Issue](https://github.com/mevers/RNAModR/issues/new) on the GitHub project site.
</strike>

[Update October 2019]

Functionality of most (all?) functions has been restored in a series of major code revisions. This has led to a new development version 0.2.0, which should become stable soon. Due to substantial changes in R/Bioconductor packages that `RNAModR` depends on and in the code base of `RNAModR` itself, reference transcriptomes from older `RNAModR` versions will not work with the current development version 0.2.0 of `RNAModR`. A full list of changes can be found in [NEWS](NEWS).

I appreciate and encourage any and all testing; for issues, please open an official [Issue](https://github.com/mevers/RNAModR/issues/new) on the GitHub project site.

---


## Installing RNAModR

### The github way (requires the [devtools](https://github.com/hadley/devtools) package)

1. Make sure you have the following R/Bioconductor packages installed

    * AnnotationDbi
    * beanplot
    * Biostrings
    * GenomeInfoDb
    * GenomicFeatures
    * GenomicRanges
    * gplots
    * RSQLite
    * rtracklayer

    You can install additional R/Bioconductor packages in the usual way:

    ```{r}
    source("http://www.bioconductor.org/biocLite.R")
    biocLite(c("AnnotationDbi", "beanplot", "Biostrings", "GenomeInfoDb", "GenomicFeatures", "GenomicRanges", "gplots", "RSQLite",    "rtracklayer"))
    ```

    Additionally, RNAModR requires two _organism-specific_ R packages to contruct a custom transcriptome. Currently, RNAModR    supports human and mouse data, based on the following reference genome versions

     * Human: hg38, hg19
     * Mouse: mm10, mm9.

    Please install the corresponding organism- and version-matching R/Bioconductor packages. For example, if genomic loci of RNA modification are based on the human GRCh38/hg38 reference genome,
    RNAModR requires the following R/Bioconductor packages:

     * BSgenome.Hsapiens.UCSC.hg38
     * org.Hs.eg.db

    which you can install in the usual way

    ```{r}
    biocLite(c("BSgenome.Hsapiens.UCSC.hg38", "org.Hs.eg.db"))
    ```

    We also offer the possibility to download pre-constructed transcriptome data, see section [Downloadable transcriptome data](#downloadTx).


2. If all package dependencies are met, install **RNAModR** with devtools

    ```{r}
    if (!require("devtools")) install.packages("devtools")
    devtools::install_github("mevers/RNAModR", build_vignettes = FALSE)
    ```

    You can force a re-install of RNAModR by adding

    ```{r}
    devtools::install_github(..., force = TRUE)
    ```

## Getting started
The following lines of R code will load the **RNAModR** library, and plot the distribution of m6A sites [[Linder et al., Nature Methods 12, 767 (2015)](http://www.nature.com/nmeth/journal/v12/n8/abs/nmeth.3453.html)] across the 5'UTR, CDS and 3'UTR of the human hg38-based transcriptome.

Note: You can also download pre-constructed transcriptome data, see the [next section](#downloadTx) for details.

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

## Downloadable transcriptome data<a name="downloadTx"></a>
You can download pre-constructed transcriptome data files for the following reference genome versions


*Homo sapiens*

1. [hg38](https://drive.google.com/open?id=1nBRsUWEq5FvoZmdYtGJWZhQCi2izajCr)
2. [hg19](https://drive.google.com/open?id=1OQnsmuieQw7KUXKPy6C5UnZGuavooW06)
3. [hg18](https://drive.google.com/open?id=18xufP2MQn39gTgkHob8dvOPXiPhKg_wc)

*Mus musculus*

1. [mm10](https://drive.google.com/open?id=17i3yHBjkL50K-o60mMFiP2SzhuW6v9nP)
2. [mm9](https://drive.google.com/open?id=1fO3BSojCb_BIE8DzmEKHw1miJOJZt0Zr)
3. [mm8](https://drive.google.com/open?id=1SqJEX0O6HL1baW8XHWkOMJr37AfAZ4q2)


----

Above files will work with the RNAModR version 0.1.1; for version 0.1.0 you may want to try the older transcriptome files. There's no guarantee that they will still work; major changes in `GenomicRanges` and other Bioconductor libraries may cause unpredictable behaviour and errors.

*Homo sapiens*

1. [hg38](https://drive.google.com/open?id=0B5_hfxBdKWHRVlBCTUlSazJfaWs)
2. [hg19](https://drive.google.com/open?id=0B5_hfxBdKWHRemRwMUtoa1ZnNVE)
3. [hg18](https://drive.google.com/open?id=0B5_hfxBdKWHRZ3ZqdXpVN0VmWlU)

*Mus musculus*

1. [mm10](https://drive.google.com/open?id=0B5_hfxBdKWHRYzliNkotN1NwSTQ)
2. [mm9](https://drive.google.com/open?id=0B5_hfxBdKWHRdkNSMmNleVJUSm8)
3. [mm8](https://drive.google.com/open?id=0B5_hfxBdKWHRMmN5WGRsRkpWcWc)

----

In order to use the transcriptome data you need to copy the RData file into your working directory.
You can check that RNAModR correctly finds the transcriptome data by running e.g.

```{r}
BuildTx("hg38");
```

Provided you have copied the file tx_hg38.RData into the working directory, this should produce the following message

```
Found existing transcriptome data. Nothing to do.
To rebuild run with force = TRUE.
```

## Documentation

The RNAModR manual can be downloaded [here](doc/RNAModR-manual.pdf).


## Contributors

Please contact [Maurits Evers](mailto:maurits.evers@anu.edu.au "Email Maurits Evers") in case of questions/suggestions.
In case of bugs/feature requests please open an issue on github.

## Licensing

The **RNAModR** R package is open source licensed under the
GNU Public License, version 3 (GPLv3).
