## Introduction

Welcome! The purpose of this R repository is to calculate escape genes. Implementation is based on a paper by Joel Berletch et al.:

```
Berletch JB, Ma W, Yang F, Shendure J, Noble WS, Disteche CM, et al. (2015)
Escape from X Inactivation Varies in Mouse Tissues. PLoS Genet 11(3): e1005079.
doi:10.1371/journal.pgen.1005079
```
Berletch's source data is located on GEO here: [GSE59777](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59777)

We have made some modifications to the code which are important to note: 

1. Our original version systematically overestimated the RPKMs (off by ~2 orders of magnitude) because it was systematically underestimating the `total_num_reads`, as it misinterpreted that quantity to be the sum of SNP-specific reads. The correct `total_num_reads` must be calculated on all mapped reads, not just SNP-specific reads, so it cannot be computed from the input into this pipeline. Rather, `total_num_reads` must be supplied a priori.
2. Our original version averages mice together, but in Disteche's chapter in the X Chromosome Inactivation Methods book published in 2018 (found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6269188/)), the authors suggest that "Biological replicates of RNA-seq experiments should be analyzed separately."
3. Finally, like the `total_num_reads`, the RPKMs must also be derived a priori, because they are calculated on all mapped reads, not just the SNP-specific reads. In the absence of that information, we impute the `total_num_reads` by estimating that the total SNP-specific reads is about 0.121 of the `total_num_reads`. Using this as a estimate brings the computed/estimated RPKM down to reasonable values.

These modifications are now included in the present version of this code.

The XCI escape gene analysis performed for mouse AT2 cells in PMID: [36638790](https://pubmed.ncbi.nlm.nih.gov/36638790) utilized our original version of this analytic pipeline prior to the incorporation of these changes present in this version available here on GitHub. Accordingly, a smaller number of total escape genes would be identified falling below the statistical cutoff than what was listed in the publication.


## Installation

Run the following in your R console. We have used R 4.3.0 with success, and any version above 3.6.0 should work fine.

```R
install.packages("devtools") # warning: takes a long time!
install.packages("reshape")
install.packages("seqinr")
install.packages("logr")
install.packages("data.table")
install.packages("this.path")  # See: https://github.com/ArcadeAntics/this.path
install.packages("optparse")
install.packages("tidyr")
install.packages("import")
install.packages("plotly")
install.packages("htmlwidgets")

install.packages("BiocManager")
BiocManager::install("rtracklayer")  # this takes ~15 minutes
BiocManager::install("GenomicFeatures")
```

Some of the data analysis was completed using plotly. To make this work, you will need to create a conda environment and link it to your R installation. In this conda environment (`r-reticulate` on my machine), install the following:

```bash
conda install plotly
pip install kaleido
```

## Downloading Data

Check the README.md files in each directory to see what to download for each folder.


## Getting Started

Setup:

```bash
cd path/to/escapegenecalculator

# Optional, since the output of these have already been generated for you
Rscript R/setup/calculate_exon_lengths.R  # ~1 minute
Rscript R/setup/extract_chromosome_lengths.R  # ~5 minutes
```

Old pipeline (output in `output-1`):

```bash
cd path/to/escapegenecalculator

# data/berletch-spleen
Rscript R/setup/generate_config.R -i data/berletch-spleen -o output-1
Rscript R/escapegenecalculator_old.R -i data/berletch-spleen
```

New pipeline with escape gene validation (output in `output-2`):

```bash
cd path/to/escapegenecalculator

# data/berletch-spleen
Rscript R/setup/generate_config.R -i data/berletch-spleen -o output-2  # default
Rscript R/escapegenecalculator.R -i data/berletch-spleen
```
```

## Contributing

Reach out to [harrison.c.wang@gmail.com](mailto:harrison.c.wang@gmail.com) and we can set up a time to discuss. You may also make your own branch or file a Github issue.
