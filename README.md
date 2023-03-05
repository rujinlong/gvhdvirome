
# GVHD virome analyses

<!-- badges: start -->
<!-- badges: end -->

This repository contains the R package `gvhdvirome`, which contains the data and code for reproducible virome analyses in the paper "..." by ...

## Installation

You can install the development version of `gvhdvirome` with:

``` r
# install.packages(c()"devtools", "renv"))
devtools::install_github("rujinlong/gvhdvirome")
renv::restore()
```

Raw data are not included in the package. To download the raw data, run the following code:

``` bash
wget https://www.dropbox.com/xxxxxx/gvhdvirome_data.tar.gz   # TODO: update link
tar -xvzf gvhdvirome_data.tar.gz
mkdir -p ${PRJDIR}/workflow/data/00-raw
mv gvhdvirome_data ${PRJDIR}/workflow/data/00-raw
```

## Reproducible analyses

The analyses in the paper can be reproduced by running code in the `workflow` directory. Each file in the `workflow` directory is a [Quarto](https://quarto.org/) document that contains the code and results for a specific analysis. The code in each file can be run in RStudio by clicking the "Run Document" button in the top right corner of the editor pane. The results of each analysis are saved in the `data/<analysis_name>` directory, where `<analysis_name>` is the file name of the analysis.