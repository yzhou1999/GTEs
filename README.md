# GTE (Group Technical Effects)

`GTE` quantifies batch effects for individual genes in single-cell data. If you have any suggestions or problems on this package, please contact Y. Zhou (yangz@stu.hit.edu.cn), or submit through this GitHub page directly.

# Installation
To run `GTE`, install from GitHub through ``devtools`` directly:
```R
install.packages('devtools')
library(devtools)
devtools::install_github("yzhou1999/GTE")
```

# Usage

For usage examples, please see the `vignettes` directory of the repository.

* [Group technical effects calculation and highly technical genes selection](https://yzhou1999.github.io/GTE/articles/GTE_usage.html)


The codes and source data to benchmark the batch separation metrics please see [here](https://github.com/yzhou1999/GTE/tree/master/source_data).

# Datasets used in the manuscript
The datasets are available freely on public databases, and can also be downloaded in the [Zenodo repository](https://doi.org/10.5281/zenodo.13358933).


# Dependencies
GTE has been successfully installed and used on Windows, Linux and Mac OS (R version >= 4.0.2).
