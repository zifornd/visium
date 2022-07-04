# visium

A quarto workflow to analyse and visualise 10x Visium Spatial Transcriptomics data 

## Contents

* [Overview](#overview)
* [Installation](#installation)
* [Usage](#usage)
* [Documentation](#documentation)
* [Contributing](#contributing)
* [Authors](#authors)
* [Tests](#tests)
* [Acknowledgements](#acknowledgements)
* [License](#license)

## Overview

Visium is a quarto workflow to process and analyse 10x Visium Spatial transcriptomic data from the 10x Genomics platform. It is compatible with both FFPE and Fresh frozen protocols and utilises the following key analysis frameworks:

* [Seurat](https://satijalab.org/seurat/articles/spatial_vignette.html)
* [OSCA](https://github.com/OSCA-source/OSCA)

## Installation

Visium and all of its dependencies can be installed via the [mamba](https://github.com/mamba-org/mamba) package manager:

1. Install Dependencies

   ```console
   $ 
   ```

## Usage

1. Run make file

   ```console
   $ make
   ```

*For more information, see the [Usage](documentation.md#usage) section of the documentation.*

## Documentation

Full documentation for Visium is available [here](documentation.md)

## Support

If you need help, open an [issue](https://github.com/zifornd/visium/issues) with one of the following labels:

- help wanted (extra attention is needed)
- question (further information is requested)

## Feedback

If you have any suggestions, open an [issue](https://github.com/zifornd/visium/issues) with one of the following labels:

- documentation (improvements or additions to documentation)
- enhancement (new feature or request)

## Contributing

To contribute to Visium, clone this repository locally and commit your code on a separate branch. Please generate unit tests for your code and run a linter before opening a pull request.

You can find more details in the [Contributing](https://github.com/zifornd/.github/blob/main/CONTRIBUTING.md) guide. 

Participation in this project is subject to a [Code of Conduct](https://github.com/zifornd/.github/blob/main/CODE_OF_CONDUCT.md).

## Authors

Visium was developed by [Ben Southgate](https://github.com/bensouthgate) and [James Ashmore](https://www.github.com/jma1991).

If you would like to be added to this list, please open a [pull request](https://github.com/zifornd/visium/pulls) with your contribution.

## Citation


## Used By

Visium is used by the following companies and institutes:

- []()

If you would like to be added to this list, please open a [pull request](https://github.com/zifornd/visium/pulls) with your information.

## Acknowledgements

This work was based primarily on the Seurat workflow for analysis of spatial transcriptomics data. Please see below for relevant citations:

> Hao Y, Hao S, Andersen-Nissen E, Mauck WM 3rd, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, Hoffman P, Stoeckius M, Papalexi E, Mimitou EP, Jain J, Srivastava A, Stuart T, Fleming LM, Yeung B, Rogers AJ, McElrath JM, Blish CA, Gottardo R, Smibert P, Satija R. Integrated analysis of multimodal single-cell data. Cell. 2021 Jun 24;184(13):3573-3587.e29. doi: 10.1016/j.cell.2021.04.048. Epub 2021 May 31. PMID: 34062119; PMCID: PMC8238499.

> Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM 3rd, Hao Y, Stoeckius M, Smibert P, Satija R. Comprehensive Integration of Single-Cell Data. Cell. 2019 Jun 13;177(7):1888-1902.e21. doi: 10.1016/j.cell.2019.05.031. Epub 2019 Jun 6. PMID: 31178118; PMCID: PMC6687398.

> Hafemeister, C., Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol 20, 296 (2019). https://doi.org/10.1186/s13059-019-1874-1

The workflow was motivated by the following projects:

- [Seurat](https://satijalab.org/seurat/articles/spatial_vignette.html)
- [OSCA](https://github.com/OSCA-source/OSCA)
- [10XGenomics/spaceranger](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger)

The documentation was informed by the following articles:

- [easiest way to create a readme](https://readme.so)
- [writing a friendly readme](https://rowanmanning.com/posts/writing-a-friendly-readme/)
- [writing well for the web](https://www.gov.uk/guidance/content-design/writing-for-gov-uk)

## License

Visium is licensed under the [MIT](LICENSE.md) license.  
Copyright &copy; 2020, Zifo rnd
