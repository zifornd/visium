# visium

A quarto workflow to analyse and visualise 10x Visium Spatial Transcriptomics data 

## Contents

* [Overview](#overview)
* [Installation](#installation)
* [Usage](#usage)
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

Please ensure you also have `make` and [`R`](https://cran.r-project.org/bin/) version 4.0 or greater installed.

1. Install Renv

   Visium and all of its dependencies can be installed via the [renv](https://rstudio.github.io/renv/articles/renv.html) package manager. In each workflow `renv::restore()` will restore the project state from the `renv.lock` file. 

   For more details on `renv` please follow this [link](https://rstudio.github.io/renv/articles/collaborating.html).

   The first thing we install is `renv` which requires `yaml` to parse dependencies within Quarto Markdown files:

   ```console
   $ R

   R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
   Copyright (C) 2022 The R Foundation for Statistical Computing
   Platform: x86_64-pc-linux-gnu (64-bit)

   R is free software and comes with ABSOLUTELY NO WARRANTY.
   You are welcome to redistribute it under certain conditions.
   Type 'license()' or 'licence()' for distribution details.

   R is a collaborative project with many contributors.
   Type 'contributors()' for more information and
   'citation()' on how to cite R or R packages in publications.

   Type 'demo()' for some demos, 'help()' for on-line help, or
   'help.start()' for an HTML browser interface to help.
   Type 'q()' to quit R.

   > install.packages("yaml")
   > install.packages("renv")
   ```

2. Install quarto

   As this workflow is written as a number of quarto documents we next we need to install Quarto by following the [link](https://quarto.org/docs/get-started/): 

   We can then install the .deb file and test quarto install:
   
   ```console
   $ sudo dpkg -i quarto-0.9.640-linux-amd64.deb
   $ quarto -V # Test version and install
   ```
   See the following [link](https://quarto.org/docs/get-started/hello/text-editor.html) for more details on rendering quarto documents `.qmd`

## Usage

1. Create sample table

   ```console
   $ vim data/sample_table.csv # containing sample meta data 
   ```

   This table should look something like the following and must contain these headers:

   | sample | image | slide | group | area | index | files | protocol |
   |--------|-------|-------|-------|------|-------|-------|----------|
   | V1_Breast_Cancer_Block_A_Section_1 | V1_Breast_Cancer_Block_A_Section_1_image | V19L29-097 | slide1 | B1 | T1T2-F10 | V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5 | FF |
   | V1_Breast_Cancer_Block_A_Section_2 | V1_Breast_Cancer_Block_A_Section_2_image | V19L29-098 | slide2 | B1 | T1T2-H10 | V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5 | FF |

2. Fill in paramaters 

   These can be filled in at the top of each workbook e.g. for `01-data-loading`:

   ```console
        ---
        title: "Data loading"
        params:
            prefix:
            - "data/"
            marker:
            - "ACTA2"
            sample.sheet:
            - "V1_Breast_Cancer.csv"
        ---
   ```
   or supplied upon render using the `-P` flag:

3. Construct your make file

   ```console
   $ vim Makefile # Edit/Comment out .qmd files not in use.
   ```

4. Test make file

   ```console
   $ make -n
   ```

5. Run make file

   ```console
   $ make
   ```

Alternatively this quarto workflow may be rendered in separate parts using [RStudio](https://quarto.org/docs/get-started/hello/rstudio.html)  

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

- [Zifornd](https://www.zifornd.com/)

If you would like to be added to this list, please open a [pull request](https://github.com/zifornd/visium/pulls) with your information.

## Acknowledgements

This work was based primarily on the Seurat workflow for analysis of spatial transcriptomics data. Please see below for relevant citations:

> Hao Y, Hao S, Andersen-Nissen E, Mauck WM 3rd, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, Hoffman P, Stoeckius M, Papalexi E, Mimitou EP, Jain J, Srivastava A, Stuart T, Fleming LM, Yeung B, Rogers AJ, McElrath JM, Blish CA, Gottardo R, Smibert P, Satija R. Integrated analysis of multimodal single-cell data. Cell. 2021 Jun 24;184(13):3573-3587.e29. doi: 10.1016/j.cell.2021.04.048. Epub 2021 May 31. PMID: 34062119; PMCID: PMC8238499.

> Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM 3rd, Hao Y, Stoeckius M, Smibert P, Satija R. Comprehensive Integration of Single-Cell Data. Cell. 2019 Jun 13;177(7):1888-1902.e21. doi: 10.1016/j.cell.2019.05.031. Epub 2019 Jun 6. PMID: 31178118; PMCID: PMC6687398.

> Hafemeister, C., Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol 20, 296 (2019). https://doi.org/10.1186/s13059-019-1874-1

> Amezquita, R.A., Lun, A.T.L., Becht, E. et al. Orchestrating single-cell analysis with Bioconductor. Nat Methods 17, 137???145 (2020). https://doi.org/10.1038/s41592-019-0654-x

The workflow was motivated by the following projects:

- [Seurat](https://satijalab.org/seurat/articles/spatial_vignette.html)
- [OSCA](https://github.com/OSCA-source/OSCA)
- [10XGenomics/spaceranger](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger)

The documentation was informed by the following articles:

- [easiest way to create a readme](https://readme.so)
- [writing a friendly readme](https://rowanmanning.com/posts/writing-a-friendly-readme/)
- [writing well for the web](https://www.gov.uk/guidance/content-design/writing-for-gov-uk)

## License

Visium is licensed under the [MIT](LICENSE) license.  
Copyright &copy; 2022, Zifo rnd
