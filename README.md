# nanopath

![](https://img.shields.io/badge/lang-python-41ab5d.svg)
![](https://img.shields.io/badge/version-0.1.1-addd8e.svg)
![](https://img.shields.io/badge/biorxiv-v0-f7fcb9.svg)

Python package to manage pipelines and data processing for the `NanoPath` submodules and command-line interface :snake:

## Overview

**`v0.1.2: beastling`**

## Install

`pip install git+https://github.com/np-core/nanopath`

### `ILLUMINA` bacterial phylogenomics, exploratory phylodynamics as `Nextflow`

TBP


### `BEAST` analyse templates for bacterial phylodynamic analyses

:whale: Containers:

`docker pull esteinig/phybeast:dev`, `singularity pull docker://esteinig/phybeast:dev`

What the module does:

`BEAST2` templates with GPU accerlaration on `BEAGLE-v3.2.1` improves compute time of models considerably. We routinely reduce compute rates from hours to minutes per million MCCMC samples on datasets with ~600 isolates represented by a non-recombinant core SNP alignment from **Snippy** and **Gubbins**. Preprocessing is implemented in a Nextflow pipeline. Updates to the `nanopath phybeast` command-line interface include tasks for creating `XML` files with sensible priors and clock configuration in `YAML` format to facilitate use of the **Coalescent Bayesian Skyline** (`np phybeast beastling xml-cosky --help`) and the **Birth-Death Skyline Serial** (`np phybeast beastling xml-bdss --help`) models. These allow us to estimate important demographic parameters of the pathogen genome pool, and with sufficient sampling resolution, can provide estimates across lineage emergence events and distinct outbreaks. In the next update we will implement the structured `Multi-Type Birth-Death` (`np phybeast beastling xml-bdmm --help`) model as well.

Features:

- `np phybeast beastling xml-cosky`
- `np phybeast beastling xml-bdss`

How to:



Publications:



### `NANOPORE` hybrid phylogenetic lineage and outbreak attribution with `Medaka`

TBP

### `NANOPORE` polishing of pure and hybrid `Illumina` genomes as `Nextflow`

TBP

### Dependencies
