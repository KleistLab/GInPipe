# G*enome-based* In*cidence Estimation* Pipe*line*

[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
![](https://img.shields.io/github/v/release/kleistlab/ginpipe)
[![DOI:10.1101/2021.05.14.21257234](http://img.shields.io/badge/DOI-10.1101/2021.05.14.21257234-blue.svg)](https://doi.org/10.1101/2021.05.14.21257234)
[![DOI](https://zenodo.org/badge/359502958.svg)](https://zenodo.org/badge/latestdoi/359502958)

This pipeline infers the trajectory of an effective population size (or incidence) for a viral pandemic from a collection of time-stamped viral sequences. The pipeline has, so far been tested for SARS-CoV-2.
In brief: Viral sequence data is placed into redundant temporal bins. For each bin, a parameter is inferred that correlates with the effective population size estimate (or incidence) of the infection. GInPipe then smoothes over all derived parameters and reconstructs continuous trajectory of the effective population size estimate (or incidence) [[1]](#1).

  -   [System requirements](#system-requirements)
      -   [Operating systems](#operating-systems)
      -   [Prerequisites](#prerequisites)
      -   [Dependencies](#dependencies)
  -   [Input](#input)
  -   [Running the pipeline](#running-the-pipeline)
      -   [Initialization](#initialization)
      -   [Execution](#execution)
  -   [Output](#output)
  -   [Demo](#demo)
  -   [Running the pipeline for SARS-CoV-2 on GISAID data](#running-the-pipeline-for-sars-cov-2-on-gisaid-data)

## System requirements

### Operating systems
This workflow was tested on macOS Mojave Version 10.14.4, macOS Catalina Version 10.15.7, as well as Linux Version 5.0.0-38.

### Prerequisites
Some tools have to be installed to run the analysis. We recommend following the steps below to set up the pipeline.

#### Install Conda/Miniconda

Conda will manage the dependencies of our pipeline. Instructions can be found here: https://docs.conda.io/projects/conda/en/latest/user-guide/install.


#### Create the working environment

Create a new environment from the given environment config in [`env.yml`](./env/env.yml), where the pipeline will be executed.
Go to the main folder of the repository and call:

```
conda env create -f env/env.yml
```

This step may take a few minutes.

To activate the environment type:

```
conda activate GInPipe
```

#### Install Snakemake

Snakemake is the workflow management system we use. Install it in your activated environment like this:

```
conda install -c conda-forge snakemake
```

NOTE: In case conda is not able to find the packages for snakemake (which was the case for the Linux version), you can install mamba in your environment

```
conda install -c conda-forge mamba
```

and download snakemake with

```
mamba install -c conda-forge -c bioconda snakemake
```

Detailed Snakemake installation instruction using mamba can be found here: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html.

#### Install R

To run R routines, R including Rscript needs to be installed for the workflow. If it is not yet, you can install it together with the needed packages in your activated environment with conda or mamba

```
mamba install -c conda-forge -c bioconda r-base r-ggplot2 r-r0 r-scales
```


### Dependencies

This workflow uses the following dependencies:

```
  - bbmap
  - seqkit
  - samtools
  - numpy
  - pysam
  - biopython
  - pandas
  - scipy
  - pyvcf
  - minimap2
  - ggplot2
  - R0
  - scales
  - pip
  - devtools
  - ginpipepy
  - ginpiper
```

They are installed automatically upon execution using the environment file [`env.yml`](./env/env.yml) and R scripts [`computeInterpolation.R`](./scripts/RScripts/splines/computeInterpolation.R) and [`computeR0.R`](./scripts/RScripts/splines/computeR0.R). *ginpipepy* is a custom package that contains functions for binning, population size calculations, and masking. It is installed within the environment using *pip* from https://github.com/KleistLab/ginpipepy. *ginpiper* is a package that contains functions for smoothing, plotting, and date transformations. It is another custom package used by the pipeline and is installed using *devtools* from https://github.com/KleistLab/ginpiper.

## Input
As an input, the pipeline requires a file containing sequences and a file with a reference consensus sequence.

For the sequences it is important that they contain a sequencing-, or better, sampling-date. The date must have the format **%YYYY-%mm-%dd**
and must be part of the sequence-name of the FASTA header, which should look like this: **'>some_name|%YYYY-%mm-%dd'**.

The header of the reference sequence can contain an arbitrary name, but importantly without white spaces:
**'>some_reference_name'**. If there are whitespaces present in the reference name string, a substring before the first whitespace will be used.

## Running the pipeline

This is a small guide on how to set up and run the pipeline.


### Initialization

As input the pipeline requires names of the sequence file, the reference genome, and binning parameters.
These variables are stored in [`config.yaml`](./config.yaml) and used as wildcards to create and link files with each other or as parameters for the binning. For more information about the YAML markup format refer to documentation: https://yaml.org

The specified paths in the config file should either be absolute, or relative to the work environment specified with -d in the snakemake call (see below in [Execution](#execution)).

**Note:** paths to files in [`config.yaml`](./config.yaml) should be either:
- absolute,
- or relative to the main pipeline Snakefile,
- or relative to the working directory (see below in Execution section).

#### 1 Sample sequences
The pipeline requires a file containing sequences, with the date in the sequence-name in GISAID format (date in format "%Y-%m-%d" at the end of header after a vertical bar).

For sequence files containing the date within the sequence-name, copy the file path and paste it into the variable **samples** of [`config.yaml`](./config.yaml):

  ```
  samples: "path/to/sequences/data"
  ```
If the headers in the sequence file do not contain the date, you can add it to headers using a custom script, given a meta table is provided along with the FASTA file.

In the field **group**, you can provide a name for the given samples, for example the country, location or any specifier for the given sequences, which is used for the plots.

  ```
  group: country
  ```


#### 2 Reported cases data file

To compare estimated population dynamics with reported active cases, provide the following parameters in the corresponding config field like this:

  ```
  reported_cases: ["path/to/reported_cases.csv","\t","date","new_cases","%m/%d/%y"]
  ```

where the first element of the list is the file name with format extension, the second element is the delimiter type in this file, date column name, active cases column name, and a format the date is stored in.

If no reported cases data is provided, leave the fields empty like this:

  ```
  reported_cases: []
  ```

#### 3 Reference sequence
Add the file path of reference/consensus sequence into the variable **consensus** of [`config.yaml`](./config.yaml).

  ```
  consensus: "path/to/consensus/sequence.fasta"
  ```

#### 4 Binning parameters
You also have to set the parameters for some of the binning methods in [`config.yaml`](./config.yaml).
You can set the number of sequences per bin, and the number of days.
Parameters can be given as an array. Additionally, minimal bin size and maximal days span should be
provided.

  ```
  number_per_bin: [20, 30]
  days_per_bin: [7, 10, 30]
  min_bin_size: 15
  max_days_span: 21
  min_days_span: 2
  ```

If parameter **number_per_bin** is an empty list, a default mode with predefined fractions of sequences (2%, 5%, 7%) is used. Alternatively, all arrays can be given in the configuration file as a list, like this:

  ```
  number_per_bin:
      - 20
      - 30
  ```

#### 5 Threshold for mutations

Low-abundance point mutations are filtered out, to avoid the consideration of sequence errors.
The cutoff for the amount of mutations at a certain sequence position in the whole dataset is set with:

  ```
  freq_cutoff: 2
  ```

#### 6 Effective reproduction number prediction

The workflow can calculate and plot the prediction of effective reproduction number. If this prediction is desired, specify it via a boolean variable in the configuration file, for example like this:

  ```
  R0: yes
  ```

If no prediction is wanted, specify it in the configuration file like this:

  ```
  R0: no
  ```

Alternative options for specifying this parameter are "on/off" and "true/false" (case insensitive).
<!--Other options for specifying this parameter also work. For examples see https://yaml.org/type/bool.html Single characters do not work. -->

#### 7 Masking using a Variant Calling File

A VCF containing sites that need to be masked can also be provided via config:

  ```
  masking: "path/to/vcf"
  ```

The pipeline will check for flag "mask" in the FILTER column to apply masking to a particular site.

If no masking file is provided, leave the field empty:

  ```
  masking: ""
  ```

### Execution

To run the pipeline activate (if not activated yet) the conda environment with

```
conda activate GInPipe
```

Go to the pipeline directory (where the Snakefile named *GInPipe* is located) and enter the following command to execute the pipeline:

```
snakemake --snakefile GInPipe --configfile path/to/config.yaml -j -d path/to/workdir
```

With parameter --configfile you can give the configuration file, described above. The -j parameter determines the number of available CPU cores to use in the pipeline. Optionally you can provide the number of cores, e.g. -j 4. With parameter -d you can set the work directory, i.e. where the results of the pipeline are written to.

## Output
The pipeline creates a folder **'results'**, containing all (intermediate) outputs, with the following structure:
```
    ????????? results                                 # Main results folder
    ???   ????????? bam                                 # sorted and indexed bam files
    ???   ????????? bins                                # binning results
    ???       ????????? cal_week                        # binned by calendar week
    ???       ????????? eq_days_10                      # binned by equal days                       
    ???       ????????? eq_size_100                     # binned by equal number of sequences
    ???       ????????? fuzzy_days_100                  # binned by equal number of sequences (fuzzy)
    ???               ????????? bin_*.bam               # binned sequences as BAM
    ???               ????????? bin_*.bai               # index files                       
    ???               ????????? header_*.tsv            # header files (seq. name & date)
    |               ????????? range_*.tsv             # range of dates of the corresponding bin
    |               ????????? list_of_files.tsv       # list of file names in the binning mode
    ???   ????????? bins_incidence                      # Individual binning results tables
    ???   ????????? meta                                # Meta information about all used sequences (name and collection date)
    ???   ????????? incidence                           # Plots and tables for final interpolated trajectory
    |       ????????? incidence.csv                   # table with interpolated population size estimates
    ???       ????????? phi_estimate.csv                # table with raw population size estimates                       
    ???       ????????? incidence.pdf                   # plot of interpolated population size
    |       ????????? wdots_incidence.pdf             # plot of interpolated population size with dot size scaled by bin size
    ???   ????????? r0                                  # Reproduction number estimate
    |       ????????? r0.csv                          # table with daily reproduction number estimates
    ???       ????????? r0.pdf                          # plot of daily reproduction number estimates                       
    ???   ????????? raw                                 # Preprocessed files
    ???   
    ????????? ...
```
The analysis results for different binning modes (plots and tables) can be found in the subfolder **results/bins_incidence** and the smoothed final trajectory can be found in subfolder **results/incidence**. The final output (depending on the setup of the pipeline) contains the following tables and plots:

- In folder *incidence*:
    - *phi_estimate.csv*: table containing estimates of population size for all binning strategies
    - *incidence.csv*: table containing interpolated final trajectory of population size; the trajectory is calculated by  combining all binning strategies
    - *incidence.pdf*: plot of interpolated trajectory, optionally overlayed with reported cases trajectory (if the corresponding table was given)
    - *wdots_incidence.pdf*: plot of interpolated trajectory with point estimate dots scaled by corresponding sub-sample size, optionally overlayed with reported cases trajectory (if the corresponding table was given)
- In folder *r0* (if option to calculate reproduction number is chosen):
    -  *r0.csv*: table containing daily reproductive number estimates; calculated from the interpolated trajectory
    -  *r0.pdf*: plot of daily reproductive number estimates with confidence interval

## Demo

A demo sequence set and a reference sequence are included in repository folders [`demo`](./demo).
The directory contains a simulated data set with

- the reference sequence (demo_reference.fasta)
- a fasta file containing the newly emerging sequences over time (demo_samples.fasta)
- the underlying true number of emerging sequences (demo_reported_cases.tsv)
- the config file to call the pipeline with (demo_config.yaml)


To run the pipeline go into the repository where the GInPipe file is located and run

```
snakemake --snakefile GInPipe --configfile demo/demo_config.yaml -j -d demo
```

It may take around 2 minutes to run the pipeline.
The result folder is created in the [`demo`](./demo) folder where you find the output files, as described above. The incidence plot of the demo sample should look like this:

![alt text](https://github.com/KleistLab/GInPipe/blob/main/demo/demo_result.png)

## Running the pipeline for SARS-CoV-2 on GISAID data

To run the pipeline on real data, such as COVID sequences from GISAID, download the according sequence data as well as the reference sequence and adapt the paths in the config.yaml as explained in [Initialization](#initialization).

## Reference
<a id="1">[1]</a>
Smith, M. R. and Trofimova, M., et al. (2021). Rapid incidence estimation from SARS-CoV-2 genomes reveals decreased case detection in Europe during summer 2020. Nature Communications  12, 6009. https://doi.org/10.1038/s41467-021-26267-y
