# G*enome-based* In*cidence Estimation* Pipe*line*

This pipeline infers the trajectory of an effective population size (or incidence) for a viral pandemic from a collection of time-stamped viral sequences. The pipeline has, so far been tested for SARS-CoV-2.
In brief: Viral sequence data is placed into redundant temporal bins. For each bin, a parameter is inferred that correlates with the effective population size estimate (or incidence) of the infection. GInPipe then smoothes over all derived parameters and reconstructs continuous trajectory of the effective population size estimate (or incidence).
  -   [Operating systems and dependencies](#operating-systems-and-dependencies)
  -   [Input](#input)
  -   [Running the pipeline](#running-the-pipeline)
      -   [1. Prerequisites](#1-prerequisites)
      -   [2. Initialization](#2-initialization)
      -   [3. Command execution](#3-command-execution)
  -   [Output](#output)
  -   [Demo](#demo)
  -   [Running the pipeline for SARS-CoV-2 on GISAID data](running-the-pipeline-for-sars-cov-2-on-gisaid-data)

## Operating systems and dependencies

This workflow was tested on macOS Mojave Version 10.14.4 and macOS Catalina 10.15.7.

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
  - matplotlib
  - minimap2
  - ggplot2 
  - mgcv
  - MASS
  - R0
  - scales
```

They are installed automatically upon execution using the environment file [`env.yml`](./env/env.yml) and R scripts [`computeInterpolation.R`](./scripts/RScripts/splines/computeInterpolation.R) and [`computeR0.R`](./scripts/RScripts/splines/computeR0.R)

Preparing the environment takes about 3 minutes.

## Input
As an input the pipeline requires a file containing sequences and a file with a reference consensus sequence.

For the sequences it is important that they contain a sequencing-, or better, sample-date. The date must have the format **%YYYY-%mm-%dd**
and can be either part of the sequence-name or provided in an additional file.
- If the date is part of the sequence-name, then the name should look like this: **'some_name | %YYYY-%mm-%dd'**.   
- If the date is provided in an additional file, add the date to corresponding FASTA headers.

## Running the pipeline

This is a small guide on how to set up and run the pipeline.

### 1. Prerequisites
Some tools have to be installed to run the analysis. We recommend following the steps below to set up the pipeline.

#### 1.1 Install Conda/Miniconda

Conda will manage the dependencies of our pipeline. Instructions can be found here: https://docs.conda.io/projects/conda/en/latest/user-guide/install


#### 1.2 Create the working environment

Create a new environment where the pipeline will be executed, for example like this:

```
conda create --name GInPipe
```

Then to activate this environment, type:

```
conda activate GInPipe
```

#### 1.3 Install Snakemake

Snakemake is the workflow management system we use. Install it like this:

```
conda install snakemake
```

### 2. Initialization

As input the pipeline requires names of the sequence file, the reference genome, and binning parameters.
These variables are stored in [`config.yaml`](./config.yaml) and used as wildcards to create and link files with each other or as parameters for the binning. For more information about the YAML markup format refer to documentation: https://yaml.org

The specified paths in the config file should either be absolute, or relative to the work environment specified with -d in the snakemake call (see 3.0).

#### 2.1 Raw sequences
The pipeline requires a file containing sequences, with the date in the sequence-name in GISAID format (date in format "%Y-%m-%d" at the end of header after a vertical bar).

For sequence files containing the date within the sequence-name, copy the file path and paste it into the variable **samples** of [`config.yaml`](./config.yaml):

  ```
  samples: "path/to/sequences/data"
  ```
If the headers in the sequence file do not contain the date, you can add it to headers using a custom script, given a meta table is provided along with the FASTA file.

#### 2.2 Reported cases data file

To compare estimated population dynamics with reported active cases, you can include a table in the folder [`reported_cases`](./reported_cases). Also provide the following parameters in the corresponding config field like this:

  ```
  reported_cases: ["path/to/reported_cases.csv","\t","date","active_cases","%m/%d/%y"]
  ```

where the first element of the list is the file name with format extension, the second element is the delimiter type in this file, date column name, active cases column name, and a format the date is stored in.

If no reported cases data is provided, leave the fields empty like this:

  ```
  reported_cases: []
  ```

#### 2.3 Reference consensus sequence
Copy and paste the file path of reference/consensus sequence into the variable **consensus** of [`config.yaml`](./config.yaml). 

  ```
  consensus: "path/to/consensus/sequence.fasta"
  ```

#### 2.4 Binning parameters
You also have to set the parameters for some of the binning methods in [`config.yaml`](./config.yaml).
You can set the number of sequences per bin, and the number of days.
Parameters can be given as an array. Additionally, minimal bin size and maximal days span should be
provided.

```
number_per_bin: [20, 30]
days_per_bin: [7, 10, 30]
min_bin_size: 15
max_days_span: 21
```

If parameter **number_per_bin** is an empty list, a default mode with predefined fractions of sequences (2%, 5%, 7%) is used. Alternatively, all arrays can be given in the configuration file as a list, like this:

```
number_per_bin:
    - 20
    - 30
```

#### 2.5 Effective reproduction number prediction

The workflow can calculate and plot the prediction of effective reproduction number. If this prediction is desired, specify it via a Boolean variable in the configuration file, for example like this:

```
R0: y
```

If no prediction is wanted, specify it in the configuration file like this:

```
R0: n
```

Other options for specifying this parameter also work. For examples see https://yaml.org/type/bool.html

### 3. Command execution

To run the pipeline, go to the pipeline directory (where the Snakefile is) and activate the conda environment created in step 1.2. Then enter the following command to execute the pipeline:

```
snakemake --use-conda --snakefile GInPipe --configfile path/to/config.yaml -j -d path/to/workdir
```

The ---use-conda parameter allows Snakemake to install packages that are listed in the environment file [`env.yml`](./env/env.yml). With parameter --configfile you can give the configuration file [`config.yml`], described above. The -j parameter determines the number of available CPU cores to use in the pipeline. Optionally you can provide the number of cores, e.g. -j 4. With parameter -d you can set the work directory, i.e. where the results of the pipeline are written to.

## Output
The pipeline creates a folder **'results'**, containing all (intermediate) outputs, with the following structure:
```
    ├── results                                 # Main results folder
    │   ├── bam                                 # sorted and indexed bam files
    │   ├── bins                                # binning results
    │       ├── cal_week                        # binned by calendar week
    │       ├── eq_days_10                      # binned by equal days                       
    │       ├── eq_size_100                     # binned by equal number of sequences
    │       └── fuzzy_days_100                  # binned by equal number of sequences (fuzzy)
    │               ├── bin_*.bam               # binned sequences as BAM
    │               ├── bin_*.bai               # index files                       
    │               ├── header_*.tsv            # header files (seq. name & date)
    |               ├── range_*.tsv             # range of dates of the corresponding bin
    |               └── list_of_files.tsv       # list of file names in the binning mode
    │   ├── bins_incidence                      # Individual binning results plots (and tables)
    │   ├── meta                                # Meta information about all used sequences (name and collection date)
    │   ├── incidence                           # Plots and tables for final interpolated trajectory
    |       ├── incidence.csv                   # table with interpolated population size estimates
    │       ├── phi_estimate.csv                # table with raw population size estimates                       
    │       ├── incidence.pdf                   # plot of interpolated population size
    |       └── wdots_incidence.pdf             # plot of interpolated population size with dot size scaled by bin size
    │   ├── r0                                  # Reproduction number estimate
    |       ├── r0.csv                          # table with daily reproduction number estimates
    │       └── r0.pdf                          # plot of daily reproduction number estimates                       
    │   └── raw                                 # Preprocessed files
    │   
    └── ...
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
snakemake --use-conda --snakefile GInPipe --configfile demo/demo_config.yaml -j -d demo
```

It may take about 3 minutes to prepare the environment and around 2 minutes to run the pipeline.
The result folder is created in the [`demo`](./demo) folder where you find the output files, as described above.

## Running the pipeline for SARS-CoV-2 on GISAID data

To run the pipeline on real data, such as COVID sequences from GISAID, download the accordings sequence data as well as the reference sequence and adapt the paths in the config.yaml as explained in 2.1.and 2.3. 
