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

#### Install Mamba/Conda

Mamba and Conda are a package manager that help to manage the package dependencies of our pipeline in environments. Mamba is a faster re-implementation of conda, hence we use it in the following.

Follow the instructions to install Mamba: 
https://mamba.readthedocs.io/en/latest/mamba-installation.html


Alternatively, you can also use Conda. Instructions can be found here: 
https://docs.conda.io/projects/conda/en/latest/user-guide/install


#### Create the working environment

Create a new environment from the given environment config in [`env.yml`](./env/env.yml), where the pipeline will be executed.
Go to the main folder of the repository and call:

```
mamba env create -f env/env.yml
```

This step may take a few minutes.

To activate the environment type:

```
mamba activate GInPipe3
```
---
*Note:* In case the conda environments are not found anymore, you may need to add the location to the environment directories.

```
conda config --append envs_dirs /path/to/conda/envs/
```


### Dependencies

This workflow uses the following dependencies:


  | python      | 3.9.18  |
  | snakemake   | 7.32.3  |
  | biopython   | 1.78    |
  | pandas      | 2.0.3   |
  | scipy       | 1.11.1  |
  | bbmap       | 38.18   |
  | numpy       | 1.24.4  |
  | matplotlib  | 3.7.2   |
  | scikit-fda  | 0.8.1   |
  | pysam       | 0.21.0  |
  | seqkit      | 2.4.0   |
  | samtools    | 1.17    |
  | minimap2    | 2.26    |



They are installed automatically upon creating the environment using the environment file [`env.yml`](./env/env.yml)

---
*Note:* In case some of the packages or tools can not automatically be installed, delete or comment out the dependency in the yml file, create the environment without it, activate the environment and install it manually.


## Running the pipeline

This is a small guide on how to set up and run the pipeline.


### Initialization

As input the pipeline requires names of the input files, binning- and smoothing parameters.
These variables are stored in [`config.yaml`](./config.yaml).

The specified paths in the config file should either be absolute, or relative to the work environment specified with -d in the snakemake call (see below in [Execution](#execution)).

*Note:* paths to files in [`config.yaml`](./config.yaml) should be either:
- absolute,
- or relative to the main pipeline Snakefile,
- or relative to the working directory (see below in Execution section).


#### Input files

As an input, the pipeline requires the samples along with a sequecing-, or better, collection date (format **%YYYY-%mm-%dd**). The file can be either given as a fasta file or a table containing the dna profile (SNVs).
The pipeline automatically chooses the workflow according to the file extension. So please name your sample file either with *.fasta* or *.csv*.


Copy the file path and paste it into the variable **samples** of [`config.yaml`](./config.yaml):

  ```
  samples: "path/to/sequences/datafile"
  ```

In the field **name**, you can provide a name for the given samples, for example the country, location or any specifier for the given sequences, which is used for the plots.

  ```
  name: country
  ```

**1. FASTA file**

If the input file is a FASTA file with the sequences,the date must be part of the header behind a vertical bar, similar to sequence-names in GISAID: **'>some_name|%YYYY-%mm-%dd'**.


Along with the FASTA file a reference sequence needs to be provided (also in FASTA format). 
The header of the reference sequence can contain an arbitrary name, but importantly without white spaces:
**'>some_reference_name'**. If there are whitespaces present in the reference name string, a substring before the first whitespace will be used.


Add the file path of reference sequence into the variable **reference** of [`config.yaml`](./config.yaml).

  ```
  reference: "path/to/consensus/sequence.fasta"
  ```


**2. CSV file**

The input file can also be a CSV file containing one column with the date and one column with the mutations which are separated by blank. The format of the each mutation is WtPositionMut. The column names must be given as *date* and *dna_profile*

```
date,dna_profile
2023-01-06,C241T T595C T670G C1931A C2790T
2023-01-06,C44T T670G T2954C
```

Here, no reference file needs to be provided.


#### Binning parameters

The pipeline assigns the sequences into consecutive bins. Parameters to create and filter those bins can be set in the [`config.yaml`](./config.yaml).
For the binning strategy you can set the number of sequences per bin, and the number of days.
Parameters can be given as an array.
  ```
  seq_per_bin: [20, 30]
  days_per_bin: [7, 10, 30]
  
  ```

Alternatively, all arrays can be given in the configuration file as a list, like this:

  ```
  number_per_bin:
      - 20
      - 30
  ```

Optionally, you can restrict which bins should be considered for the phi estimatation, by setting the minimal bin size (default 1), as well as the  minimal and maximal days spanning the bin (default 1 and 21).


  ```
  min_bin_size: 200
  min_days_span: 2
  max_days_span: 21
  ```

#### Parameters for phi estimation

Low-abundance point mutations can be filtered out (optionally), to avoid the consideration of sequence errors.
The cutoff for the amount of mutations at a certain sequence position in the whole dataset is set with:

  ```
  freq_cutoff: 2
  ```

If only a certain region within the sequence is of intereset or regions are prone to errors, positions can be optionally masked by setting

  ```
  masking: 1-10,15,200-210
  ```

The mutations are separated by ",". Positions ranges are given with pos1-pos2.

If no positions are masked, leave the field empty:


  ```
  masking: ""
  ```


A line is smoothed through the phi point estimates with a kernel smoother. The smoothing bandwidth is set with


  ```
  smoothing_bandwidth_phi: 3
  ```

The parameter is optional, with default 7.


#### Estimate the minimal number of infected 

The smoothed phi estimates are a proxy for the underlying true number of infected. 
If the reported new cases are available for the same time horizon, we can estimate the minimal number of truely infected. 

The reported cases table can be provided with the following parameters in the [`config.yaml`](./config.yaml):

  ```
  reported_cases: ["path/to/reported_cases.csv",",","date","new_cases","%m/%d/%y"]
  ```

where the first element of the list is the file name with format extension, the second element is the delimiter type in this file, date column name, active cases column name, and a format the date is stored in.

If no reported cases data is provided, leave the fields empty like this:

  ```
  reported_cases: []
  ```

Before the minimal true incidence is calculated the phi estimates and reported cases can be smoothed to prevent the normalisation by an outlier. The smoothing bandwith is set with

```
smoothing_bandwidth_mi: 7
```

The default value is 7.

Also the time frame to be considered can be optionally set with parameters

```
from_date: "2022-01-01"
to_date:"2022-12-31"
```

If all dates are considered, leave the fields empty: 

```
from_date: ""
to_date:""
```

### Execution

To run the pipeline activate (if not activated yet) the conda environment with

```
conda activate GInPipe3
```

Go to the pipeline directory (where the Snakefile named *GInPipe* is located) and enter the following command to execute the pipeline:

```
snakemake --snakefile GInPipe --configfile path/to/config.yaml -j -d path/to/workdir
```

With parameter --configfile you can give the configuration file, described above. The -j parameter determines the number of available CPU cores to use in the pipeline. Optionally you can provide the number of cores, e.g. -j 4. With parameter -d you can set the work directory, i.e. where the results of the pipeline are written to.

## Output
The pipeline creates a folder **'results'**, containing all outputs, with the following structure:
```
    ├── results                                   # Main results folder
    │   ├── phi_estimates                         # phi estimation results
    │       ├── phi_estimates_per_bin_*.csv       # binning results tables containing the point estimates per bin
    │       ├── smoothed_phi_estimates_*.csv      # table containing the smoothed phi value per day             
    │       ├── sequence_stats_per_day_*.csv      # table with sequence statistics per day
    │       └── plot_smoothed_phi_estimates_*.pdf # plot with point estimates and smoothed line
    │   ├── incidence                             # minimal incidence results
            ├── minimal_incidence_*.csv           # table with minimal incidence
    │       └── plot_minimal_incidence_*.pdf      # plot with minimal incidence and reported cases
```

## Demo

As a demo, we provide a set of German SARS-CoV-2 sequences from 2022[^1] along with reported cases[^2] in the folder [`demo`](./demo).

The directory comprises a data set with

- a csv file containing the dna profiles of SARS-CoV-2 sequences over time (demo_samples.csv)
- the reported cases (demo_reported_cases.csv)
- the config file to call the pipeline with (demo_config.yaml)


To run the pipeline go into the repository where the GInPipe file is located and run

```
snakemake --snakefile GInPipe --configfile demo/demo_config.yaml -j -d demo
```

or 

```
snakemake --snakefile GInPipe --configfile demo_csv/demo_config.yaml -j -d demo
```

respectively.


The result folder is created in the demo folder where you find the output files, as described above. 
<!-- The incidence plot of the demo sample should look like this: -->

<!-- ![alt text](https://github.com/KleistLab/GInPipe/blob/main/demo/demo_result.png) -->

[^1]: https://github.com/robert-koch-institut/SARS-CoV-2-Sequenzdaten_aus_Deutschland
[^2]: https://github.com/robert-koch-institut/SARS-CoV-2-Infektionen_in_Deutschland

## Running the pipeline for SARS-CoV-2 

To run the pipeline on real data, such as COVID sequences from GISAID or dna profiles from a [covSonar](https://github.com/rki-mf1/covsonar) database, download the according sequence data as well as the reference sequence and adapt the paths in the config.yaml as explained in [Initialization](#initialization).

## Reference
<a id="1">[1]</a>
Smith, M. R. and Trofimova, M., et al. (2021). Rapid incidence estimation from SARS-CoV-2 genomes reveals decreased case detection in Europe during summer 2020. Nature Communications  12, 6009. https://doi.org/10.1038/s41467-021-26267-y
