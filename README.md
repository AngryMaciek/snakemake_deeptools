# Snakemake pipeline for RNA-Seq data analysis with deepTools
*Maciej_Bak  
Swiss_Institute_of_Bioinformatics*

[deepTools](https://deeptools.readthedocs.io/en/develop/) is a very nice toolset for exploring RNA-Seq data.  
This repository is a snakemake workflow that is based on the example usage from the deepTools manual:  
https://deeptools.readthedocs.io/en/develop/content/example_usage.html  
My aim was to develop an automatized and reproducible pipeline for my research which I would now happily share with the community :) 

## Snakemake pipeline execution
Snakemake is a workflow management system that helps to create and execute data processing pipelines. It requires Python 3 and can be most easily installed via the bioconda package from the anaconda cloud service.

### Step 1: Download and installation of Miniconda3
Linux:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  source .bashrc
  ```

macOS:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh
  source .bashrc
  ```

### Step 2: Pandas and Snakemake installation

To execute the workflow one would require pandas python library and snakemake workflow menager.  
Unless a  specific snakemake version is specified explicitly it is most likely the best choice to install the latest versions:
  ```bash
  conda install -c conda-forge pandas
  conda install -c bioconda snakemake
  ```

In case you are missing some dependancy packages please install them first (with `conda install ...` as well).

### Step 3: Pipeline execution
Specify all the required information (input/output/parameters) in the config.yaml  
The main input to the pipeline is a simple design table which has to have the following format:

sample  bam  
[sample_name] [path_to_bam]  
...

Where:
* Each row is a sequencing sample.
* All the bam files need to have a different name regardless of their location.
* Design table might have more columns than these above.

Apart from the design table the pipeline requires a FASTA-formatted genome file.

Once the metadata are ready write a DAG (directed acyclic graph) into dag.pdf:
  ```bash
  bash snakemake_dag_run.sh
  ```

There are two scripts to start the pipeline, depending on whether you want to run locally or on a SLURM computational cluster. In order to execute the workflow snakemake automatically creates internal conda virtual environments and installs software from anaconda cloud service. For the cluster execution it might be required to adapt the 'cluster_config.json' and submission scripts before starting the run.
  ```bash
  bash snakemake_local_run_conda_env.sh
  bash snakemake_cluster_run_conda_env.sh
  ```

## License

Apache 2.0
