# Setting up environment to run snakemake 
(https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

1. Download miniforge
   * https://github.com/conda-forge/miniforge#mambaforge
   * Linux:

     `wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"`
   * MacOS:

     `curl -fsSLo Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-$(uname -m).sh"`
   * Windows: download installation package from https://github.com/conda-forge/miniforge#mambaforge
  
2. Create new environment installed with snakemake (enter code below into miniforge cli)

`mamba create -c conda-forge -c bioconda -n snakemake snakemake`

3. Activate the environment (enter code below into miniforge cli)

`mamba activate snakemake`
   - check if snakemake is installed by entering
     `snakemake --version`
4. Clone the workflow to a local directory

`git clone https://github.com/mervynzwkoh/snakemake.git`

# Executing the workflow

5. Toggle to the snakemake directory in miniforge interface

  `cd {path to directory where snakemake workflow was cloned`
  
6.  Enter the following command to run the pipeline

`snakemake res/7seqTrack.done --cores 4 --sdm conda --configfile config/config.yaml --config src={path of directory containing input data}`
   - path name should contain backslashes "/" instead of forward slash "\"

# Input Data Preparation

All input data should be placed under one directory. The following input data are required (the naming of the files/directories should be strictly followed):
1. ***msa.fasta*** file
2. ***MasterList.csv*** file
3. ***fastavcfloc*** directory
   - contains two contains two more directories, ***fasta*** and ***vcf***, within those directories are the fasta and vcf files
   - the names of the fasta and vcf files must be the same as their respective accession numbers in the ***MasterList.csv*** file ie. ***{accessionnum}.fasta*** or ***{accessionnum}.vcf***
