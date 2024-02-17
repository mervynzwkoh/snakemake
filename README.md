Setting up environment to run snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

1. Download miniforge
   - https://github.com/conda-forge/miniforge#mambaforge
2. Create new environment installed with snakemake (enter code below into miniforge cli)
   - "mamba create -c conda-forge -c bioconda -n snakemake snakemake"
3. Activate the environment (enter code below into miniforge cli)
   - "mamba activate snakemake"
   - check if snakemake is installed by entering "snakemake --version"
4. Git pull the workflow to a local directory
   - https://github.com/mervynzwkoh/snakemake.git

Executing the workflow

5. Edit the config file in the workflow
   - enter path of directory that contains all the fasta/vcf/msa file
   - the directory should include (ensure files and directories are renamed accordingly):
     - msa.fasta
     - MasterList: csv with accession numbers and release dates
     - fastavcfloc: directory that contains two more directories, "fasta" and "vcf", within those directories are the fasta and vcf files
       - the names of the fasta and vcf files should just be the accession number ie. "accessionnum.fasta" or "accessionnum.vcf"
6. Toggle to the snakemake directory in miniforge interface
7. Enter the following command to run the pipeline
   - "snakemake res/7seqTrack.done --cores 4 --sdm conda"
