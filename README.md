Setting up environment to run snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

1. Download miniforge
   - https://github.com/conda-forge/miniforge#mambaforge
2. Create new environment installed with snakemake (enter code below into miniforge cli)
   - "mamba create -c conda-forge -c bioconda -n snakemake snakemake"
3. Activate the environment (enter code below into miniforge cli)
   - "mamba activate snakemake"
   - check if snakemake is installed by entering "snakemake --version"
4. Clone the workflow to a local directory
   - https://github.com/mervynzwkoh/snakemake.git

Executing the workflow

5. Toggle to the snakemake directory in miniforge interface
6. Enter the following command to run the pipeline
   - "snakemake res/7seqTrack.done --cores 4 --sdm conda --configfile config/config.yaml --config src={path of directory containing all fasta/vcf/msa files}"
   - path name should contain backslashes "/" instead of forward slash "\"
   - the directory should include (ensure files and directories are renamed accordingly):
     - msa.fasta
     - MasterList: csv with accession numbers and release dates
       - must contain column headers "Accession Number" and "Release Date"
     - fastavcfloc: directory that contains two more directories, "fasta" and "vcf", within those directories are the fasta and vcf files
       - the names of the fasta and vcf files should just be the accession number ie. "accessionnum.fasta" or "accessionnum.vcf"
