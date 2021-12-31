#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=18:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --mem=100G   # maximum memory per node
#SBATCH --job-name="blastMC30fHalfHIV"
#SBATCH --mail-user=blitt@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load parallel 
cd /work/LAS/jernigan-lab/Ben/

base_path=/work/LAS/jernigan-lab/Ben
rm -r $base_path/ProjectHIV/rawOutput
mkdir $base_path/ProjectHIV/rawOutput 

ls $base_path/ProjectHIV/HIVProteinsFirstHalf/ | parallel "$base_path/blastMC30V2/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp -query $base_path/ProjectHIV/HIVProteinsFirstHalf/{} -db $base_path/ProjectHIV/swissprot/swissprot -outfmt '6 qacc sacc evalue qstart qend sstart send' -matrix BLOSUM90 -out $base_path/ProjectHIV/rawOutput/{.}MC30BlastOut.tsv -max_target_seqs 25000"

