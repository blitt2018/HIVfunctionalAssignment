#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=10:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --mem=40G   # maximum memory per node
#SBATCH --job-name="NRtoUniIDs_1"
#SBATCH --mail-user=blitt@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge 
module load python/3.7.7-b5s6jni
module load py-numpy/ 
module load py-scipy/ 
module load py-pandas/ 

#paths to accessory/reference files needed 
base_path=/work/LAS/jernigan-lab/Ben
uni_go_mapping_file=$base_path/ProjectHIV/swissprot_ID_GOs.tab

#paths to input/output directories (output for one thing is input into the next) 
blast_raw_output_dir=$base_path/ProjectHIV/rawOutput/
uni_01_dir=$base_path/ProjectHIV/uniIDOutput01/ 
uni_001_dir=$base_path/ProjectHIV/uniIDOutput001/ 
GO_01_dir=$base_path/ProjectHIV/ProtGOs01/ 
GO_001_dir=$base_path/ProjectHIV/ProtGOs001/
enriched_01_dir=$base_path/ProjectHIV/enrichedList01/ 
enriched_001_dir=$base_path/ProjectHIV/enrichedList001/
enriched_01_CSV=$base_path/ProjectHIV/enrichedList01/01lengths.csv
enriched_001_CSV=$base_path/ProjectHIV/enrichedList001/001lengths.csv

#paths to python files 
evalCutScript=$base_path/ProjectHIV/makeEvalCut.py 
uniToGOScript=$base_path/ProjectHIV/uniIDsToGOsSwissprot.py
fisherTestScript=$base_path/ProjectHIV/applyFisherTestSwissprot.py 
enrichmentCSVScript=$base_path/ProjectHIV/getEnrichedCSV.py

 
python3 $evalCutScript .01 $blast_raw_output_dir $uni_01_dir
python3 $evalCutScript .001 $blast_raw_output_dir $uni_001_dir

python3 $uniToGOScript $uni_go_mapping_file  $uni_01_dir $GO_01_dir
python3 $uniToGOScript $uni_go_mapping_file  $uni_001_dir $GO_001_dir
 
python3 $fisherTestScript $uni_go_mapping_file $GO_01_dir
python3 $fisherTestScript $uni_go_mapping_file $GO_001_dir

mkdir $enriched_01_dir
mkdir $enriched_001_dir
mv $GO_01_dir*GOLIST* $enriched_01_dir 
mv $GO_001_dir*GOLIST* $enriched_001_dir 

python3 $enrichmentCSVScript $enriched_01_dir $enriched_01_CSV
python3 $enrichmentCSVScript $enriched_001_dir $enriched_001_CSV  
