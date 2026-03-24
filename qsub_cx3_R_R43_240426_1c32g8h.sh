#!/bin/bash
#PBS -lselect=1:ncpus=16:mem=256gb:ompthreads=32
#PBS -lwalltime=8:00:00
#PBS -o std_out.txt
#PBS -e std_err.txt


# Change to the submission directory

cd "$PBS_O_WORKDIR"

#submission notification
MESSAGE="##### Submitting R Script $SCRIPT... "
echo $MESSAGE

#load R environment
module load miniforge/3

eval "$(~/miniforge3/bin/conda shell.bash hook)"

#activate environment
conda activate R43_240426

echo "CONDA_PREFIX=$CONDA_PREFIX"
which R
which Rscript

#run R script
~/miniforge3/envs/R43_240426/bin/Rscript "$SCRIPT" >> "${SCRIPT}_out_err.txt" 2>&1


#completion notification
echo "R Script submission completed"