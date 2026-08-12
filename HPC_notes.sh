# HPC login and job-submission notes.

# Off campus, connect through Zscaler.

# login
ssh lvd25@login.cx3.hpc.imperial.ac.uk

# Initial Miniforge setup
module load miniforge/3

# check modules that are loaded
ml list

ml tools/prod

# subsequent loading of miniforge
eval "$(~/miniforge3/bin/conda shell.bash hook)"

# Setup environment
conda create -n projectx
conda activate R43_240426
conda install r -c r-png
conda deactivate

# connect to RDS via Finder on Mac: Finder -> go -> connect to server -> use lvd25 as login
# copying files from and to the cluster in lecture: Demystifying HPC

# Submit from the script directory using the PBS wrapper
cd $HOME/AST_scRNAseq_TREM2/LD_X_Thesis_Presentatioon_Figures_scripts


qsub -v SCRIPT="LD_X04_B_characterisation_plots.R" ../qsub_cx3_R_R43_240426_1c32g8h.sh



# general commands
ls # see files within the directory
qstat -u lvd25 # check jobs
qdel job_ID # delete a job

exit # exit the cluster