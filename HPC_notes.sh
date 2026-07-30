# this file is to store notes for login and job submission in HPC cluster/RDS (research data store)

# when off-campus use zscaler

# login
ssh lvd25@login.cx3.hpc.imperial.ac.uk

# To begin using Anaconda, you need to install miniforge
module load miniforge/3

# check modules that are loaded
ml list

ml tools/prod #always has to be loaded, however, I think its there on default, so no need to do sth

# subsequent loading of miniforge
eval "$(~/miniforge3/bin/conda shell.bash hook)"

# activate or create new environment
conda create -n projectx #create new environment
conda activate R43_240426 #activate env. (conda/source)
conda install r -c r-png # install R in the new environment, I guess this only has to be done in the very beginning
conda deactivate #deactivate current environment

# connect to RDS via Finder on Mac: Finder --> go --> connect to server --> use lvd25 as login
# copying files from and to the cluster in lecture: Demystifying HPC

# job script
# i am using this script as a wrapper for the run: qsub_cx3_R_R43_240426_1c32g8h.sh
# before is use this script, i have to hae the directory defined as e.g.
cd $HOME/AST_scRNAseq_TREM2/LD_X_Thesis_Presentatioon_Figures_scripts


qsub -v SCRIPT="LD_X04_B_characterisation_plots.R" ../qsub_cx3_R_R43_240426_1c32g8h.sh



# general commands
ls # see files within the directory
qstat -u lvd25 # check jobs
qdel job_ID # delete a job

exit # exit the cluster