#!/bin/bash
#SBATCH --nodes=1                                                     # Number of nodes
#SBATCH --ntasks-per-node=3                                           # Run a # of tasks per node
#SBATCH --cpus-per-task=1                                             # Number of CPU cores per task
#SBATCH --mem=122gb                                                   # Job memory request
#SBATCH --time=15-00:00:00                                            # Time limit days-hrs:min:sec
#SBATCH --partition=kbs                                               # Must specify your partition
#SBATCH --mail-type=ALL                                               # Mail events
#SBATCH --mail-user=burrill.haley@ku.edu                           # Where to send mail
#SBATCH --output=/home/h128b273/work/globus_climush-sequences/soil-litt-redo/Analyses/indicatorsp.log                 # Standard outpu$
#SBATCH --job-name=indicatorsp                               # Job name

cd /home/h128b273/work/globus_climush-sequences/soil-litt-redo/Analyses/
module load R/4.4

Rscript ~/work/globus_climush-sequences/soil-litt-redo/Analyses/indicatorsp.R