#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06-00:00:00
#SBATCH --mem=14G
#SBATCH --array=1-10
#SBATCH --mail-user=julianna.d.ianni@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="nosvd_0pt01_1"
#SBATCH --output=nosvd_0pt01_1_%A_%a.out

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
 
setpkgs -a matlab_r2015a

matlab -nodisplay -nosplash -r "POCSRR100sub($SLURM_ARRAY_TASK_ID,1,0.01,3,1,40000,0,3,1,1,1,0,0,'both',0,Inf,0)"

