#!/bin/bash
#SBATCH --job-name=sawtooth
#SBATCH --output=output.log 
#SBATCH --error=error.log   
#SBATCH --time=24:00:00           
#SBATCH --partition=standard      
#SBATCH --account=ec288-wetting
#SBATCH --qos=lowpriority
#SBATCH --ntasks=1              
#SBATCH --cpus-per-task=250      


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_ENV=TRUE   


srun ./run.exe
