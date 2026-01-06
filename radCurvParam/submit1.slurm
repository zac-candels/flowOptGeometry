#!/bin/bash
#SBATCH --job-name=two-six
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=ec288-wetting
#SBATCH --qos=lowpriority

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1    
#SBATCH --cpus-per-task=288         
#SBATCH --exclusive

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=spread            # Still crucial

srun --cpu-bind=cores ./run.exe