#!/bin/bash
#SBATCH --job-name=sawtooth
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=ec288-wetting
#SBATCH --qos=lowpriority

#SBATCH --ntasks=4
#SBATCH --cpus-per-task=9

# Clean inherited OpenMP affinity (CRITICAL)
unset OMP_PLACES
unset OMP_PROC_BIND

# Let Slurm define the threading
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_DYNAMIC=FALSE

echo "SLURM_NTASKS=$SLURM_NTASKS"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"

srun --cpu-bind=cores ./run.exe
