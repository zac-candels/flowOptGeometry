#!/bin/bash
#SBATCH --job-name=noSizeGrad  # Job name
#SBATCH --output=output.log # Standard output log file
#SBATCH --error=error.log   # Standard error log file
#SBATCH --time=24:00:00           # Time limit hrs:min:sec
#SBATCH --partition=standard      # Partition name
#SBATCH --account=ec288-wetting
#SBATCH --qos=lowpriority
#SBATCH --ntasks=5              # Number of tasks (processes)
#SBATCH --cpus-per-task=280      # Number of CPU cores per task

# ===============================
# OpenMP configuration
# ===============================
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_ENV=TRUE   # Optional: print OpenMP environment info at startup

echo "===================================="
echo "Running LBM simulation with OpenMP"
echo "Allocated CPUs per task: $SLURM_CPUS_PER_TASK"
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "===================================="

# Run the executable
srun ./run.exe
