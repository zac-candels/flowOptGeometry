#!/bin/bash

# Loop through all directories
for dir in test_CA_60 test_CA_120; do
  # Check if the directory exists and contains an executable file
  if [[ -d "$dir" && -x "$dir/run.exe" ]]; then
    # Create a SLURM batch script
    cat <<EOF > "${dir}job.slurm"
#!/bin/bash
#SBATCH --job-name=${dir%/}_job   # Job name
#SBATCH --output=${dir}output.log # Standard output log file
#SBATCH --error=${dir}error.log   # Standard error log file
#SBATCH --time=48:00:00           # Time limit hrs:min:sec
#SBATCH --partition=standard      # Partition name
#SBATCH --account=sc139-wetting
#SBATCH --qos=lowpriority
#SBATCH --ntasks=1                # Number of tasks (processes)
#SBATCH --cpus-per-task=8        # Number of CPU cores per task

# Change to the directory
cd $dir

# Run the executable
srun ./run.exe
EOF

    # Submit the job script
    sbatch "${dir}job.slurm"
  else
    echo "Skipping $dir: Not a directory or no executable found."
  fi
done

echo "Job submission complete."

