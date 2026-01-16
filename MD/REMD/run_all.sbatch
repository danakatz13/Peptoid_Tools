#!/bin/bash

#SBATCH --job-name=amber
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=72:00:00

module purge
module load amber/openmpi/intel/22.00

echo "Job started on $(date)"
echo "Running in directory: $(pwd)"

echo "=== MINIMIZATION ==="

sander -O \
  -i min.in \
  -p design_1.parm7 \
  -c design_1.rst7 \
  -o min1.out \
  -r min1.rst

if [ $? -ne 0 ]; then
  echo "Minimization failed"
  exit 1
fi

echo "=== EQUILIBRATION ==="

srun sander.MPI \
  -ng 12 \
  -groupfile equilibrate.groupfile

if [ $? -ne 0 ]; then
  echo "Equilibration failed"
  exit 1
fi

echo "=== REMD PRODUCTION ==="

srun sander.MPI \
  -ng 12 \
  -groupfile remd.groupfile

if [ $? -ne 0 ]; then
  echo "REMD failed"
  exit 1
fi

echo "Job finished successfully on $(date)"

