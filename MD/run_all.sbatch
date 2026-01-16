#!/bin/bash
#SBATCH --job-name=amber_md
#SBATCH --output=pipeline_%j.log
#SBATCH --error=pipeline_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --time=96:00:00
#SBATCH --account=torch_pr_489_general

module purge
module load amber/24.03
SYSNAME=$(basename "$PWD")

echo "Running AMBER pipeline for system: ${SYSNAME}"

PARM_FILE="${SYSNAME}.parm7"
COOR_FILE="${SYSNAME}.rst7"

if [[ ! -f "$PARM_FILE" || ! -f "$COOR_FILE" ]]; then
    echo "Error: Could not find $PARM_FILE or $COOR_FILE in current directory."
    exit 1
fi
echo "Starting Minimization..."

sander -O \
    -i min.in \
    -p "${PARM_FILE}" \
    -c "${COOR_FILE}" \
    -o min1.out \
    -r min1.rst \
    -inf min1.info

if [ $? -ne 0 ]; then
    echo "Minimization FAILED. Stopping pipeline."
    exit 1
fi
pmemd.cuda -O \
    -i heat.in \
    -p "${PARM_FILE}" \
    -c min1.rst \
    -o heat.mdout \
    -x heat.nc \
    -r heat.rst7

if [ $? -ne 0 ]; then
    echo "Heating FAILED. Stopping pipeline."
    exit 1
fi

echo "Starting Production MD..."

pmemd.cuda -O \
    -i md.in \
    -p "${PARM_FILE}" \
    -c heat.rst7 \
    -o md1.mdout \
    -x md1.nc \
    -r md1.rst7

if [ $? -ne 0 ]; then
    echo "Production MD FAILED."
    exit 1
fi

echo "Pipeline Finished Successfully for ${SYSNAME}!"

