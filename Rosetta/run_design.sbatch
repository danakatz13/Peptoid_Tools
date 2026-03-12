#!/bin/bash

#SBATCH --job-name=hexamer_design
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=48:00:00
#SBATCH --mem=100GB

module load openmpi/gcc/4.0.5
module load gcc/10.2.0

srun ../../../../srv/release-prep/temp_src/rosetta.source.release-296/main/source/bin/rosetta_scripts.mpi.linuxgccrelease -parser:protocol design_only_complex.xml -nstruct 1000 -in:file:s ../../../inputs/hexamer_scaffold_docked/*pdb 
