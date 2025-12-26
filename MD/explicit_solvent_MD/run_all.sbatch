#!/bin/bash

#SBATCH --job-name=all_md
#SBATCH --gres=gpu:1
#SBATCH --mem=50GB
#SBATCH --time=48:00:00

module load amber/openmpi/intel/22.03

pmemd -O -i 1min.in\
            -o 1min.out -p dss.parm7 -c dss.rst7 -r 1min.rst7\
            -inf 1min.info -ref dss.rst7 -x mdcrd.1min

pmemd.cuda -O -i 2mdheat.in\
            -o 2mdheat.out -p dss.parm7 -c 1min.rst7 -r 2mdheat.rst7\
            -inf 2mdheat.info -ref 1min.rst7 -x mdcrd.2mdheat

pmemd.cuda -O -i 3md.in\
            -o 3md.out -p dss.parm7 -c 2mdheat.rst7 -r 3md.rst7\
            -inf 3md.info -ref 2mdheat.rst7 -x mdcrd.3md

pmemd.cuda -O -i 4md.in\
            -o 4md.out -p dss.parm7 -c 3md.rst7 -r 4md.rst7\
            -inf 4md.info -ref 3md.rst7 -x mdcrd.4md

pmemd -O -i 5min.in\
            -o 5min.out -p dss.parm7 -c 4md.rst7 -r 5min.rst7\
            -inf 5min.info -ref 4md.rst7 -x mdcrd.5min

pmemd.cuda -O -i 6md.in\
            -o 6md.out -p dss.parm7 -c 5min.rst7 -r 6md.rst7\
            -inf 6md.info -ref 5min.rst7 -x mdcrd.6md

pmemd.cuda -O -i 7md.in\
            -o 7md.out -p dss.parm7 -c 6md.rst7 -r 7md.rst7\
            -inf 7md.info -ref 6md.rst7 -x mdcrd.7md

pmemd.cuda -O -i 8md.in\
            -o 8md.out -p dss.parm7 -c 7md.rst7 -r 8md.rst7\
            -inf 8md.info -ref 7md.rst7 -x mdcrd.8md

pmemd.cuda -O -i 9md.in\
            -o 9md.out -p dss.parm7 -c 8md.rst7 -r 9md.rst7\
            -inf 9md.info -ref 8md.rst7 -x mdcrd.9md
~
