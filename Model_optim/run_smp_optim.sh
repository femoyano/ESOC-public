#!/bin/sh
#BSUB -a openmp
#BSUB -q fat
#BSUB -W 30:00
#BSUB -o RUN-%J.out
#BSUB -n 30
#BSUB -R span[hosts=1]
Rscript --slave "run_optim_smp.R"

