#!/bin/bash
#SBATCH --partition=purley-cpu
#SBATCH --time=23:00:00
#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --mem=300gb
#SBATCH --output=/home/jywang02/multest/make_output.out
#SBATCH --error=/home/jywang02/multest/make_error.err
##SBATCH --export=all
#SBATCH --exclude=purley-x86-cpu01



source /home/jywang02/.bashrc



cd /home/jywang02/multest/   
srun --mpi=pmi2 -N 1 python -m dlg.deploy.start_dlg_cluster -l /home/jywang02/multest -L /home/jywang02/multest/graph/makemodel3.graph --part-algo mysarkar --ssid make1-`date --iso-8601=seconds` --algo-param max_cpu=1 --pg-modifiers "" -d 

#iso-8601
#iso-8859-1 
