#!/bin/bash
##SBATCH --partition=purley-cpu   
#SBATCH --partition=all-x86-cpu
#SBATCH --time=23:00:00
#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --mem=6gb
#SBATCH --output=/home/jywang02/multest/fix_output.out
#SBATCH --error=/home/jywang02/multest/fix_error.err
#SBATCH --export=all
#SBATCH --exclude=purley-x86-cpu01

#module use /home/app/modulefiles
#module load askapsoft/cpu-cloud-dingo-jacal-master

source /home/jywang02/.bashrc

#module use /home/app/modulefiles
#module load python/cpu-3.6.5
#module load casacore/cpu-py3.6.5-3.1.0

cd /home/jywang02/multest/   
srun --mpi=pmi2 --partition=all-x86-cpu -N 1 python -m dlg.deploy.start_dlg_cluster -l /home/jywang02/multest -L /home/jywang02/multest/graph/fix61.graph --part-algo mysarkar --ssid HelloWorld-`date --iso-8601=seconds` --algo-param max_cpu=1 --pg-modifiers "" -d 

#iso-8601
#iso-8859-1 
