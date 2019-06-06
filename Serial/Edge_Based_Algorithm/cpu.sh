#!/bin/bash
#SBATCH -n 1 # Number of cores 
#SBATCH --ntasks-per-node=1 # Ensure that all cores are on one machine 
#SBATCH -t 0-02:55 # Runtime in D-HH:MM 
#SBATCH -p cpu # Partition to submit to
./main ../../../DataSet/movie_taste_net.col
