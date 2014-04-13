#!/bin/bash

#SBATCH -J abc_dice                     # Job Name
#SBATCH -o abc_dice.o%j                 # Output file name (%j expands to jobID)
#SBATCH -e abc_dice.e%j                 # Error file name (%j expands to jobID)
#SBATCH -n 256                          # Total number of mpi tasks requested
#SBATCH -p development                  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 00:20:00                     # Run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=tjhladish@gmail.com # Email notification address (UNCOMMENT)
#SBATCH --mail-type=all                 # Email notification address (UNCOMMENT)
ibrun ./abc_mpi abc_config.json             # Run the MPI executable named a.out
