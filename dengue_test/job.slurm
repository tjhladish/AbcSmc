#!/bin/bash

#SBATCH -J abc                     # Job Name
#SBATCH -o abc.o%j                 # Output file name (%j expands to jobID)
#SBATCH -e abc.e%j                 # Error file name (%j expands to jobID)
#SBATCH -n 1024                         # Total number of mpi tasks requested
#SBATCH -p normal                       # Queue (partition) name -- normal, development, etc.
#SBATCH -t 24:00:00                     # Run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=tjhladish@gmail.com # Email notification address (UNCOMMENT)
#SBATCH --mail-type=all                 # Email notification address (UNCOMMENT)
module load gsl
ibrun tacc_affinity ./abc_mpi abc_config.json             # Run the MPI executable named a.out
