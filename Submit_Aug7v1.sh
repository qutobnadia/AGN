#!/bin/bash
#SBATCH -n 48              # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 02-08:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared
#SBATCH --mem 100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nadia.qutob@cfa.harvard.edu
#SBATCH -o myoutput.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors.err  # File to which STDERR will be written, %j inserts jobid

module load gcc/10.2.0-fasrc01
module load hdf5/1.10.7-fasrc01
module load python/3.10.9-fasrc01

python GP_July10_t95_3000_tor4_mask5_Mg2.py & 
python GP_July10_t95_3000_tor4_mask200_Mg2.py & 
python GP_July10_t95_3000_tor4_lower_mask5_Mg2.py &
python GP_July10_t95_3000_tor4_lower_mask200_Mg2.py & 

wait 
