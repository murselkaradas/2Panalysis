#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=05:00:00
#SBATCH --mem=60GB
#SBATCH --job-name=SVDCompression
#SBATCH --mail-type=END
#SBATCH --mail-user=karadm01@nyumc.org
#SBATCH --output=slurm_%j.out
module purge
module load matlab/R2018a


matlab -nodisplay -r "SVDCompression('/gpfs/scratch/karadm01/2Pdata/MK18881/210707/aligned/MK18881_210707_rightglom8dodors_00001_00001.tif')"

exit

