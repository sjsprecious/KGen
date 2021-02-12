#!/bin/csh
#SBATCH  --job-name=play
#SBATCH  --nodes=1
#SBATCH  --ntasks-per-node=1
#SBATCH  --ntasks=1
#SBATCH  --gres=gpu:v100:1
#SBATCH --output=out.%j
#SBATCH --time=00:10:00
#SBATCH --account=NTDD0004
#SBATCH --reservation=TDD_4xV100
#SBATCH  --exclusive                        

cd $cwd

module purge
module load pgi/20.4
module load netcdf/4.7.4
module load cuda/10.1
module load openmpi/4.0.5

### Run program
make clean
make
