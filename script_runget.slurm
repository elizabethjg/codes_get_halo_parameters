#!/bin/bash
#SBATCH --job-name=halo_par
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=3-00:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

. /etc/profile
source $HOME/.bashrc

module load gcc gsl hdf5 openmpi

OMP_NUM_THREADS=28 srun --ntasks=1 --cpus-per-task=28 ./get_halo_props_PIC_v3 /mnt/simulations/SIDM_simus/Lentes/V2/CDM/Snap39/ 74 0.96 ../SIDM/halo_props/halo_props_cdm_z96 &
OMP_NUM_THREADS=28 srun --ntasks=1 --cpus-per-task=28 ./get_halo_props_PIC_v3 /mnt/simulations/SIDM_simus/Lentes/V2/SIDM/Snap39/ 74 0.96 ../SIDM/halo_props/halo_props_cdm_z96 &
wait
