#!/bin/bash
#SBATCH --job-name=hpar_z0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=3-00:00:00
#SBATCH --array=1-2
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

. /etc/profile
source $HOME/.bashrc

module load gcc gsl hdf5 openmpi

if [ $SLURM_ARRAY_TASK_ID -eq 1 ];
then
    srun ./get_halo_props_PIC_v3 /mnt/projects/lensing/SIDM_project/Lentes/Eli_Agus/snapshot_050/CDM/ 214 0.0 ../SIDM/halo_props/halo_props_fof_cdm_z0
    ##srun ./get_halo_props_PIC_v3 /mnt/simulations/SIDM_simus/Lentes/V2/CDM/Snap39/ 75 0.96 ../SIDM/halo_props/halo_props_cdm_z96
fi
if [ $SLURM_ARRAY_TASK_ID -eq 2 ];
then
    srun ./get_halo_props_PIC_v3 /mnt/projects/lensing/SIDM_project/Lentes/Eli_Agus/snapshot_050/SCDM/ 214 0.0 ../SIDM/halo_props/halo_props_fof_sidm1_z0
    ##srun ./get_halo_props_PIC_v3 /mnt/simulations/SIDM_simus/Lentes/V2/SIDM1/Snap39/ 75 0.96 ../SIDM/halo_props/halo_props_sidm1_z96
fi
