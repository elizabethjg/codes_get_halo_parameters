#!/bin/bash
#SBATCH --job-name=hpar_z5
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

for i in `seq $SLURM_ARRAY_TASK_ID 1 $SLURM_ARRAY_TASK_MAX`;
do
    if [ $i -eq 1 ];
    then
        srun ./get_halo_props_PIC_v3 /mnt/simulations/SIDM_simus/Lentes/V2/CDM/Snap43/ 143 0.51 ../SIDM/halo_props/halo_props_cdm_z51
        ##srun ./get_halo_props_PIC_v3 /mnt/simulations/SIDM_simus/Lentes/V2/CDM/Snap39/ 75 0.96 ../SIDM/halo_props/halo_props_cdm_z96
    else
        srun ./get_halo_props_PIC_v3 /mnt/simulations/SIDM_simus/Lentes/V2/SIDM1/Snap43/ 143 0.51 ../SIDM/halo_props/halo_props_cdm_z51
        ##srun ./get_halo_props_PIC_v3 /mnt/simulations/SIDM_simus/Lentes/V2/SIDM1/Snap39/ 75 0.96 ../SIDM/halo_props/halo_props_cdm_z96
    fi
done
