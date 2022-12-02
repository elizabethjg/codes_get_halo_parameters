#!/bin/bash
#SBATCH --job-name=test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=3-00:00:00
#SBATCH --array=1-2
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=lpereyra@unc.edu.ar

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

. /etc/profile
###source $HOME/.bashrc

module load gcc gsl hdf5 openmpi

if [ $SLURM_ARRAY_TASK_ID -eq 1 ];
then
  	srun ./get_halo_props_PIC_v3 /mnt/projects/lensing/SIDM_project/Lentes/Eli_Agus/snapshot_050/rockstar/matcheados/CDM/ 242 0.0 extend_halo_propsv2_rock2_match_cdm_z0
fi
if [ $SLURM_ARRAY_TASK_ID -eq 2 ];
then
  	srun ./get_halo_props_PIC_v3 /mnt/projects/lensing/SIDM_project/Lentes/Eli_Agus/snapshot_050/rockstar/matcheados/SIDM1/ 242 0.0 extend_halo_propsv2_rock2_match_sidm1_z0
fi
