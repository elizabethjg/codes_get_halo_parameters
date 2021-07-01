#!/bin/bash
#SBATCH --job-name=getparam
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=6:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


. /etc/profile
source $HOME/.bashrc

cd /mnt/projects/lensing/HALO_SHAPE/MICEv1.0/codes_get_halo_parameters/

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load gcc gsl

srun ./get_halo_props_PIC_v2 '../catalogs/halo_detail_8_5' '../catalogs/halo_props2_8_5_2'

