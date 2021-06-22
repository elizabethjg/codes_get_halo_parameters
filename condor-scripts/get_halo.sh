#!/bin/bash

# load local profile
. /etc/profile

# setup gsl

export OMP_NUM_THREAD=4

"PATH_TO_EXECUTABLE"/get_halo_props_PIC_v2 '/pnfs/pic.es/data/vo.paus.pic.es/mice/tape/raw/production/N4096_L3072_LC3/HaloCat/MICE_N4096_L3072_LC_b0p2/lightconedir_129/halo_detail_8_5' 'OUTPUT_PATH'
