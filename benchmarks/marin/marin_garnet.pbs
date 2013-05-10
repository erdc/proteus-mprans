#!/bin/bash
#PBS -A ERDCV00898ENQ
#PBS -l application=proteus
##PBS -l walltime=012:00:00
##PBS -l ncpus=1440
#PBS -l walltime=024:00:00
#PBS -l ncpus=512
#PBS -q standard
#PBS -N marin
#PBS -j oe
#PBS -m eba
#PBS -M cekees@gmail.com
# nnodes = ncpus/32
source /opt/modules/default/etc/modules.sh
source ${PROTEUS}/envConfig/garnet.gnu.bash
cd $PBS_O_WORKDIR
mkdir $WORKDIR/marin.$PBS_JOBID
# n = N*nnodes = N*ncpus/32
#aprun -n 512 parun marin_so.py -l 7 -v --cacheArchive -O petsc.options --profile -D $WORKDIR/marin.$PBS_JOBID
aprun -n 256 -N 16 -S 4 parun marin_so.py -M 3.875 -l 7 -v -O ../inputTemplates/petsc.options.schur_upper_a11_asm_boomeramg -D $WORKDIR/marin.$PBS_JOBID
