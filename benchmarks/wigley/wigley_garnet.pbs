#!/bin/bash
#PBS -A ERDCV00898ENQ
#PBS -q standard
#PBS -l walltime=072:00:00
#PBS -l ncpus=4096
##PBS -q debug
##PBS -l walltime=001:00:00
##PBS -l ncpus=512
#PBS -N wigley
#PBS -j oe
#PBS -l application=proteus
#PBS -m eba
#PBS -M cekees@gmail.com
source /opt/modules/default/etc/modules.sh
source ${PROTEUS}/envConfig/garnet.gnu.bash
cd $PBS_O_WORKDIR
export MPICH_RANK_REORDER_METHOD=2
mkdir $WORKDIR/wigley.$PBS_JOBID
#note, garnet has about 62 GB of usable memory per node.if using -N 32 pass -M 1.9 to proteus, if -N 16 pass -M 3.8, etc.
#aprun -n 256 -N 16 -d 2 parun wigley_so.py -l 7 -v -M 3.75 -O ../inputTemplates/petsc.options.asm -D $WORKDIR/wigley.$PBS_JOBID
aprun -n 1024 -N 8 -S 2 parun wigley_so.py -l 7 -v -m -O ../inputTemplates/petsc.options.schur_upper_a11_asm_boomeramg -D $WORKDIR/wigley.$PBS_JOBID
