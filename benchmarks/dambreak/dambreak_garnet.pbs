#!/bin/bash
#PBS -A ERDCV00898R40
##PBS -l walltime=024:00:00
##PBS -l ncpus=512
#PBS -l walltime=001:00:00
#PBS -l select=4:ncpus=32:mpiprocs=16
##PBS -l walltime=048:00:00
##PBS -l ncpus=1024
##PBS -q standard
#PBS -q debug
#PBS -N dambreak
#PBS -j oe
#PBS -l application=proteus
#PBS -V
#PBS -m eba
#PBS -M cekees@gmail.com
source /opt/modules/default/etc/modules.sh
module swap PrgEnv-pgi  PrgEnv-gnu
module swap xt-libsci acml
#source ${PROTEUS}/envConfig/garnet.gnu.bash
cd $PBS_O_WORKDIR
mkdir $WORKDIR/dambreak.$PBS_JOBID
#-n total number of cores, -N cores per node, -S cores per numa node (4 numa nodes per node)
#aprun -n 32 parun dambreak_so.py -l 7 -v -G -m -O ../inputTemplates/petsc.options.schur_upper_a11_asm_boomeramg -D $WORKDIR/dambreak.$PBS_JOBID
aprun -n 64  parun dambreak_so.py -F -l 7 -v -O ../inputTemplates/petsc.options.schur_upper_a11_asm_boomeramg -D $WORKDIR/dambreak.$PBS_JOBID
#aprun -n 256 -N 8 -S 2 parun dambreak_so.py -l 7 -v -G -M 8.0 -O ../inputTemplates/petsc.options.schur_upper_a11_asm_boomeramg -D $WORKDIR/dambreak.$PBS_JOBID
