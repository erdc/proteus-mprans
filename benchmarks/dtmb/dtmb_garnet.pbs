#!/bin/bash
#PBS -A ERDCV00898ENQ
#PBS -q standard
#PBS -l walltime=072:00:00
#PBS -l ncpus=4096
#PBS -N dtmb
#PBS -j oe
#PBS -l application=proteus
#PBS -m eba
#PBS -M cekees@gmail.com
source /opt/modules/default/etc/modules.sh
source ${PROTEUS}/envConfig/garnet.gnu.bash
cd $PBS_O_WORKDIR
export MPICH_RANK_REORDER_METHOD=2
mkdir $WORKDIR/dtmb.$PBS_JOBID
#aprun -n 64  parun dtmb_so.py -l 7 -m -O ../inputTemplates/petsc.options.schur_upper_a11_asm_boomeramg -D $WORKDIR/dtmb.$PBS_JOBID
#aprun -n 128  parun dtmb_so.py -l 7 -m -O ../inputTemplates/petsc.options.schur_upper_a11_asm_boomeramg -D $WORKDIR/dtmb.$PBS_JOBID
#aprun -n 256 -N 16 -S 4  parun dtmb_so.py -l 7 -v -m -O ../inputTemplates/petsc.options.schur_upper_a11_asm_boomeramg -D $WORKDIR/dtmb.$PBS_JOBID
aprun -n 1024 -N 8 -S 2  parun dtmb_so.py -l 7 -v -m -O ../inputTemplates/petsc.options.schur_upper_a11_asm_boomeramg -D $WORKDIR/dtmb.$PBS_JOBID
