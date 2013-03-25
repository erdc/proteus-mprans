#!/bin/bash
#PBS -A ERDCV00898ENQ
#PBS -l walltime=048:00:00
#PBS -l ncpus=32
#PBS -q R368598
##PBS -q standard
#PBS -N flume2D
#PBS -j oe
#PBS -l application=proteus
#PBS -V
#PBS -m eba
#PBS -M cekees@gmail.com
##PBS -q R1537762 
cd $PBS_O_WORKDIR
mkdir $WORKDIR/flume2D.$PBS_JOBID
#aprun -n <proc>               when fully using all the (16) cores
#aprun -N 12 -S 3 -n <procs>   when putting 12 processes on the node
#aprun -N 8 -S 2 -n <procs>    when putting 8 processes on the node
#aprun -N 4 -S 1 -n <procs>    when putting 4 processes on the node
aprun -n 32 parun flume_so.py -l 7 -v -G -O ../inputTemplates/petsc.options\
.schur_upper_a11_asm_boomeramg -D $WORKDIR/flume2D.$PBS_JOBID
