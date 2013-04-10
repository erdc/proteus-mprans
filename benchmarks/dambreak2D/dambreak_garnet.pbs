#!/bin/bash
#PBS -A ERDCV00898ENQ
#PBS -l walltime=048:00:00
#PBS -l ncpus=32
#PBS -q standard
#PBS -N dambreak2D
#PBS -j oe
#PBS -l application=proteus
#PBS -V
#PBS -m eba
#PBS -M cekees@gmail.com
cd $PBS_O_WORKDIR
mkdir $WORKDIR/dambreak2D.$PBS_JOBID
#aprun -n <proc>               when fully using all the (16) cores
#aprun -N 12 -S 3 -n <procs>   when putting 12 processes on the node
#aprun -N 8 -S 2 -n <procs>    when putting 8 processes on the node
#aprun -N 4 -S 1 -n <procs>    when putting 4 processes on the node
#aprun -n 32 parun dambreak_so.py -l 7 -v -G -O ../inputTemplates/petsc.options.asm -D $WORKDIR/dambreak2D.$PBS_JOBID
aprun -n 32 parun dambreak_so.py -l 7 -v -G -M 1.9 -B -O ../inputTemplates/petsc.options.schur_upper_a11_asm_boomeramg -D $WORKDIR/dambreak2D.$PBS_JOBID
