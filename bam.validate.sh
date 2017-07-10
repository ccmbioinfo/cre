#!/bin/bash

#PBS -l walltime=5:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

if [ -z $bam ]
then
    bam=$1
fi

picard ValidateSamFile I=$bam > $bam.check