#!/bin/bash
#SBATCH --job-name=bcl2fastq
#SBATCH --time=23:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH --output=%x-%j.out

#usage: sbatch ~/cre/bcl2fastq.sh <run_folder> <SampleSheet.csv>

#run: base run-folder where RunInfo.xml , SampleSheet.csv are present along with Data/Intensities/BaseCalls/L<lane folders>
#FASTQ output will be placed in $run/Data/Intensities/BaseCalls/

run=$1
sample_sheet=$2

if [ -z $run ]; then
	run=$1;
fi;

#pass this if the sample_sheet name is not under run folder/named differently
if [[ -z $sample_sheet || -n $2 ]]; then
	sample_sheet=$2;
fi;

arg=" -R ${run} ";
if [ ! -z $sample_sheet ]; then 
	arg=${arg}" --sample-sheet ${sample_sheet}";
fi;



#2.19 needs 20G of RAM for 10 threads
# sometimes 20G is not enough - I've got an error
# bcl2fastq is designed up to 32G
module load bcl2fastq

#I use no-lane splitting, and the result fastq files could not be uploaded to basespace
echo "CMD: bcl2fastq $arg -p  10 --no-lane-splitting"
bcl2fastq $arg -p  10 --no-lane-splitting 

