#!/bin/bash 
#SBATCH -n 16 
#SBATCH -N 1
#SBATCH -o M1.10.ST.out      # File to which STDOUT will be written
#SBATCH -e M1.10.ST.err      # File to which STDERR will be written
#SBATCH -p sched_mit_hill
#SBATCH --mem 60G
#SBATCH --mail-type END 
#SBATCH --mail-user rhitchin@mail.einstein.yu.edu

module add engaging/clc-assembly-cell
module load engaging/prodigal/2.6.1 

#QUALITY CHECKING
#ls whatever.ncbi.output | awk '{print "clc_quality_trim -f 33 -o "$1"_trimmed_singletons.fq -p  "$1"_trimmed_paired.fq -r  "$1}' 
clc_quality_trim -f 33 -o  m1.10.ST_trimmed_singletons.fq -p  m1.10.ST_trimmed_paired.fq -r  m1.10.ST

#RUN CLC OVERLAP
#ls whatever.ncbi.output | awk '{print "clc_overlap_reads -r "$1"_trimmed_paired.fq -j  "$1".clc.overlapped.paired.fq -n "$1".clc.nonoverlapped.paired"}'
clc_overlap_reads -r m1.10.ST_trimmed_paired.fq -j  m1.10.ST.clc.overlapped.paired.fq -n m1.10.ST.clc.nonoverlapped.paired 

#ASSEMBLE OVERLAPPED CLC OUTPUT
#ls whatever.ncbi.output | awk '{print "clc_assembler -q  "$1".clc.overlapped.paired.fq "$1".clc.nonoverlapped.paired.fq "$1"_trimmed_singletons.fq  -o "$1"_trimmed_overlapped.paired.byclc"}'  
clc_assembler -q  m1.10.ST.clc.overlapped.paired.fq m1.10.ST.clc.nonoverlapped.paired.fq m1.10.ST_trimmed_singletons.fq     -o m1.10.ST_trimmed_overlapped.paired.byclc

#RUN Prodigal protein prediction on clc output
#ls whatever.ncbi.output | awk '{print "prodigal -p meta -i "$1"_trimmed_overlapped.paired.byclc.clc -o  "$1"_trimmed_overlapped.paired.byclc.clc.prodigal -a "$1"_trimmed_overlapped.paired.byclc.clc.faa"}'
prodigal -p meta -i m1.7.ST_trimmed_overlapped.paired.byclc.clc -o  m1.7.ST_trimmed_overlapped.paired.byclc.clc.prodigal -a m1.7.ST_trimmed_overlapped.paired.byclc.clc.faa 
