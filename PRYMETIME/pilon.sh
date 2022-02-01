#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=160G
#SBATCH -o pilon_%j.out
#SBATCH -e pilon_%j.err
#SBATCH -J pilon


# map illumina reads for pilon
#bwa mem -t 14 ~/nanopore_6/FEY_48/FEY_48_Flye_All_med_rac.fasta ~/illumina_2/FEY48_ill_com.fastq | samtools view - -Sb | samtools sort - -@14 -o ~/nanopore_6/FEY_48/FEY_48_Flye_All_med_rac_pil.bam
bwa mem -t ${N_THREADS} "$1" "$2" | samtools view - -Sb | samtools sort - -@ ${N_THREADS} -o "$3"
#samtools index ~/nanopore_6/FEY_48/FEY_48_Flye_All_med_rac_pil.bam
samtools index -@ ${N_THREADS} "$3"

# pilon polish
#pilon -Xmx160G --genome ~/nanopore_6/FEY_48/FEY_48_Flye_All_med_rac.fasta --bam ~/nanopore_6/FEY_48/FEY_48_Flye_All_med_rac_pil.bam --output ~/nanopore_6/FEY_48/FEY_48_Flye_All_med_rac_pil
#--threads listed as experimental for pilon
pilon -Xmx160G --genome "$1" --bam "$3" --output "$4"
