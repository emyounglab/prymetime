#!/bin/bash
#SBATCH -J eng_sig_figure

source activate /home/jhcollins/anaconda2/envs/r_env

cp eng_sig_figure.R "$1"
cd "$1"

cp "$1"_final.bed genome.bed
cp "$1"_blastn.gff blastn.gff
Rscript eng_sig_figure.R "$1"
