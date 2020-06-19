#!/bin/bash
#SBATCH -J eng_sig_figure

EXECDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"
PREFIX=$(basename "$1")

cp "$EXECDIR"/eng_sig_figure.R "$1"
cd "$1"

cp "$PREFIX"_final.bed genome.bed
cp "$PREFIX"_blastn.gff blastn.gff
Rscript eng_sig_figure.R "$1"
