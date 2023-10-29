#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J eng_sig_alitv

EXECDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"
PREFIX=$(basename "$1")

cp "$EXECDIR"/alitv.yml "$1"
cp "$2" "$1"
cd "$1"
cp "$2" reference.fasta
cp "$PREFIX"_final.fasta prymetime_assembly.fasta

perl ~/AliTV-perl-interface-1.0.6/bin/alitv.pl --project "$PREFIX" alitv.yml
