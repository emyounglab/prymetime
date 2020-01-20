#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=160G
#SBATCH -o illumina_merge_%j.out
#SBATCH -e illumina_merge_%j.err
#SBATCH -J illumina_merge

EXECDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"
python "${EXECDIR}/racon_merge.py" "$1" "$2" > "$3"
