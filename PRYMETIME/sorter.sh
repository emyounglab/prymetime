#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o sorter_%j.out
#SBATCH -e sorter_%j.err
#SBATCH -J sorter_merge

EXECDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"
cp "${EXECDIR}/sorter.py" "$1"
cd "$1"
python3 sorter.py
