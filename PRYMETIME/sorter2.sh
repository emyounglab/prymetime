#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16

EXECDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"
cp "${EXECDIR}/sorter2.py" "$1"
cd "$1"
python3 sorter2.py
