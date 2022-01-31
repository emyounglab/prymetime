#!/bin/bash

#fail if there's a typo in variable names
set -u
#fail if any command fails
set -e

MIN_COVERAGE_DEPTH=10

log() {
    echo "[$( date '+%Y-%m-%d %H:%M:%S' ) $( basename ${BASH_SOURCE[0]} )]: $1"
}

EXECDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"
cd "$4"

log "Building index"
bowtie2-build --threads ${N_THREADS} "$1" refgenome
echo "Done"

log "Mapping reads"
bowtie2 -x refgenome --threads ${N_THREADS} --no-unal -1 "$2" -2 "$3" -S - \
    | samtools view --threads ${N_THREADS} -b - \
    | samtools sort --threads ${N_THREADS} -m 5G - -o mapping_result_sorted.bam
log "Done"

log "Building index"
samtools index -@ ${N_THREADS} mapping_result_sorted.bam
log "Done"
log "Filtering contigs"
samtools mpileup mapping_result_sorted.bam | python3 "${EXECDIR}/filter_contigs.py" accept.fasta reject.fasta ${MIN_COVERAGE_DEPTH}
log "Done"

rm -f refgenome*.bt2
