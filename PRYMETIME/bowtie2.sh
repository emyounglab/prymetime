#!/bin/bash

log() {
    echo "[$( date '+%Y-%m-%d %H:%M:%S' ) $( basename ${BASH_SOURCE[0]} )]: $1"
}

cd "$4"

log "Building index"
bowtie2-build --threads ${N_THREADS} "$1" medaka-idx
log "Done"

log "Mapping reads"
bowtie2 -x medaka-idx --threads ${N_THREADS} --no-unal -U "$2" -S "$3"
log "Done"

log "Removing index"
rm medaka-idx*.bt2
