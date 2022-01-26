#!/usr/bin/env python3

import sys
from Bio import SeqIO

def main():
    # ensure we can write outputs and that all inputs are available
    with open(sys.argv[1], 'w') as f:
        pass
    with open(sys.argv[2], 'w') as f:
        pass
    THRESHOLD = int(sys.argv[3])

    # compute coverage using data from stdin
    contig_reads = {}
    for line in sys.stdin:
        contig_name, _, _, read_count, _, _ = line.strip().split('\t')
        reads = contig_reads.get(contig_name, 0)
        contig_reads[contig_name] = reads + int(read_count)

    # accept contigs if they have average read count >= threshold
    # reject otherwise
    with open(sys.argv[1], 'w') as accept:
        with open(sys.argv[2], 'w') as reject:
            with open('assembly.fasta') as f:
                accept_batch = []
                reject_batch = []
                for record in SeqIO.parse(f, 'fasta'):
                    contig_name = record.id
                    reads = contig_reads.get(contig_name, None)
                    if reads is None:
                        print(f'Contig {contig_name} is in assembly.fasta but there is no read count entry for it', file=sys.stderr)
                        sys.exit(1)

                    average_reads = reads / len(record.seq)
                    if average_reads >= THRESHOLD:
                        print(f'Accepted {contig_name} because it has an average of {average_reads} reads')
                        accept_batch.append(record)
                        if len(accept_batch) >= 100:
                            SeqIO.write(accept_batch, accept, 'fasta')
                            accept_batch = []
                    else:
                        print(f'Rejected {contig_name} because it has an average of {average_reads} reads')
                        reject_batch.append(record)
                        if len(reject_batch) >= 100:
                            SeqIO.write(reject_batch, reject, 'fasta')
                            reject_batch = []
                if accept_batch:
                    SeqIO.write(accept_batch, accept, 'fasta')
                if reject_batch:
                    SeqIO.write(reject_batch, reject, 'fasta')

if __name__ == '__main__':
    main()
