# -*- coding: utf-8 -*-
"""
Title: Sorting Circular and Linear Contigs
Created on Tue Aug 13 2019

@author: Eric
@email: ericyoung7@gmail.com
"""
import glob, os
import pandas as pd
from Bio import SeqIO

### Make a dictionary of circular or linear from Flye output info file ###
##########################################################################
df = pd.read_csv("assembly_info.txt",sep='\t')

circ_D = {}

for x in range(0, len(df.index)):
    circ_D[df.loc[x,"#seq_name"]] = df.loc[x,"circ."]

### Make a fasta file of only the circular contigs ###
######################################################
fasta_sequences = SeqIO.parse(open("accept.fasta"),'fasta')

cir_seqs = [x for x in fasta_sequences if circ_D[x.id] is "Y"]

SeqIO.write(cir_seqs , "cir_contigs.fasta", "fasta")

### Make a fasta file of only the linear contigs ###
####################################################
fasta_sequences = SeqIO.parse(open("accept.fasta"),'fasta')

lin_seqs = [x for x in fasta_sequences if circ_D[x.id] is 'N']

SeqIO.write(lin_seqs , "lin_contigs.fasta", "fasta")
