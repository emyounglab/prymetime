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

### Make a dictionary of repetitive or linear from Flye output info file ###
##########################################################################
df = pd.read_csv("assembly_info.txt",sep='\t')

rep_D = {}

for x in range(0, len(df.index)):
    rep_D[df.loc[x,"#seq_name"]] = df.loc[x,"repeat"]

### Make a fasta file of only the circular and repetitive contigs ###
######################################################
fasta_sequences = SeqIO.parse(open("assembly.fasta"),'fasta')

cir_seqs = [x for x in fasta_sequences if circ_D[x.id] is "+" or rep_D[x.id] is "+"]

SeqIO.write(cir_seqs , "cir_contigs.fasta", "fasta")

### Make a fasta file of only the linear contigs ###
####################################################
fasta_sequences = SeqIO.parse(open("assembly.fasta"),'fasta')

lin_seqs = [x for x in fasta_sequences if circ_D[x.id] is '-' and rep_D[x.id] is '-']

SeqIO.write(lin_seqs , "lin_contigs.fasta", "fasta")


