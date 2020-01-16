# -*- coding: utf-8 -*-
"""
Title: Sending Contigs to Nucmer
Created on Tue Aug 13 2019

@author: Eric
@email: ericyoung7@gmail.com
"""
import glob, os
import pandas as pd
from Bio import SeqIO
from pymummer import nucmer

short_contigs = []
contigs = []

for x in SeqIO.parse(open("pilon.fasta"),'fasta'):

    if len(x.seq) < 50000:
        
        short_contigs.append(x)
        
        SeqIO.write(x, "%(x)s.fasta" % {'x':x.id}, 'fasta')
        
    else:

        contigs.append(x)
        #print("long", x.id)
        
for pathname in glob.glob("*.fasta"):
    
    basename = os.path.basename(pathname)
    
    for x in short_contigs:

        if x.id in basename :
        
            runner = nucmer.Runner(basename, basename, "%(x)s_out.coords" % {'x':x.id},
                                   maxmatch=True, simplify=False, mincluster=2000,
                                   min_id=99, min_length=2000, coords_header=True)  
            runner.run()

# The below lines are for saving fasta files of the contigs if desired
#SeqIO.write(short_contigs , "short_contigs.fasta", "fasta")
#SeqIO.write(lin_contigs , "lin_contigs.fasta", "fasta")

# The below lines are for visually checking which files are repetitive or not
'''
for pathname in glob.glob("*.coords"):

    basename = os.path.basename(pathname)
    name = basename.split(".")

    df = pd.read_csv(basename)

    print(df)

    if len(df.index) > 1 :

        print(name[0], "morethan 1")
'''

cir_rep_contigs = [x for x in SeqIO.parse(open("cir_contigs.fasta"), 'fasta')]

for x in short_contigs:
    
    if len(pd.read_csv("%(x)s_out.coords" % {'x': x.id}).index) > 4 :
        
        cir_rep_contigs.append(x)
        
    else:
        #print(x.id)
        contigs.append(x)


for x in cir_rep_contigs :
    
    SeqIO.write(cir_rep_contigs, "cir_rep_contigs.fasta", "fasta")
    
SeqIO.write(contigs, "polished_contigs.fasta", "fasta")
