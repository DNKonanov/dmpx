from ast import parse
import pandas as pd
import numpy as np
from argparse import ArgumentParser
import os
import subprocess
from adapter_templates import adapter_r, adapter_l

from Bio.Seq import reverse_complement
from Bio.SeqIO import parse

parser = ArgumentParser()

parser.add_argument('--reads', type=str, help='multiplexed fastq file', required=True)
parser.add_argument('--barcode_kit_file', type=str, default='barcodes.tsv', help='barcodes table')

args = parser.parse_args()

barcodes = pd.read_csv(args.barcode_kit_file, sep=',')

try:
    os.mkdir('tmpfiles')
    os.mkdir('tmpfiles/barcodes_db')
    os.mkdir('tmpfiles/results')
except:
    pass

# create db

with open('tmpfiles/barcodes_db/adapters_i7.fasta', 'w') as f:
    read_cnt = 0
    for i in barcodes.index:
        
        f.write('>{}\n'.format(barcodes.i7[i]))
        f.write(
            adapter_l.replace('[i7]', reverse_complement(barcodes.i7[i])))
        f.write('\n')
        

with open('tmpfiles/barcodes_db/adapters_i5.fasta', 'w') as f:
    read_cnt = 0
    for i in barcodes.index:
        
        f.write('>{}\n'.format(barcodes.i5[i]))
        f.write(
            adapter_r.replace('[i5]', str(barcodes.i5[i]))
        )
        f.write('\n')

# indexing
subprocess.run('bwa index tmpfiles/barcodes_db/adapters_i7.fasta', shell=True)
subprocess.run('bwa index tmpfiles/barcodes_db/adapters_i5.fasta', shell=True)

# mapping
subprocess.run('bwa mem -t 20 tmpfiles/barcodes_db/adapters_i7.fasta {} > tmpfiles/barcodes_i7.sam'.format(args.reads), shell=True)
subprocess.run('bwa mem -t 20 tmpfiles/barcodes_db/adapters_i5.fasta {} > tmpfiles/barcodes_i5.sam'.format(args.reads), shell=True)


sam_i7 = pd.read_csv(
    'tmpfiles/barcodes_i7.sam', 
    sep='\t', 
    comment='@', 
    header=None, 
    usecols=range(6),
    names=['read', 'flag', 'refseq', 'pos', 'mapq', 'cigar']
)

sam_i7 = sam_i7[sam_i7.flag.isin([0,16])]


sam_i5 = pd.read_csv(
    'tmpfiles/barcodes_i5.sam', 
    sep='\t', 
    comment='@', 
    header=None, 
    usecols=range(6),
    names=['read', 'flag', 'refseq', 'pos', 'mapq', 'cigar']
)

sam_i5 = sam_i5[sam_i5.flag.isin([0,16])]

merged_df = sam_i7.merge(sam_i5, on='read')
merged_df.dropna(inplace=True)
merged_df.sort_values(by=['refseq_x', 'refseq_y'], inplace=True)



variants = [(barcodes.i7[i], barcodes.i5[i]) for i in barcodes.index]

matches = 0
mismatches = 0
drop_index = []
for i in merged_df.index: 
    if (merged_df.refseq_x[i], merged_df.refseq_y[i]) in variants:
        matches += 1
    else:
        drop_index.append(i)
merged_df = merged_df[~merged_df.index.isin(drop_index)]


# save trimmed reads as one file

f = open('tmpfiles/passed_reads.fastq', 'w')

fin = parse(args.reads, format='fastq')


passed_reads = set(merged_df.read)

overhang = 20

for rec in fin:
    desc = rec.description.split(' ')[0]
    if desc not in passed_reads:
        continue
    
    vls = merged_df[merged_df.read == desc].values[0]
    
    cigar1, cigar2 = vls[5], vls[10]
    
    len_l = int(cigar1.split('M')[-1][:-1])
    len_r = int(cigar2.split('M')[-1][:-1])
    
    if vls[1] == 16:
        len_l, len_r = len_r, len_l
    
    f.write(rec[-len_l+overhang:-(len(rec.seq) - len_r)-overhang].format('fastq'))
    
f.close()


# demultiplex reads
trimmed_reads = parse('tmpfiles/passed_reads.fastq', format='fastq')

files = {}

for rec in trimmed_reads:
    desc = rec.description.split(' ')[0]
    vls = merged_df[merged_df.read == desc].values[0]
    
    if (vls[2], vls[7]) not in files:
        files[(vls[2], vls[7])] = open('tmpfiles/results/{}_{}_barcode.fastq'.format(vls[2], vls[7]), 'w')
    
    files[(vls[2], vls[7])].write(
        rec.format('fastq')
    )

for pair in files:
    files[pair].close()