from ast import parse
from heapq import merge
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
parser.add_argument('--output_dir', type=str, default='Results', help='output directory name. Default is Results')
parser.add_argument('-t', type=int, default=20, help='number of threads used')

args = parser.parse_args()

barcodes = pd.read_csv(args.barcode_kit_file, sep=',')

barcode_ids = {}

for i in barcodes.index:
    barcode_ids[(barcodes.i7[i], barcodes.i5[i])] = barcodes.id[i]


print(barcodes)

try:
    os.mkdir('{}'.format(args.output_dir))
    os.mkdir('{}/barcodes_db'.format(args.output_dir))
    os.mkdir('{}/results'.format(args.output_dir))
except:
    pass

# create db

with open('{}/barcodes_db/adapters_i7.fasta'.format(args.output_dir), 'w') as f:
    read_cnt = 0
    for i in barcodes.index:
        
        f.write('>{}\n'.format(barcodes.i7[i]))
        f.write(
            adapter_l.replace('[i7]', reverse_complement(barcodes.i7[i])))
        f.write('\n')
        

with open('{}/barcodes_db/adapters_i5.fasta'.format(args.output_dir), 'w') as f:
    read_cnt = 0
    for i in barcodes.index:
        
        f.write('>{}\n'.format(barcodes.i5[i]))
        f.write(
            adapter_r.replace('[i5]', str(barcodes.i5[i]))
        )
        f.write('\n')

print('-----------------Indexing-----------------')

# indexing
subprocess.run('bwa index {}/barcodes_db/adapters_i7.fasta'.format(args.output_dir), shell=True)
subprocess.run('bwa index {}/barcodes_db/adapters_i5.fasta'.format(args.output_dir), shell=True)

print()
print('-----------------Mapping-----------------')

# mapping
subprocess.run('bwa mem -t {t} {fol}/barcodes_db/adapters_i7.fasta {reads} > {fol}/barcodes_i7_full.sam'.format(t = args.t, fol = args.output_dir, reads=args.reads), shell=True)
subprocess.run('bwa mem -t {t} {fol}/barcodes_db/adapters_i5.fasta {reads} > {fol}/barcodes_i5_full.sam'.format(t = args.t, fol = args.output_dir, reads=args.reads), shell=True)

awk_command = "'" + "{" + "print $1{tab}$2{tab}$3{tab}$4{tab}$5{tab}$6".format(tab = '"' + "\t" + '"') + "}" + "'"
subprocess.run('cat {fol}/barcodes_i7_full.sam | awk -F "\t" {awk_command} > {fol}/barcodes_i7.sam'.format(fol = args.output_dir, awk_command = awk_command), shell=True)
subprocess.run('cat {fol}/barcodes_i5_full.sam | awk -F "\t" {awk_command} > {fol}/barcodes_i5.sam'.format(fol = args.output_dir, awk_command = awk_command), shell=True)
#subprocess.run('bwa mem -t {t} {fol}/barcodes_db/adapters_i5.fasta {reads} > {fol}/barcodes_i5.sam'.format(t = args.t, fol = args.output_dir, reads=args.reads), shell=True)


print()
print('-----------------SAM files processing-----------------')

sam_i7 = pd.read_csv(
    '{}/barcodes_i7.sam'.format(args.output_dir), 
    sep='\t', 
    comment='@', 
    header=None, 
    usecols=range(6),
    names=['read', 'flag', 'refseq', 'pos', 'mapq', 'cigar']
)

sam_i7 = sam_i7[sam_i7.flag.isin([0,16])]


sam_i5 = pd.read_csv(
    '{}/barcodes_i5.sam'.format(args.output_dir), 
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



print()
print('-----------------Getting barcode variants-----------------')
variants = [(barcodes.i7[i], barcodes.i5[i]) for i in barcodes.index]

matches = 0
mismatches = 0
drop_index = []

print()
print('-----------------Filtering barcode variants-----------------')
for i in merged_df.index: 
    if (merged_df.refseq_x[i], merged_df.refseq_y[i]) in variants:
        matches += 1
    else:
        drop_index.append(i)
merged_df = merged_df[~merged_df.index.isin(drop_index)]

merged_df.reset_index(drop=True, inplace=True)

# save trimmed reads as one file


print()
print('-----------------Reads filtering-----------------')

f = open('{}/passed_reads.fastq'.format(args.output_dir), 'w')

fin = parse(args.reads, format='fastq')


passed_reads = set(merged_df.read)

overhang = 20


vls_data = {}


for df_i in merged_df.index:

    _vls = list(merged_df.iloc[df_i])
    vls_data[merged_df.read[df_i]] = _vls




for rec in fin:
    desc = rec.description.split(' ')[0]
    if desc not in passed_reads:
        continue
    
    vls = vls_data[desc]
    
    cigar1, cigar2 = vls[5], vls[10]
    
    try:
        len_l = int(cigar1.split('M')[-1][:-1])
        len_r = int(cigar2.split('M')[-1][:-1])

    except:
        print(desc, cigar1, cigar2)
        continue
        
    if vls[1] == 16:
        len_l, len_r = len_r, len_l
    
    f.write(rec[-len_l+overhang:-(len(rec.seq) - len_r)-overhang].format('fastq'))
    
f.close()


print()
print('-----------------Reads demultiplexing-----------------')
# demultiplex reads
trimmed_reads = parse('{}/passed_reads.fastq'.format(args.output_dir), format='fastq')

files = {}

for rec in trimmed_reads:
    desc = rec.description.split(' ')[0]
    
    vls = vls_data[desc]
    
    if (vls[2], vls[7]) not in files:
        files[(vls[2], vls[7])] = open('{}/results/{}.fastq'.format(args.output_dir, barcode_ids[(vls[2], vls[7])]), 'w')
    
    files[(vls[2], vls[7])].write(
        rec.format('fastq')
    )

for pair in files:
    files[pair].close()

print('Done!')
