#!/usr/bin/python3

import sys, argparse, os
import collections
import bamnostic as bs
import seqhandlers as sq
import ioprocessing as io

parser = argparse.ArgumentParser()

#Mandatory args
parser.add_argument('-b', '--bamfile', help = "input bam file")
parser.add_argument('-v', '--vcffile', help = "input vcf file")
parser.add_argument('-p', '--out_positive',
                    help = 'filename for positive hit output')

#Optional
parser.add_argument('--out_negative',
                    action = 'store_true',
                    help = 'output negative results, too')
parser.add_argument('--minstretch', type = int, default = 2,
                    help = 'Minimum length of type X to score positively')
parser.add_argument('--opposite', nargs = '*', type = str,
                    default = ['14189', '14366'],
                    help = 'SNPs called, but more prominent in default type')
parser.add_argument('--minscore', type = int, default = 2,
                    help = 'Minimum score to consider a read a positive hit')

parser.add_argument('--plotdata', action = 'store_true',
                    help = 'Output per site data for plotting')

args = parser.parse_args()

print('Reading data...')
#cigar tags to exclude to get reference coordinates
ref_exclude = [4, 1, 5, 6, 8]

#Get data from files
bamfile = bs.AlignmentFile(args.bamfile, "rb")
read_info = [ (read.query_name,
             read.reference_start,
             read.reference_end,
             #read.query_alignment_sequence,
             read.query_alignment_length,
             read.cigar)
             for read in bamfile ]

var = sq.get_variants(args.vcffile)

print('Input files:')
print(f'\t{args.bamfile}')
print(f'\t{args.vcffile}')
print(f'Output file:\n\t {args.out_positive}')
print(f'Switched SNPs: {args.opposite}')
print(f'Minimum length of type to score positively: {args.minstretch}')
print(f'Outputting negative results: {args.out_negative}')

print('Processing...')


fullreaddic = io.create_full_dic(read_info, var, args.opposite)
filtdic = io.filter_readdic(fullreaddic, args.minstretch, args.minscore)
delta = sq.supporting_reads(read_info)
#print(delta)

print('Writing output...')
with open(args.out_positive, 'w') as out:
    out.write(f'read\tcigar\tresult\tscore\ttype\n')
    for read, info in filtdic.items():
        if read in delta:
            out.write(read+'\t'+''.join(info[1])+'\t'+'\t'.join(info[2])+'\tdelta\n')
        else:
            out.write(read+'\t'+''.join(info[1])+'\t'+'\t'.join(info[2])+'\twt\n')
            
if args.out_negative:
    with open('negative.out', 'w') as neg:
        neg.write(f'read\tcigar\tresult\tscore\ttype\n')
        for read, readdic in fullreaddic.items():
            ps_cigar = sq.read_summary(list(readdic.values()))
            res = sq.rec_possible(ps_cigar, args.minstretch, args.minscore)

            if res[0] == 'False':
                if read in delta:
                    neg.write(f'{read}\t' + ''.join(ps_cigar) + '\t' + '\t'.join(res) + '\tdelta\n')
                else:
                    neg.write(f'{read}\t' + ''.join(ps_cigar) + '\t' + '\t'.join(res) + '\twt\n')

if args.plotdata:
    dat = io.plot_data(fullreaddic)
    with open('plotdata.tsv', 'w') as pld:
        for line in dat:
            pld.write("\t".join(str(i) for i in line))
            if line[0] in delta:
                pld.write('\tdelta\n')
            else:
                pld.write('\twt\n')


