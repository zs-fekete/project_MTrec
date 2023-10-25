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
parser.add_argument('--out_negative', #momentarily doesn't do anything
                    action = 'store_true',
                    help = 'output negative results, too')
parser.add_argument('--minstretch', type = int, default = 2,
                    help = 'Minimum length of type X to score positively')
parser.add_argument('--opposite', nargs = '*', type = str,
                    default = ['14189', '14366'],
                    help = 'SNPs called, but more prominent in default type')
parser.add_argument('--minscore', type = int, default = 2,
                    help = 'Minimum score to consider a read a positive hit')

#options to add
#parser.add_argument('--return_score') - 

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

#complete_reads = collections.defaultdict(lambda : [])

fullreaddic = io.create_full_dic(read_info, var, args.opposite)

print('Writing output...')
with open(args.out_positive, 'w') as out:
    out.write(f'read\tcigar\tresult\n')
    with open('negative.out', 'w') as neg:
        for read, readdic in fullreaddic.items():
            ps_cigar = sq.read_summary(list(readdic.values()))
            res = sq.rec_possible(ps_cigar, args.minstretch, args.minscore)
            
            if res[0] == 'True':
                out.write(f'{read}\t' + ''.join(ps_cigar) + '\t' + '\t'.join(res) + '\n')
            else:
                neg.write(f'{read}\t' + ''.join(ps_cigar) + '\t' + '\t'.join(res) + '\n')

