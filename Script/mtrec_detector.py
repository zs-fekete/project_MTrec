#!/usr/bin/python3

import sys, argparse, os
import collections
import bamnostic as bs
import seqhandlers as sq

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

#options to add
#parser.add_argument('--return_score') - 

args = parser.parse_args()

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

print(args.opposite)
print(args.out_positive)
print(args.bamfile)
print(args.vcffile)
print(args.minstretch)
print(args.out_negative)

complete_reads = collections.defaultdict(lambda : [])

readnames = list(set([read[0] for read in read_info]))
fullreaddic = {readname : [] for readname in readnames}

for read in read_info:
    ref_snps = sq.snp_coordinates(read, ref_exclude)
    nonsnps = sq.nonsnp_coordinates(var, read, ref_snps)
        
    readdic = sq.opposite_snps(args.opposite,
                            sq.create_readdic(var, read, ref_snps, nonsnps))

    if len(fullreaddic[read[0]]) != 0:
        fullreaddic[read[0]] = fullreaddic[read[0]] | readdic
    else:
        fullreaddic[read[0]] = readdic
    #print(readdic)

with open(args.out_positive, 'w') as out:
    out.write(f'read\tcigar\tresult\n')
    with open('negative.out', 'w') as neg:
        for read, readdic in fullreaddic.items():
            ps_cigar = sq.read_summary(list(readdic.values()))
            res = sq.rec_possible(ps_cigar, args.minstretch)
            
            if res == 'True':
                out.write(f'{read}\t' + ''.join(ps_cigar) + f'\t{res}\n')
            else:
                neg.write(f'{read}\t' + ''.join(ps_cigar) + f'\t{res}\n')

