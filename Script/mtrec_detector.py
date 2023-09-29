#!/usr/bin/python3

import sys, argparse, os
import collections
import bamnostic as bs
import seqhandlers as sq

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bamfile', help = "input bam file")
parser.add_argument('-v', '--vcffile', help = "input vcf file")
parser.add_argument('-p', '--out_positive',
                    help = 'filename for positive hit output')
parser.add_argument('--out_negative',
                    action = 'store_true',
                    help = 'output negative results, too')

args = parser.parse_args()

#Data to consider - turn into option later
opposite = ['14189', '14366']
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

vars = sq.get_variants(args.vcffile)

#Primitive detection of possible recombination
with open('negative.out', 'w') as neg:
    with open(args.out_positive, 'w') as out:
        for read in read_info:
            ref_snps = sq.snp_coordinates(read, ref_exclude)
            nonsnps = sq.nonsnp_coordinates(vars, read, ref_snps)
                
            readdic = sq.opposite_snps(opposite,
                                    sq.create_readdic(vars, read, ref_snps, nonsnps))
            ps_cigar = sq.read_summary(list(readdic.values()))
            res = sq.rec_possible(ps_cigar)
        
            if res == 'True':
                out.write(f'{read[0]}\t{read[1]}\t' + ''.join(ps_cigar) + f'\t{res}\n')                        
            elif args.out_negative:
                neg.write(f'{read[0]}\t{read[1]}\t' + ''.join(ps_cigar) + f'\t{res}\n')                        
