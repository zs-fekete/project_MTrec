#!/usr/bin/python3

import sys
import os
import argparse
import subprocess
import collections
import bamnostic as bs

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--infile", help = 'input bam file')
parser.add_argument('-o', '--outfile', help = 'textfile listing pure MT readnames',
                    default = 'out.txt')
parser.add_argument('-b', '--bamout', help = 'pure MT read bam output',
                    default = 'puremt.bam')
args = parser.parse_args()

#Read input bam
bam = bs.AlignmentFile(args.infile, 'rb')

#Create a list of not pure MT reads
def get_nonpure_ids(bam):
    """Find reads with supplementary alignments not on MT"""
    nonmt_list = []
    for read in bam:
        try:
            sa_tag = read.get_tag('SA')
            rname = read.query_name
    
            for tag in sa_tag.strip().split(';'):
                if not tag.startswith('MT') and len(tag) > 2:
                    nonmt_list.append(read.query_name)
                    print(tag.split(',')[0])

        except KeyError:
            pass
        
    return(list(set(nonmt_list)))

nonpure = get_nonpure_ids(bam)

bam = bs.AlignmentFile(args.infile, 'rb')
readidlist = [read.query_name for read in bam
              if read.query_name not in nonpure]

#Write pure MT read ids to file
with open(args.outfile, 'w') as out:
    for i in readidlist:
        out.write(f'{i}\n')

#Run samtools to get pure MT reads in a bam file
command = f'samtools view -b -o {args.bamout} -N {args.outfile} {args.infile}'
subprocess.call(command, shell = True)
