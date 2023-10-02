#!/usr/bin/python3.6

import argparse, sys

parser=argparse.ArgumentParser()

parser.add_argument('--infile', help='Input file name (in blastn outfmt 7)')
parser.add_argument('--outfile', help='Name of output file to write')
parser.add_argument('--minidentity', help='Minimum accepted identity level of sequences (float)', type = float, default = 50.0)
parser.add_argument('--maxevalue', help='Maximum accepted e-value (float)', type = float, default = 1e-10)

args=parser.parse_args()

infile = args.infile
outfile = args.outfile
max_evalue = args.maxevalue
min_identity = args.minidentity

try:
    prev_id = ""
    prev_evalue = 1
    queries = []
    with open(outfile, 'w') as out:
        out.write("FILTERED BLAST RESULTS\n")
        out.write(f"Printing queries and unique gene names with\n")
        out.write(f"\t - minimum identity of {min_identity}\n")
        out.write(f"\t - maximum e-value of {max_evalue}\n\n")

    with open(infile) as f:
        for row in f:
            if row.startswith( "# Query" ) or "hits found" in row:
                with open(outfile, 'a') as out:
                    out.write(row)
            elif not row.startswith("#"):
                array = row.split("\t")
                query = array[0]
                current_id = array[1]
                evalue = float(array[10])
                identity = float(array[2])
                if current_id != prev_id and evalue <= max_evalue and identity >= min_identity and query not in queries:
                    queries.append(query)
                    first_evalue = evalue
                    with open(outfile, 'a') as out:
                        out.write(f"{current_id}\t{identity}\t{evalue}\n")
                elif current_id != prev_id and evalue <= max_evalue and identity >= min_identity and query in queries and evalue <= first_evalue:
                    with open(outfile, 'a') as out:
                        out.write(f"{current_id}\t{identity}\t{evalue}\n")
                prev_id = current_id
except TypeError:
    print("Missing arguments! For help please use the -h option")
