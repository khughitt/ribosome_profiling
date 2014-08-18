#!/usr/bin/env python3
"""
Ribosome profiling ORF detection
2014/07/31

Authors
-------
Trey Belew (abelew@gmail.com)
Keith Hughitt (khughitt@umd.edu)

Overview
--------
This script attempts to detect novel open reading frames (ORFs) using evidence
of translation from a ribosome profiling experiment, and a list of known coding
regions.

Currently, this script expects input data in the form of a collection of
(gzipped) CSV files for each strand and each chromosome of the target organism.
Each row in the CSV files includes an identifier, start, and stop position of a
translated region detected in the ribosome profiling data, e.g.:

    > head ch1_forward.csv.gz

    "","start","stop"
    "1",2745,3018
    "2",4413,4632
    "3",4794,5127
    ....

Usage
-----

./orf_detection.py --gff=annotations.gff \
        "/path/to/riboseq/forward-strand/data/*.forward.csv.gz"
        "/path/to/riboseq/reverse-strand/data/*.reverse.csv.gz"

TODO
----
- Everything...
"""
import io
import os
import csv
import sys
import gzip
import glob
import argparse
from Bio import Seq, SeqIO

# @Trey: initially I tend to skip the "if __name__ == '__main__':" Python
# idiom. This makes it easier to debug since the state of the program is
# availble in the global namespace after execution in IPython.

#########################################
# Main
#########################################

# Parse input
parser = argparse.ArgumentParser(description='Ribo-Seq ORF detection.')
parser.add_argument('-f', '--fasta', required=True,
                   help='Location of genome.')
parser.add_argument('-g', '--gff', required=True,
                   help='Location of GFF annotation file to use.')
parser.add_argument('forward', metavar='PATH',
                   help=('Wildcard string specifying path to forward-strand'
                         'input files'))
parser.add_argument('reverse', metavar='PATH',
                   help=('Wildcard string specifying path to reverse-strand'
                         'input files'))
args = parser.parse_args()

# Check for incorrect or missing input
if ((not os.path.isfile(args.gff)) or
    (len(args.forward) == 0) or
    (len(args.reverse) == 0)):
        print("Invalid input!\n")
        parser.print_help()
        sys.exit()

# Expand glob string
forward_strand = glob.glob(args.forward)
reverse_strand = glob.glob(args.reverse)

# Load genome sequence
chromosomes = {x.id:x for x in SeqIO.parse(args.fasta, format='fasta')}

# Load annotations

# Iterate over ribosome profiling coordinates
# Positive strand
for infile in forward_strand:
    # determine chromosome number
    ch_num = re.match('.*(LmjF\.\d+).*', infile).groups()[0]
    ch = chromosomes[ch_num]

    # load file
    if infile.endswith('.gz'):
        fp = gzip.open(infile)
    else:
        fp = open(infile)

    reader = csv.DictReader(io.TextIOWrapper(fp),
                            fieldnames=['id', 'start', 'stop'])

    # Iterate over putative coding sequences
    for i, region in enumerate(reader):
        # skip header row
        if region['id'] == '':
            continue

        # Determine region to look for overlapping orfs
        start = int(region['start']) - 100
        stop = int(region['stop']) + 100

        # Pull sequence for the range
        seq = ch[start:stop].seq

        # Find all possible ORFs in the above region
        orfs = find_orfs(seq, 50)

        # Find the best match
        best_match = orfs[0]

        # Score = len(orf) * percent_overlap
        best_score =

        if len(orfs) > 1:

            for orf in orfs[1:]


# Reverse strand

# Filter out known CDSs

# Generate a list of possible


def find_orfs(seq, min_protein_length, strand=1, trans_table=1):
    """
    Finds ORFs of a specified minimum length in a SeqRecord.

    Based on: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec360
    """
    answer = []
    seq_len = len(seq)

    # Get sequence associated with the specified location and strand
    if strand == 1:
        dna_seq = seq
    else:
        dna_seq = seq.reverse_complement()

    for frame in range(3):
        trans = str(dna_seq[frame:].translate(trans_table))
        trans_len = len(trans)
        aa_start = 0
        aa_end = 0

        # Iterate through ORFS in reading frame
        while aa_start < trans_len:
            # Set end counter to position of next stop codon
            aa_start = trans.find("M", aa_start)
            aa_end = trans.find("*", aa_start)

            # If no start or stop codons found, stop here
            if aa_start == -1 or aa_end == -1:
                break

            # extend stop codon until ORF is of sufficient length
            while (aa_end - aa_start < min_protein_length) and aa_end > -1:
                aa_end = trans.find("*", aa_end + 1)

            # If no ORFs of sufficent size found, stop here
            if aa_end == -1:
                break

            # Compute coordinates of ORF
            if strand == 1:
                start = frame + aa_start * 3
                end = min(seq_len, frame + aa_end * 3 + 3)
            else:
                start = seq_len - frame - aa_end * 3 - 3
                end = seq_len - frame - aa_start * 3

            # Add to output
            str_strand = "+" if strand == 1 else '-'
            answer.append((start, end, str_strand))

            # increment start counter and continue
            aa_start = aa_end + 1
    answer.sort()
    return answer

