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
import os
import sys
import glob
import argparse

# @Trey: initially I tend to skip the "if __name__ == '__main__':" Python
# idiom. This makes it easier to debug since the state of the program is
# availble in the global namespace after execution in IPython.

#########################################
# Main
#########################################

# Parse input
parser = argparse.ArgumentParser(description='Ribo-Seq ORF detection.')
parser.add_argument('-g', '--gff', dest='annotations', required=True,
                   help='Location of GFF annotation file to use.')
parser.add_argument('forward-strand', metavar='PATH',
                   help=('Wildcard string specifying path to forward-strand'
                         'input files'))
parser.add_argument('reverse-strand', metavar='PATH',
                   help=('Wildcard string specifying path to reverse-strand'
                         'input files'))
args = parser.parse_args()

# Check for incorrect or missing input
if ((not os.path.isfile(args.annotations)) or
    (len(args.forward_strand) == 0) or
    (len(args.reverse_strand) == 0)):
        print("Invalid input!\n")
        parser.print_help()
        sys.exit()

# Load annotations

# Iterate over ribosome profiling coordinates
# Positive strand
for infile in args.forward_strand:
    if infile.endswith('.gz'):
        # opn
        pass
