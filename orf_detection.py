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
"""
import io
import os
import re
import csv
import sys
import gzip
import glob
import argparse
from Bio import Seq, SeqIO

def main():
    """Main"""
    args = parse_args()

    # Load genome sequence
    chromosomes = {x.id:x for x in SeqIO.parse(args.fasta, format='fasta')}

    # Load annotations

    # Iterate over ribosome profiling coordinates
    # Positive strand
    for infile in args.forward_strand:
        # determine chromosome number
        ch_num = re.match('.*(LmjF\.\d+).*', infile).groups()[0]
        ch = chromosomes[ch_num]

        result = find_matching_orfs(infile, args.outdir, ch, '+',
                                    args.window_size, args.min_size)

def find_matching_orfs(infile, outdir, ch, strand, window_size, min_size):
    """Parses the specified file of putative CDSs and return a list of
       the best matching ORFs for each translated region."""
    # load file
    fp = gzip.open(infile) if infile.endswith('.gz') else open(infile)

    reader = csv.DictReader(io.TextIOWrapper(fp),
                            fieldnames=['id', 'start', 'stop'])

    # Output filenamoe
    outfilename = (os.path.basename(infile).replace('.csv', '')
                                           .replace('.gz', '') + ".bed")
    outfile = os.path.join(outdir, outfilename)

    # Open CSV writer
    outfp = open(outfile, 'w')
    writer = csv.writer(outfp, delimiter='\t')

    # Skip header row
    _ = next(reader)

    # Iterate over putative coding sequences
    for i, cds in enumerate(reader):
        # Length of putative CDS
        cds_length = int(cds['stop']) - int(cds['start'])

        # Skip over very small translated regions
        #if cds_length < 30:
        #    continue

        # Determine region to look for overlapping orfs
        start = int(cds['start']) - window_size
        stop = int(cds['stop']) + window_size

        # Pull sequence for the range
        cds_padded = ch[start:stop].seq

        # Find all possible ORFs in the above region
        orfs = find_orfs(cds_padded, min_size)
        print("Checking CDS %d (%d possible ORFs)" % (i, len(orfs)))

        # If no ORFs found in specified region, continue
        if len(orfs) == 0:
            print("No ORFs found for current CDS (len: %d)" % cds_length)
            continue

        # Find the best match
        best_match = orfs[0]
        best_score = compute_score(cds_padded, orfs[0], window_size)

        if len(orfs) > 1:
            for orf in orfs[1:]:
                score = compute_score(cds_padded, orf, window_size)
            if score > best_score:
                best_match = orf
                best_score = score

        print("Score: %f (len: %d)" % (best_score, cds_length))

        # Write entry in bed file
        orf_start = start + best_orf[0]
        orf_stop = start + best_orf[1]
        row = [ch.id, orf_start, orf_stop, "CDS%d" % i,
               best_score * 1000, best_orf[2]]
        writer.writerow(row)

    outfp.close()

    # Filter out known CDSs
    # Filter out low-scoring matches?

def compute_score(cds, orf, window_size):
    """
    Computes the match score of a putative ORF for a putative CDS.

    The goal of this function is to determine which of the ORFs contained
    within a given sequence range is most likely to correspond to the "true"
    CDS indicated by a translated region in the ribosome profiling data.

    To two criteria used to determine this are:
        1. Similarity between ribosome profiling region and ORF size
        2. Amount of overlap between the two regions

    Thus, a likely candidate ORF will be very close in size to the translated
    region, and will share a high overlap with the putative CDS.

    Arguments
    ---------
    cds : Bio.Seq
        A BioPython Seq instance representing the padded region containing a
        putative CDS.
    orf : tuple
        A triplet containing the start, stop, and strand information for an
        ORF in the region of interest.
    window_size: int
        The amount of bases added to either side of the CDS when searching for
        ORFs.

    Returns
    -------
    score : float
        Returns a normalized score in the range of 0-1.
    """
    START = 0
    END = 1

    # CDS and ORF lengths
    cds_length = len(cds) - (2 * window_size)
    orf_length =  orf[END] - orf[START]

    # Number of overlapping bases
    overlap = (min(orf[END], cds_length + window_size) -
               max(window_size, orf[START]))

    # Overlap score
    overlap_score = overlap / orf_length

    # Number of non-overlapping bases
    #non_overlap = orf_length - overlap

    # How similar ORF and CDS sizes are (1 = same size)
    fit_score = 1 - (abs(orf_length - cds_length) /
                     max(orf_length, cds_length))

    # QUESTION:
    # Should we penalize ORFs that are smaller than the CDSs?
    # Currently, an ORF that falls entirely within the CDS will achieve the
    # maximum overlap score (0.5).

    # For now, we will weight each component equally
    return (0.5 * fit_score) + (0.5 * overlap_score)

def find_orfs(seq, min_protein_length, strand=1, trans_table=1):
    """
    Finds ORFs of a specified minimum length in a SeqRecord.

    Based on: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec360
    """
    answer = []
    seq_len = len(seq)

    # nucleotide sequence
    nuc = seq if strand == 1 else seq.reverse_complement()

    # start at the begining of the translated sequence
    aa_start = 0

    # iterate over reading frames
    for frame in range(3):
        # translate sequence
        trans = str(nuc[frame:].translate(trans_table))
        trans_len = len(trans)

        # find the next methionine
        aa_start = trans.find("M", aa_start)
        ## The only thing I could see was the scope of aa_end looked
        ## odd.  So I moved it here.
        ## Other than that, I don't bother keeping the list of answers
        ## and sorting it, instead I simply print a bed entry.
        aa_end = 0

        while aa_start < trans_len:
            # find nearest stop codon
            aa_end = trans.find("*", aa_start)

            # if no stop codon is found, go until end of sequence
            if aa_end == -1:
                aa_end = trans_len

            # check to make sure that ORF meets size requirements
            if aa_end - aa_start >= min_protein_length:
                if strand == 1:
                    start = frame + aa_start * 3
                    end = min(seq_len, frame + (aa_end * 3) + 3)
                else:
                    start = seq_len - frame - (aa_end * 3) - 3
                    end = seq_len - frame - (aa_start * 3)

                str_strand = "+" if strand == 1 else '-'
                if (start > 0 and end < seq_len):
                    answer.append((start, end, str_strand))

            # Find next Methionine and continue search
            aa_start = trans.find("M", aa_start + 1)

            if aa_start == -1:
                break

    # sort and return list of ORFs
    answer.sort()
    return answer

def parse_args():
    """Parses input arguments"""
    # Parse input
    parser = argparse.ArgumentParser(description='Ribo-Seq ORF detection.')
    parser.add_argument('-f', '--fasta', required=True,
                    help='Location of genome.')
    parser.add_argument('-g', '--gff', required=True,
                    help='Location of GFF annotation file to use.')
    parser.add_argument('-o', '--outdir', default='output',
                    help='Directory to save output to.')
    parser.add_argument('-w', '--window-size', default=100,
                    help=('Additional bases on either side of putative CDS to '
                            'include when looking for ORFs.'))
    parser.add_argument('-m', '--min-size', default=50,
                        help='Minimum ORF amino acid size to search for.')
    parser.add_argument('forward', metavar='PATH',
                    help=('Wildcard string specifying path to forward-strand'
                            'input files'))
    parser.add_argument('reverse', metavar='PATH',
                    help=('Wildcard string specifying path to reverse-strand'
                            'input files'))
    args = parser.parse_args()

    # Check for incorrect or missing input
    if not os.path.isfile(args.gff):
        raise IOError("Invalid GFF filepath")
    if not os.path.isfile(args.fasta):
        raise IOError("Invalid FASTA filepath")

    # Expand glob string
    args.forward_strand = glob.glob(args.forward)
    args.reverse_strand = glob.glob(args.reverse)

    if len(args.forward_strand) == 0 or len(args.reverse_strand) == 0:
        raise IOError("Invalid input expression for forward or reverse strand")

    # Create output directory
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    return args

if __name__ == "__main__":
    sys.exit(main())
