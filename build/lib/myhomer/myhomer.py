#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import os

def main():

    #----------------------------------------------------------------
    #
    # Input, Output, Arguments
    #
    #----------------------------------------------------------------
    
    parser = argparse.ArgumentParser(
        prog="myhomer",
        description="Command-line script to perform motif finding against set of known motifs"
    )

    # Input
    parser.add_argument("peaks_file", help="HOMER Peak file", type=str, metavar="FILE")
    
    parser.add_argument("fasta_ref", \
                        help="faidx indexed reference sequence file", \
                        type=str, metavar="FILE")

    # Output
    parser.add_argument("output_dir", help="Output directory", type=str)
    
    # Parse args
    args = parser.parse_args()

    #print("Current working directory:", os.getcwd())
    #print("Absolute path to peaks file:", os.path.abspath(args.peaks_file))

    # Load Fasta
    if not os.path.exists(args.fasta_ref):
        raise FileNotFoundError(f"{args.fasta_ref} does not exist")

    fasta_sequences = SeqIO.to_dict(SeqIO.parse(args.fasta_ref, 'fasta'))

    # Get the absolute path of the peaks file
    peaks_file = args.peaks_file

    # Read HOMER peaks file, ignoring the first 39 lines and treating line 40 as column headers
    peaks = pd.read_csv(peaks_file, sep='\t', skiprows=39, header=0)

    #----------------------------------------------------------------
    #
    # Extract Sequences from the start to end positions of each peak
    #
    #----------------------------------------------------------------

    sequences = []
    for index, row in peaks.iterrows():
        seq_id = row['chr']  # Now reading the 'chr' column
        start = row['start']
        end = row['end']
        sequences.append(str(fasta_sequences[seq_id].seq[start:end]))
    
    #----------------------------------------------------------------
    #
    # Create PFM
    #
    #----------------------------------------------------------------

    # Initialize position frequency matrix with zeros
    pfm = np.zeros((4, len(sequences[0])))

    # DNA base to index mapping
    base_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    for sequence in sequences:
        for i, base in enumerate(sequence):
            pfm[base_to_index[base]][i] += 1
    
    # Print PFM
    print("Position Frequency Matrix (PFM):")
    print(pfm)

    #----------------------------------------------------------------
    #
    # Create PWM
    #
    #----------------------------------------------------------------

    pwm = np.log2((pfm+1) / (len(sequences)+4))  # Add pseudocounts to avoid taking log of zero

    # Print PWM
    print("Position Weight Matrix (PWM):")
    print(pwm)

if __name__ == "__main__":
    main()