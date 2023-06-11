#!/usr/bin/env python3
import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from myutils import compute_pfm, compute_pwm, ScoreSeq, extract_motifs, GetThreshold, RandomSequence, compute_pvalue, compute_nucleotide_frequencies
import pandas as pd
import logomaker
import datetime
import time
import numpy as np
import os

def main():
    print(f"Current working directory: {os.getcwd()}")
    print(sys.argv)

    parser = argparse.ArgumentParser(
        prog="myhomer",
        description="Command-line script to perform motif finding against set of known motifs"
    )

    parser.add_argument("peaks_file", help="HOMER Peak file", type=str, metavar="FILE")
    
    parser.add_argument("fasta_ref", \
                        help="faidx indexed reference sequence file", \
                        type=str, metavar="FILE")

    parser.add_argument("output_dir", help="Output directory", type=str)
    
    args = parser.parse_args()

    if not os.path.exists(args.fasta_ref):
        raise FileNotFoundError(f"{args.fasta_ref} does not exist")

    fasta_sequences = SeqIO.to_dict(SeqIO.parse(args.fasta_ref, 'fasta'))

    # Print variable settings
    peaks_file = args.peaks_file
    fasta_ref = args.fasta_ref
    output_dir = args.output_dir
    print("\n====================Variables setting====================")
    print(f"\nHOMER Peaks File: {peaks_file}")
    print(f"\nReference Genome: {fasta_ref}")
    print(f"\nOutput Directory: {output_dir}")
    print("=========================================================\n")

    # Read HOMER peaks file, ignoring the first 39 lines and treating line 40 as column headers
    peaks = pd.read_csv(peaks_file, sep='\t', skiprows=39, header=0)

    # MEME-formatted motif database
    motif_db = "../motif_databases/MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme"
    motifs = extract_motifs(motif_db)

    # -------------------------------------------------------------
    #
    # Run
    #
    # -------------------------------------------------------------

    # Extract 100bp sequences centered on the midpoint of each peak
    sequences = []
    for index, row in peaks.iterrows():
        seq_id = row['chr']
        start = row['start']
        end = row['end']
        mid_point = start + ((end - start) // 2)
        start_seq = max(mid_point - 50, 0)
        end_seq = min(mid_point + 50, len(fasta_sequences[seq_id].seq))
        sequences.append(str(fasta_sequences[seq_id].seq[start_seq:end_seq]))

    output_filename = "../benchmark/Klf4/sea_results/sequences.txt"

    # Write the sequences to the file
    with open(output_filename, 'w') as f:
        for sequence in sequences:
            f.write("%s\n" % sequence)

    nucleotide_frequencies = compute_nucleotide_frequencies(sequences)
    
    # Compute PFM, PWM
    pfm = compute_pfm(sequences)
    pwm = compute_pwm(pfm, len(sequences))

    # Generate a null distribution of scores
    numsim = 10000
    n = len(next(iter(motifs.values()))[0])
    null_scores = [ScoreSeq(pwm, RandomSequence(n, nucleotide_frequencies)) for i in range(numsim)]

    # Calculate score threshold
    pval_threshold = 0.01  # You can adjust this based on your desired level of significance
    thresh = GetThreshold(null_scores, pval_threshold)

    # Find the top motifs and append them to the motifs_sorted list
    motifs_scores = []
    for motif_name, motif_pwm in motifs.items():
        # Transpose the motif_pwm
        motif_pwm = np.transpose(motif_pwm)

        scores = [ScoreSeq(motif_pwm, seq) for seq in sequences]
        avg_score = np.mean(scores)
        pval = compute_pvalue(scores, null_scores)

        motifs_scores.append((motif_name, avg_score, pval))

    # Sort the motifs by score in descending order
    motifs_scores.sort(key=lambda x: x[1], reverse=True)

    # Get the top 10 motifs
    top_motifs = motifs_scores[:10]

    # Create output directory if it does not exist
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(os.path.join(output_dir, "motifs_enrichment_results.txt"), "w") as pwm_file:
        # Print the top 10 motifs to the file
        for i, (motif_name, avg_score, pval) in enumerate(top_motifs):
            pwm_file.write(f"{motif_name}, Score: {avg_score}, P-value: {pval}\n")

            # Get the PWM for the current motif
            motif_pwm = motifs[motif_name]

            # Transpose the motif_pwm
            motif_pwm = np.transpose(motif_pwm)

            # Generate and save sequence logo
            pwm_df = pd.DataFrame(motif_pwm, columns=['A', 'C', 'G', 'T'])
            logo = logomaker.Logo(pwm_df)
            logo.fig.savefig(os.path.join(output_dir, f"{motif_name}_logo.png"))

            # Write motif_pwm to the file
            np.savetxt(pwm_file, motif_pwm, delimiter='\t', fmt='%.3f')

            pwm_file.write('\n')  # Add a newline between motifs

if __name__ == "__main__":
    print("\nRunning myhomer ..... ")
    current_datetime = datetime.datetime.now()
    start_time = time.time()

    main()

    end_time = time.time()
    elapsed_time = end_time - start_time
    formatted_datetime = current_datetime.strftime("%Y-%m-%d %H:%M:%S")
    formatted_date, formatted_time = formatted_datetime.split(" ")

    print("\nDone!!!")
    print("\nDate:", formatted_date)
    print("Time:", formatted_time)
    m, s = divmod(elapsed_time, 60)
    h, m = divmod(m, 60)
    print(f"Elapsed Time: {int(h)}h:{int(m)}m:{s:.2f}s")
