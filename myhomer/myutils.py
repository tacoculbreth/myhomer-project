import numpy as np
import random
from scipy.stats import pearsonr, percentileofscore, norm
import sys
import re

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def ERROR(msg):
	"""
	Print an error message and die

	Parameters
	----------
	msg : str
	   Error message to print
	"""
	sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + "{msg}\n".format(msg=msg) )
	sys.exit(1)

def extract_motifs(pfm_file):
    """Extract motifs from MEME-formatted motifs database

    Parameters
    ----------
    pfm_file: MEME-formatted motifs database
        Set of known motifs from HOCOMOCO

    Returns
    -------
    all motifs from database 
    """
    motifs = {}

    with open(pfm_file, 'r') as f:
        lines = f.readlines()
        motif_name = ''
        motif_data = []
        for line in lines:
            if line.startswith('MOTIF'):
                # Save previous motif
                if motif_name:
                    # Transpose PFM
                    motifs[motif_name] = list(map(list, zip(*motif_data)))
                    motif_data = []

                # Extract new motif name
                motif_name = line.split()[1]

            # Save motif data
            if re.match(r'^\d', line):
                motif_data.append(list(map(float, line.strip().split())))

        # Save last motif
        if motif_name and motif_data:
            # Transpose the PFM
            motifs[motif_name] = list(map(list, zip(*motif_data)))  

    return motifs


def compute_pfm(sequences):
    """Compute Positional Frequency Matrix

    Parameters
    ----------
    sequences : list of str
        List of sequences (e.g. binding_sites)
    
    Returns
    -------
        pfm : 2d np.array
    """
    pfm = np.zeros((4, len(sequences[0])))
    base_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    for sequence in sequences:
        for i, base in enumerate(sequence):
            pfm[base_to_index[base]][i] += 1

    return pfm

def compute_pwm(pfm, num_sequences, pseudocount=0.001):
    """Compute Positional Weight Matrix

    Parameters
    ----------
    pfm: 2d-matrix
        Positional Frequency Matrix
    num_sequencies: int
        Number of sequences
    pseudocount: double
        Pseudocount to add to 0 to avoid taking log of 0

    Returns
    -------
        pwm : 2d np.array
    """
    # Normalize PFM to create a Position Probability Matrix (PPM)
    ppm = pfm / num_sequences

    # Add pseudocounts and compute PWM
    pwm = np.log2(np.where(ppm==0, pseudocount, ppm + pseudocount))

    return pwm

def ScoreSeq(pwm, sequence):
    """ Score a sequence using a PWM
    
    Parameters
    ----------
    pwm : 2d np.array
       Position weight matrix
    sequence : str
       Sequence of nucleotides to be scored
       
    Returns
    -------
    score : float
       PWM score of the sequence
    """
    pwm_width = pwm.shape[1]
    max_score = -np.inf
    base_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}  # Add this dictionary to map base to index

    for start in range(len(sequence) - pwm_width + 1):
        subseq = sequence[start:start+pwm_width]
        score = 0
        for i in range(pwm_width):
            score += pwm[base_to_index[subseq[i]], i]  # Use the dictionary to index pwm

        if score > max_score:
            max_score = score

    return max_score

def compute_nucleotide_frequencies(sequences):
    """Compute the frequencies of each nucleotide (A, C, G, T) in a list of sequences.

    Parameters
    ----------
    sequence : str
        Sequence of nucleotides
    
    Returns
    -------
    dictionary of nucleotide frequencies
    """
    count_A = count_C = count_G = count_T = 0
    total_length = 0

    for seq in sequences:
        count_A += seq.count('A')
        count_C += seq.count('C')
        count_G += seq.count('G')
        count_T += seq.count('T')
        total_length += len(seq)

    freq_A = count_A / total_length
    freq_C = count_C / total_length
    freq_G = count_G / total_length
    freq_T = count_T / total_length

    return [freq_A, freq_C, freq_G, freq_T]

def RandomSequence(length, nucleotide_frequencies):
    """ Generate random sequence from given frequencies

    Parameters
    ----------
    length: int
        Length of random sequence to generate
    nucleotide_frequencies: dictionary
        Dictonary of nucleotide frequencies for A, C, G, T
    
    Returns
    -------
    Random sequence with given nucleotide frequencies
    """
    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(np.random.choice(nucleotides, p=nucleotide_frequencies) for _ in range(length))


def GetThreshold(null_dist, pval):
    """ Find the threshold to achieve a desired p-value
    
    Given a null distribution (list of values),
    find the threshold to achieve a desired p-value
    
    Parameters
    ----------
    null_dist : list of float
       Null distribution of scores
    pval : float
       % of null_dist that should be above the threshold returned
       
    Returns
    -------
    thresh : float
       Threshold to achieve the desired p-value    
    """
    thresh = 0 
    null_dist_sorted = sorted(null_dist, reverse=True)
    idx = int(len(null_dist_sorted) * pval)
    thresh = null_dist_sorted[idx]
    
    return thresh

def compute_pvalue(scores, null_scores):
    """Calculate the p-value for the scores compared to the null distribution.
    
    Parameters
    ----------
    scores: list of scores
    null_scores: list of null scores

    Returns
    -------
    p-value
    """
    # Remove non-finite values from null_scores and scores
    null_scores = np.array(null_scores)
    null_scores = null_scores[np.isfinite(null_scores)]
    
    scores = np.array(scores)
    scores = scores[np.isfinite(scores)]

    # Fit a normal distribution to the null scores
    mu, std = norm.fit(null_scores)
    
    # Compute the p-value for each score and take the average
    pvals = norm.sf(scores, loc=mu, scale=std)  # sf = 1 - cdf
    pval = np.mean(pvals)
    
    return pval

def search_motifs(sequences, motifs, pwm_func=compute_pwm):
    """ Search sequences against a database of motifs

    Parameters
    ----------
    sequences : list of str
        The sequences to search against the motifs.
    motifs : dict
        The database of motifs, in the form of a dict mapping motif names to PFMs.
    pwm_func : callable
        Function to convert a PFM to a PWM. Default is `compute_pwm`.

    Returns
    -------
    dict
        A dict mapping motif names to average scores across all sequences and their enrichment p-values.
    """
    scores = {motif_name: 0 for motif_name in motifs.keys()}
    pvals = {motif_name: 0 for motif_name in motifs.keys()}
    
    for sequence in sequences:
        for motif_name, pfm in motifs.items():
            # Convert PFM (list of lists) to a NumPy array
            pfm_np = np.array(pfm)
            pwm = pwm_func(pfm_np, len(sequence))  # Convert PFM to PWM
            score = ScoreSeq(pwm, sequence)  # Score the sequence
            scores[motif_name] += score  # Add to the motif's total score
    
    # Normalize scores by number of sequences
    for motif_name in scores.keys():
        scores[motif_name] /= len(sequences)
        
    # Generate null distribution of scores for each motif and calculate p-values
    for motif_name in motifs.keys():
        pfm = motifs[motif_name]
        pfm_np = np.array(pfm)
        pwm = pwm_func(pfm_np, len(sequences[0]))  # Convert PFM to PWM
        null_scores = [ScoreSeq(pwm, RandomSequence(len(sequences[0]), nucleotide_frequencies)) for _ in range(10000)]  # Generate null distribution
        if np.std(null_scores) > 0:
            pvals[motif_name] = compute_pvalue([scores[motif_name]], null_scores)  # Calculate p-value
        else:
            pvals[motif_name] = np.nan

    
    # Sort motifs by p-value
    motifs_sorted = sorted(zip(scores.keys(), scores.values(), pvals.values()), key=lambda x: x[2])
    
    return motifs_sorted

def motif_similarity(pwm1, pwm2):
    """Calculate similarity between two PWMs using Pearson Correlation

    Parameters
    ----------
    pwm1 : numpy.ndarray
    pwm2 : numpy.ndarray

    Returns
    -------
    float
        Similarity score
    """
    # Check if PWMs have the same length
    if pwm1.shape[1] != pwm2.shape[1]:
        # If not, pad the shorter one with zeros
        if pwm1.shape[1] < pwm2.shape[1]:
            pwm1 = np.pad(pwm1, ((0,0), (0, pwm2.shape[1] - pwm1.shape[1])))
        else:
            pwm2 = np.pad(pwm2, ((0,0), (0, pwm1.shape[1] - pwm2.shape[1])))

    # Flatten the PWMs and compute the correlation
    return pearsonr(pwm1.flatten(), pwm2.flatten())[0]

def find_similar_motifs(query_pwm, motifs_db, n_top=10):
    """Find the most similar motifs to the query PWM

    Parameters
    ----------
    query_pwm : numpy.ndarray
        The query PWM.
    motifs_db : dict
        A dictionary of PFMs where keys are motif names and values are PFMs.
    n_top : int, optional
        The number of top similar motifs to return. The default is 10.

    Returns
    -------
    list
        A list of tuples where the first element is a motif name and the second is a similarity score.
    """
    # Convert PFMs to PWMs and score all motifs in the database
    scores = {name: motif_similarity(query_pwm, np.array(pwm)) for name, pwm in motifs_db.items()}

    # Sort the scores and take the top n
    top_motifs = sorted(scores.items(), key=lambda x: -x[1])[:n_top]

    return top_motifs
