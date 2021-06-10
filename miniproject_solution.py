"""Given a FASTA file, produces k-mer frequency histograms and k-mer indexfor 
all sequences in the file"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import design_process_example as handle_sequences
import argparse_example as handle_kmers
from modded_script_example import reverse_complement 


def extract_all_kmers(sequence, k):
    """
    Get all k-mers in sequence, forward and reverse complement

        Parameters:
            sequence (string): sequence to decompose into k-mers

        Returns:
            all_kmers (list): all k-mers in sequence (both forward and reverse 
                strand)
    """

    # Get k-mers of forward sequence
    all_kmers = handle_kmers.extract_kmers(sequence, k)

    # Add k-mers of reverse complement
    reverse_complement_seq = reverse_complement(sequence)
    all_kmers += handle_kmers.extract_kmers(reverse_complement_seq, k)

    return all_kmers


def get_unique_kmers(sequences, k):
    """
    Extracts all unique k-mers among sequences

        Parameters:
            sequences (dict): Sequences keyed by their FASTA record name
            k (int): k-mer length

        Returns:
            unique_kmers (list): the set of all unique k-mers in our sequences
    """

    # First extract all the k-mers in our sequences and put them into a giant 
    # list
    all_kmers = []
    for record_name in sequences:
        sequence = sequences[record_name]
        kmers = extract_all_kmers(sequence, k)
        all_kmers += kmers

    # Then turn the list into a set to get only the unique k-mers
    unique_kmers = set(all_kmers)

    return unique_kmers


def count_kmers(kmer_index, kmers):
    """
    Counts the occurence of every k-mer in a set of k-mers

        Parameters:
            kmer_index (list): the set of all unique k-mers in our sequences
            kmers (list): k-mers extracted from the sequence

        Returns:
            kmer_count (dict): Count of each k-mer keyed by the k-mer string 
            (example element: 'AAA': 1)
    """
    # Initialize empty dictionary keyed by k-mer index
    kmer_counts = {}
    for kmer in kmer_index:
        kmer_counts[kmer] = 0

    # Iterate through kmers in set and increment counts in dict
    for kmer in kmers:
        kmer_counts[kmer] += 1

    return kmer_counts


def kmer_count_hist(kmer_counts, filename):
    """
    Saves a histogram of the k-mer counts

        Parameters:
            kmer_counts (dict): Count of each k-mer keyed by the k-mer string 
                (example element: 'AAA': 1)
            filename (string): Path to output figure

        Returns: k-mer frequency histogram saved to file
    """
    # Get the k-mer length (we use this for the x-label but this is optional)
    kmer_length = len(list(kmer_counts.keys())[0])

    # Get the counts only which we'll use to make the histogram
    counts = list(kmer_counts.values())

    # Initialize figure object
    figure, axes = plt.subplots()

    # Make histogram
    axes.hist(counts)

    # Add some labels 
    plt.ylabel('Count')
    plt.xlabel('{}-mer Frequency'.format(kmer_length))

    # Save and close
    figure.savefig(filename, dpi=300)
    plt.close(figure)


def build_kmer_freq_index(all_kmer_counts):
    """
    Builds a k-mer frequency index for all sequences
        
        Parameters:
            all_kmer_counts (dict): Dictionaries of k-mer counts for every 
            sequence keyed by their sequence name

        Returns:
            kmer_freq_index (DataFrame): Normalized k-mer frequencies for each 
                sequence. k-mers are row labels; sequences are column labels 
    """
    kmer_freq_index = pd.DataFrame.from_dict(all_kmer_counts).sort_index()

    # Normalize by column sum
    kmer_freq_index = kmer_freq_index.div(kmer_freq_index.sum(axis=0), axis=1)

    return kmer_freq_index


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-k", help="the k-mer length", type=int, default=3)
    parser.add_argument("fasta", 
                        help="the path to the FASTA file with your sequences", 
                        type=str)

    args = parser.parse_args()

    sequences = handle_sequences.get_sequences(args.fasta)
    unique_kmers = get_unique_kmers(sequences, args.k)

    all_kmer_counts = {}
    for record_name in sequences:
        sequence = sequences[record_name]
        kmers = extract_all_kmers(sequence, args.k)
        kmer_counts = count_kmers(unique_kmers, kmers)

        sequence_name = record_name.split(' ')[0][1:]
        all_kmer_counts[sequence_name] = kmer_counts

        filename = '{}.png'.format(sequence_name)
        kmer_count_hist(kmer_counts, filename)

    kmer_freq_index = build_kmer_freq_index(all_kmer_counts)
    kmer_freq_index.to_csv('coronaviruses_freq_index.csv')


if __name__ == "__main__":
    main()


