import sys


def compute_gc(seq):
    """
    Computes the GC content of a sequence

        Parameters:
            seq (str): A DNA sequence

        Returns:
            gc_content (float): The frequency of GCs in the sequence
    """

    # Initialize variable to store GC frequency
    gc_count = 0

    # Count number of GC
    for base in seq:
        if base in ['G', 'C']:
            gc_count += 1

    # Normalize
    gc_content = gc_count/len(seq)

    return gc_content


def get_sequences(filename):
    """
    Returns the sequences in a fasta file

        Parameters:
            filename (str): Name of the fasta file

        Returns:
            sequences (dict): Dictionary of sequences keyed by fasta record
    """

    # Initialize empty dict to store sequences
    sequences = {}

    # Read fasta file and store in dictionary
    with open(filename, 'r') as fasta:
        record_name = ''
        for line in fasta:
            if line.startswith('>'):
                record_name = line.strip()
                sequences[record_name] = ''
            else:
                sequence = line.strip()
                sequences[record_name] += sequence

    return sequences


def main():
    """Prints out GC-contents of all sequences in fasta file"""
    # Get the filename from the command line argument
    filename = sys.argv[1]

    sequences = get_sequences(filename)

    # Initialize dictionary to store GC contents
    gc_contents = {}

    # Iterate through sequences and compute GC content for each, then store in
    # dictionary
    for record_name in sequences:
        sequence = sequences[record_name]
        gc_content = compute_gc(sequence)
        gc_contents[record_name] = gc_content

    print(gc_contents)


if __name__ == "__main__":
    main()
