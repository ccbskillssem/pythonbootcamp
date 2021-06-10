def base_frequencies(seq):
    """Compute base frequencies of a sequence"""

    # Get the length of the sequence
    sequence_len = len(seq)

    # Initialize base frequencies
    base_frequencies = {
        'A': 0,
        'C': 0,
        'T': 0,
        'G': 0
    }

    # Count bases
    for base in seq:
        base_frequencies[base] += 1

    # Normalize count
    for base in base_frequencies:
        base_frequencies[base] = base_frequencies[base]/sequence_len

    return base_frequencies


def reverse_complement(seq):
    """Compute reverse complement of a sequence."""

    # Initialize dict of complements
    complements = {
        'A': 'T',
        'C': 'G',
        'T': 'A',
        'G': 'C',
        'Y': 'R',
        'R': 'Y',
        'W': 'S',
        'S': 'W'
    }

    # Initialize reverse complement
    rev_seq = ''

    # Loop through and populate list with reverse complement
    for base in reversed(seq.upper()):
        rev_seq += complements[base] 

    return rev_seq


if __name__ == "__main__":
    # Executed only if the script is run directly. Will not execute if imported.
    # Prints the name of the function in this script
    print("Stuff this script can do: ", reverse_complement.__doc__, 
          base_frequencies.__doc__)
