import argparse

def extract_kmers(seq, k):
    """Decomposes sequence into a list of k-mers"""
    kmers = []
    for start_pos in range(len(seq)):
        end_pos = start_pos + k
        if end_pos <= len(seq):
            kmer = seq[start_pos:end_pos]
            kmers.append(kmer)

    return kmers


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-k", help="the k-mer length", type=int, default=3)
    parser.add_argument("seq", help="the sequence to be decomposed", type=str)

    args = parser.parse_args()

    kmers = extract_kmers(args.seq, args.k)
    print(kmers)
