# alphabets
def get_alphabet(alph: str = "DNA"):
    """Returns the alphabet as a list."""
    DNA_alphabet = ["A", "C", "G", "T"]
    RNA_alphabet = ["A", "C", "G", "U"]

    if alph == "DNA":
        return DNA_alphabet
    elif alph == "RNA":
        return RNA_alphabet
    else:
        raise ValueError("Invalid alphabet. Please choose 'DNA' or 'RNA'.")

def reverse_complement(seq: str, alph: str = "RNA"):
    """Expects a DNA or RNA sequence.
    Returns the reverse complement of the sequence."""
    if alph == "DNA":
        complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
    elif alph == "RNA":
        complements = {"A": "U", "C": "G", "G": "C", "U": "A"}
    else:
        raise ValueError("Invalid alphabet. Please choose 'DNA' or 'RNA'.")

    seq = seq[::-1]
    seq = "".join([complements[c] for c in seq])

    return seq

def complement(seq: str, alph: str = "RNA"):
    """Expects a DNA or RNA sequence.
    Returns the complement (not the reverse complement!) of the sequence."""
    if alph == "DNA":
        complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
    elif alph == "RNA":
        complements = {"A": "U", "C": "G", "G": "C", "U": "A"}
    else:
        raise ValueError("Invalid alphabet. Please choose 'DNA' or 'RNA'.")

    seq = "".join([complements[c] for c in seq])

    return seq

def one_point_mutation(seq: str, alph: str = "RNA"):
    """Expects a DNA or RNA sequence.
    Returns a generator object that yields all possible one point mutations."""
    nucleotides = get_alphabet(alph)
    for index, base in enumerate(seq):
        for nt in nucleotides:
            if base != nt:
                yield seq[:index] + nt + seq[index+1:]
                
                
def GC_content(seq: str):
    """Expects a DNA or RNA sequence.
    Returns the GC content of the sequence as a percentage."""
    seq = seq.upper()
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)