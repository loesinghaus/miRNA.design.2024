regions = [(0, 1), (1, 8), (8, 11), (11, 14), (14, 17), (17, 20), (20, 21)]

def get_mismatches(target, query):
    """Returns a list of mismatches and wobbles between two sequences.
    Assumes that both target and query are miRNA sequences of the same length.
    DOES NOT WORK FOR AN ACTUAL TARGET SEQUENCE, WHICH IS THE REVERSE COMPLEMENT."""
    mismatch_positions = [0 for i in range(len(target))]
    wobble_positions = [0 for i in range(len(target))]
    for i in range(len(target)):
        if target[i] != query[i]:
            mismatch_positions[i] = 1
        if target[i] == 'A' and query[i] == 'G':
            wobble_positions[i] = 1
            mismatch_positions[i] = 0
        elif target[i] == 'C' and query[i] == 'T':
            wobble_positions[i] = 1
            mismatch_positions[i] = 0

    # the RISC likes an "A" in the target sequence at position 1
    if target[0] == 'U':
        mismatch_positions[0] = 0
        wobble_positions[0] = 0

    return mismatch_positions, wobble_positions

def sum_mismatches_in_regions(mismatch, wobble, regions=regions):
    """Returns a list of mismatches and wobbles in each region."""
    mismatch_counts = []
    wobble_counts = []
    for region in regions:
        mismatch_counts.append(sum(mismatch[region[0]:region[1]]))
        wobble_counts.append(sum(wobble[region[0]:region[1]]))

    try:
        assert sum(mismatch_counts) + sum(wobble_counts) == sum(mismatch) + sum(wobble)
    except AssertionError:
        print(mismatch_counts,wobble_counts,mismatch,wobble)

    return mismatch_counts, wobble_counts

def count_mismatches_in_region(target, query, regions=regions):
    """Assumes that both target and query are miRNA sequences of the same length,
    NOT the reverse complement.
    Returns a list of mismatches and wobbles in each region."""
    mismatch, wobble = get_mismatches(target, query)
    mismatch_counts, wobble_counts = sum_mismatches_in_regions(mismatch, wobble, regions)

    return mismatch_counts, wobble_counts


def get_mismatches_with_reverse_complement(sequence, target):
    """Returns a list of mismatches and wobbles between two sequences.
    Assumes that sequence is a miRNA and target an actual target sequence,
    e.g. the reverse complement of a miRNA sequence."""
    mismatch_positions = [1 for i in range(len(target))]
    wobble_positions = [0 for i in range(len(target))]

    complements = {"A": "U", "C": "G", "G": "C", "U": "A"}
    wobbles = {"A": "none", "C": "none", "G": "U", "U": "G"}
    
    # make the target the complement
    target = target[::-1]

    for i in range(len(target)):
        if complements[target[i]] == sequence[i]:
            mismatch_positions[i] = 0
        if wobbles[target[i]] == sequence[i]:
            wobble_positions[i] = 1
            mismatch_positions[i] = 0

    return mismatch_positions, wobble_positions