import numpy as np
regions = [(0, 1), (1, 8), (8, 11), (11, 14), (14, 17), (17, 20), (20, 21)]

def region_split(sequence, regions=regions):
    """Splits a sequence into a list of regions, each of which is a string.
    Assumes that the sequence is a miRNA sequence."""
    
    return [sequence[region[0]:region[1]] for region in regions]

def get_mismatches_two_mirnas(intended_target, query):
    """Returns a list of individual mismatches and wobbles between the target sequence of a microRNA intended_target
    and a second attacking microRNA query. This functions does not count an A in the first position of the
    actual target sequence as a mismatch regardless of the query sequence.
    
    Assumes that both target and query are miRNA sequences of the same length.
    DOES NOT WORK FOR AN ACTUAL TARGET SEQUENCE, WHICH IS THE REVERSE COMPLEMENT."""
    
    mismatch_positions = [0 for i in range(len(intended_target))]
    wobble_positions = [0 for i in range(len(intended_target))]
    
    for i in range(len(intended_target)):
        if intended_target[i] != query[i]:
            mismatch_positions[i] = 1
            
        # it's a wobble if the intended target is an A and the query is an G
        # since the G in the query can bind to the U in the actual target sequence
        if intended_target[i] == 'A' and query[i] == 'G':
            wobble_positions[i] = 1
            mismatch_positions[i] = 0
        
        # it's a wobble if the intended target is a C and the query is a T
        # since the U in the query can bind to the G in the actual target sequence
        elif intended_target[i] == 'C' and query[i] == 'T':
            wobble_positions[i] = 1
            mismatch_positions[i] = 0

    # the RISC likes an "A" in the actual target sequence at position 1
    if intended_target[0] == 'T':
        mismatch_positions[0] = 0
        wobble_positions[0] = 0

    return mismatch_positions, wobble_positions

def sum_mismatches_in_regions(mismatch, wobble, regions=regions):
    """Accepts a list of individual mismatches and wobbles between the target sequence of a microRNA intended_target
    and a second attacking microRNA query.
    
    Returns a list of mismatches and wobbles summed over each region in regions."""
    
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

def count_mismatches_in_region(intended_target, query, regions=regions):
    """Returns a list of the counts of mismatches and wobbles in each region between
    the target sequence of a microRNA intended_target and a second attacking microRNA query.
    Assumes that both target and query are miRNA sequences of the same length,
    NOT the reverse complement.
    
    Returns a list of the counts of mismatches and wobbles in each region."""
    
    mismatch, wobble = get_mismatches_two_mirnas(intended_target, query)
    mismatch_counts, wobble_counts = sum_mismatches_in_regions(mismatch, wobble, regions)

    return mismatch_counts, wobble_counts

def get_mismatches_with_reverse_complement(sequence1, sequence2):
    """Returns a list of mismatches and wobbles between two DNA sequences.
    Assumes that sequence is a miRNA and target an actual target sequence of the same length,
    i.e. the reverse complement of a miRNA sequence."""
    
    # initialize the mismatch and wobble positions
    # by assuming all mismatches and no wobbles
    mismatch_positions = [1 for i in range(len(sequence2))]
    wobble_positions = [0 for i in range(len(sequence2))]

    complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
    wobbles = {"A": "none", "C": "none", "G": "T", "T": "G"}
    
    # make the target the complement
    sequence2 = sequence2[::-1]

    for i in range(len(sequence2)):
        if complements[sequence2[i]] == sequence1[i]:
            mismatch_positions[i] = 0
        if wobbles[sequence2[i]] == sequence1[i]:
            wobble_positions[i] = 1
            mismatch_positions[i] = 0

    return mismatch_positions, wobble_positions

def merge_identical_mirnas(df_expression, mirbase):
    """Merges identical miRNAs in the expression dataframe.
    Assumes log10 expression levels."""
    
    # First, we merge identical miRNAs
    mirbase["sequence_norm_short"] = mirbase["sequence_norm"].str[:18]

    # find out which miRNAs are identical
    mirna_groups = mirbase.groupby("sequence_norm_short").groups
        
    # add their expression levels
    for key, group in mirna_groups.items():
        # filter the group to those that are in the expression dataframe
        group = [x for x in group if x in df_expression.index]
        if len(group) > 1:
            # get the expression levels
            expression_levels = df_expression.loc[group, :]
            # make them linear
            expression_levels = 10**expression_levels
            # sum them
            expression_levels = np.log10(expression_levels.sum(axis=0))
            # set all mirnas in the group to the summed expression levels
            for mirna in group:
                for column in df_expression.columns:
                    # is the original value NaN? skip!
                    if np.isnan(df_expression.loc[mirna, column]):
                        continue
                    df_expression.loc[mirna, column] = expression_levels[column]
                
    return df_expression, mirna_groups