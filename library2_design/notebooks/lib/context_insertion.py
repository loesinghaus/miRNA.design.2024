import re
import pandas as pd

def distance_to_start_codon(pos, var_region_len=164, three_p_seq_len=55):
    """This function determines the distance from the insertion position to the start codon."""
    return three_p_seq_len + var_region_len - pos

def determine_generic_ins_positions(
        var_region_len=164, three_p_seq_len=55, no_of_inserts=1, target_len=21, dist_between=6):
    """This function determines the positions of the generic inserts sites in the oligo.
    The optimal distance for cooperativity is between 17 and 35 nt for seed sites. In the Gam paper,
    the distance between full-length target sites is literally 0. I will use
    a fixed distance of 6 nt between target sites that are assumed to be 21 nt long.
    This makes sure that the distance between insert sites is a multiple of 3.

    var_region_len: length of variable region in oligo
    three_p_seq_len: length of 3' sequence until the start codon after the variable region
    """

    # determine the total length of the insert into the oligo
    insert_len = no_of_inserts*target_len + (no_of_inserts-1)*dist_between

    # set the insertion position such that there are equal number of nucleotides on both sides of the insert
    insertion_pos = int((var_region_len-insert_len)/2)

    # check if the distance from the first insertion position to the start codon is a multiple of 3
    dist_to_start = distance_to_start_codon(insertion_pos, var_region_len, three_p_seq_len)  
    # if not, move the insertion position to the next multiple of 3
    if dist_to_start%3 != 0:
        insertion_pos += dist_to_start%3
    # make sure I did that right
    dist_to_start = distance_to_start_codon(insertion_pos, var_region_len, three_p_seq_len) 
    if dist_to_start%3 != 0:
        raise ValueError("The distance to the start codon is not a multiple of 3.")
    
    # determine the positions of the other insert sites
    insert_positions = [insertion_pos]
    for i in range(1, no_of_inserts):
        insert_positions.append(insert_positions[i-1] + target_len + dist_between - (target_len+dist_between)%3)

    return insert_positions

def insert_miRNA_sites(context, insert_positions, target_df, miRNA_names):
    """This function inserts the miRNA sites into the oligo.
    context: context sequence
    insert_positions: positions of the insert sites
    mirbase: mirbase dataframe
    miRNA_names: miRNA sequence
    """

    # get the miRNA target sites
    miRNA_targets = target_df.loc[miRNA_names, 'target'].values
    # get the length of the target sequences
    target_len = len(miRNA_targets[0])
    
    # get the position of the first ATG in each target site mod 3
    ATG_pos_mod3 = target_df.loc[miRNA_names, 'ATG_pos_mod3'].values
    ATG_pos_mod3 = [sublist[0] for sublist in ATG_pos_mod3 if sublist]
    
    # get the position of all ATGs in each target site
    ATG_pos_total = target_df.loc[miRNA_names, 'ATG_pos'].values
    # get the ATG count
    ATG_count = target_df.loc[miRNA_names, 'ATG_count'].values

    # shift the insert position by the mod3 of the ATG position
    insert_positions = [int(pos) - int(ATG_pos_mod3[i]) for i, pos in enumerate(insert_positions)]

    # insert the miRNA sites into the oligo
    for i, pos in enumerate(insert_positions):
        context = context[:pos] + miRNA_targets[i] + context[pos+target_len:]

    # find parts that are okay to edit afterwards
    # get all positions
    editable_positions = list(range(len(context)))
    # delete all positions that are part of the insertions
    remove_entries = [list(range(pos, pos+target_len)) for pos in insert_positions]
    # flatten the list
    remove_entries = [item for sublist in remove_entries for item in sublist]
    # remove the positions from the list of editable positions
    editable_positions = [pos for pos in editable_positions if pos not in remove_entries]

    # print(insert_positions)
    # print(editable_positions)

    # make sure that no additional ATGs were inserted into the context
    # find all ATGs in the oligo
    ATG_pos_oligo = [match.start() for match in re.finditer('ATG', context)]
    # find the position of all ATGs in the target site in the reference frame of the context
    ATG_pos_targets = [[ATG+pos for ATG in ATG_pos_total[i] if ATG_count[i] > 0] for i, pos in enumerate(insert_positions)]
    ATG_pos_targets = [item for sublist in ATG_pos_targets for item in sublist]

    # # if the lists are not identical, print both lists
    # if ATG_pos_oligo != ATG_pos_targets:
    #     print("Before fixing context:")
    #     print("-----------")
    #     print(ATG_pos_oligo)
    #     print(ATG_pos_targets)

    # find the position of all ATGs that were introduced by the insertions
    ATG_pos_inserted = [pos for pos in ATG_pos_oligo if pos not in ATG_pos_targets]
    # iterate over the positions of the inserted ATGs
    for pos in ATG_pos_inserted:
        # # if the A of the ATG is within two nucleotides of the beginning of an insertion site, replace it by a T
        # if any([(insert_pos-pos <= 2) and (pos < insert_pos) for insert_pos in insert_positions]):
        #     context = context[:pos] + 'T' + context[pos+1:]
        # # if the G of the ATG is within two nucleotides of the end of an insertion site, replace it by a C
        # if any([(pos+2-insert_pos-target_len <= 2) and (pos+2 >= insert_pos + target_len) for insert_pos in insert_positions]):
        #     context = context[:pos+2] + 'C' + context[pos+3:]
        # check if any of the ATG nucleotides are within editable positions
        # replace the first nucleotide that is not part of a miRNA target site by a C
        for i in range(pos, pos+3):
            if i in editable_positions:
                context = context[:i] + 'C' + context[i+1:]
                break
    
    ATG_pos_oligo = [match.start() for match in re.finditer('ATG', context)]
    # if ATG_pos_oligo != ATG_pos_targets:
    #     print("After fixing context:")
    #     print("-----------")
    #     print(ATG_pos_oligo)
    #     print(ATG_pos_targets)

    # assert that the lists are identical
    assert ATG_pos_oligo == ATG_pos_targets

    # check if 'AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA', 'ACTAAA' was inserted
    AATAAA = [match.start() for match in re.finditer(r'(?=(AATAAA|ATTAAA|AGTAAA|TATAAA|ACTAAA))', context)]
    # if so, replace a nucleotide that is not part of the miRNA target site by a C
    if AATAAA:
        for pos in AATAAA:
            # find the first position in the AATAAA sequence that is not part of a miRNA target site
            for i in range(pos, pos+6):
                if i in editable_positions:
                    if context[i] != 'C':
                        context = context[:i] + 'C' + context[i+1:]
                        break

    # check if restriction sites for BsaI were inserted
    BsaI = [match.start() for match in re.finditer(r'(?=(GAGACC|GGTCTC))', context)]
    change_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if BsaI:
        for pos in BsaI:
            # find the first position in the BsaI sequence that is not part of a miRNA target site
            for i in range(pos, pos+6):
                if i in editable_positions:
                    context = context[:i] + change_dict[context[i]] + context[i+1:]
                    break

    # verfiy that after all the modifications, the contet still contain all miRNA target sequences
    for target in miRNA_targets:
        assert target in context
            # # write targets and context to a text file
            # with open('targets.txt', 'w') as f:
            #     for target in miRNA_targets:
            #         f.write(target + '\n')
            #     f.write(context)

    return context