import datetime
import os
import pandas
import numpy as np
import pysam
import concurrent.futures
import re
import traceback
import parasail
para_gap_open_penalty = 5
para_gap_extension_penalty = 1

align_input_dir = '3_alignments'
align_output_dir = '4_filtered_alignments'
create_alignment_df = False
create_alignment_overview_df = True

# This assumes the bwa aligner. If you use a different aligner, you may need to change this.
dist_threshold = 5
insertion_threshold = 2
deletion_threshold = 2
read_1_length = 150
read_2_length = 140

lib_names = [
    'HEK293T_3UTR_r1',
]

ref_filename = '2_references/references/references.fasta/'
ref_3UTR = pysam.FastaFile(ref_filename)

# Initialize dictionaries to hold the sequences from each FASTA file
ref_seqs_3UTR = {}

# Read sequences from the 3' UTR FASTA file
for ref in ref_3UTR.references:
    ref_seqs_3UTR[ref] = ref_3UTR.fetch(ref)

if not os.path.exists(align_output_dir):
    os.mkdir(align_output_dir)

def MD_to_edit_distance(MD_tag):
    # Extracting edit distance from MD tag
    # Example: MD:Z:7C3C6C7T105
    # Find all numbers in the MD:Z string, which represent lengths of matches
    match_lengths = re.findall(r'\d+', MD_tag)
    
    # Convert the found strings to integers and sum them up to get the total number of matches
    total_matches = sum(int(length) for length in match_lengths)
    
    return total_matches

def compute_edit_distance(query, quality, reference, inverse=False):
    if inverse:
        query = query[::-1]
        reference = reference[::-1]
        quality = quality[::-1]
    
    # Perform a global alignment
    result = parasail.nw_trace_striped_32(query, reference, para_gap_open_penalty, para_gap_extension_penalty, parasail.dnafull)
    aligned_query, aligned_reference = result.traceback.query, result.traceback.ref
    
    quality_mask = [q > 20 for q in quality]
    
    # is there a point where there are three consecutive parts with poor read quality (false in quality mask)?
    # if so, we will stop the alignment there
    stop_point = None
    for i in range(len(quality_mask) - 2):
        if not quality_mask[i] and not quality_mask[i + 1] and not quality_mask[i + 2]:
            stop_point = i
            break
    if stop_point is None:
        stop_point = len(quality_mask)
    
    # Initialize edit distance
    edit_distance = 0
    # Track the current position in the original query sequence
    query_pos = 0
    
    for q, r in zip(aligned_query, aligned_reference):
        if query_pos == stop_point or query_pos == len(query):
            break
        # Check for alignment gaps
        if q == "-":
            # Increase edit distance for indels if they occur in high-quality regions
            if quality_mask[query_pos]:
                edit_distance += 1
        else:
            # For mismatches (or gaps in the reference), consider the quality score
            if quality_mask[query_pos] and q != r and q != "N" and r != "N":
                edit_distance += 1
            # Increment query_pos only if there's no gap in the query
            query_pos += 1
    
    print(f"edit_distance: (inverse: {inverse}) {edit_distance}")
    
    return edit_distance

def get_insertions_and_deletions(CIGAR_tag):
    # Extracting insertions and deletions from CIGAR tag
    # Example: CIGAR: 150M1I140M
    # Find all numbers in the CIGAR string, which represent lengths of matches, insertions, and deletions
    match_lengths = re.findall(r'\d+', CIGAR_tag)
    
    # Find all letters in the CIGAR string, which represent operations
    operations = re.findall(r'[A-Z]', CIGAR_tag)
    
    # Get the indices of insertions and deletions
    insertion_indices = [i for i, operation in enumerate(operations) if operation == 'I']
    deletion_indices = [i for i, operation in enumerate(operations) if operation == 'D']
    
    # Get the lengths of insertions and deletions
    insertion_lengths = [int(match_lengths[i]) for i in insertion_indices]
    deletion_lengths = [int(match_lengths[i]) for i in deletion_indices]
    
    total_insertions = sum(insertion_lengths)
    total_deletions = sum(deletion_lengths)
    
    return total_insertions, total_deletions

def print_alignment_details(alignments, label):
    print(f"Details of {label}:")
    if alignments:
        for i, alignment in enumerate(alignments):
            print(f"Alignment {i+1}:")
            print(f"  Query Name: {alignment.query_name}")
            print(f"  Reference Name: {alignment.reference_name}")
            print(f"  Query Sequence: {alignment.query_sequence}")
            print(f"  Query Qualities: {alignment.query_qualities}")
            print(f"  Alignment Flag: {alignment.flag}")
            print(f"  Mapping Quality: {alignment.mapping_quality}")
            print(f"  Position: {alignment.reference_start}")
            print(f"  CIGAR String: {alignment.cigarstring}")
    else:
        print("No alignments found.")

def interpret_sam_flag(flag_value):
    flags = {
        0x1: "The read is part of a pair",
        0x2: "Both reads of the pair are properly aligned",
        0x4: "The read is unmapped",
        0x8: "The mate is unmapped",
        0x10: "The read is mapped to the reverse strand",
        0x20: "The mate is mapped to the reverse strand",
        0x40: "This is the first read in the pair",
        0x80: "This is the second read in the pair",
        0x100: "The read is not primary (secondary alignment)",
        0x200: "The read fails platform/vendor quality checks",
        0x400: "The read is PCR or optical duplicate"
    }

    print(f"FLAG: {flag_value} (binary: {bin(flag_value)})")
    for bitmask, description in flags.items():
        if flag_value & bitmask:
            print(f"- {description}")

# Note: Structure of a SAM alignment line:
# 0: QNAME ; Query template NAME: The name of the read or the read pair that was aligned. Read pairs have the same QNAME.
# 1: FLAG ; Bitwise FLAG: Flags indicating the alignment state of the read pair. See section 1.4 and section 2.1.
# 2: RNAME ; Reference sequence NAME: Name of the sequence in the reference genome.
# 3: POS ; 1-based leftmost mapping POSition: The leftmost coordinate of the clipped sequence.
# 4: MAPQ ; MAPping Quality: MAPQ = −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer. I can ignore this for this type of analysis.
# 5: CIGAR ; CIGAR string: Compact Idiosyncratic Gapped Alignment Report. Describes the alignment in detail. Look at this later to see if I need it.
# 6: RNEXT ; Ref. name of the mate/next read: Refers to the name of the reference sequence in the next read in the template.
# This field is set as ‘=’ when the other read in the template has the same RNAME value (see also RNEXT and PNEXT).
# This field is set as ‘*’ when the information is unavailable, or when the next read does not map to the same reference.
# 7: PNEXT ; Position of the mate/next read: The position of the primary alignment of the next read in the template. Should be 1 for my library.
# 8: TLEN ; observed Template LENgth: The observed template length.
# The sign of this field is determined by the relative strand orientations. (Should be 204 and -204 for my library.)
# 9: SEQ ; segment SEQuence: The primary read sequence. This field can be a ‘*’ when the sequence is not stored.
# 10: QUAL ; ASCII of Phred-scaled base QUALity+33: The base quality scores. This field can be a ‘*’ when the quality is not stored.
# 11+: OPT ; variable OPTional fields in the format TAG:VTYPE:VALUE: Optional fields in the format TAG:VTYPE:VALUE.
# 11: NM ; Edit distance to the reference: The edit distance (NM tag) equals the minimum number of edits needed to transform the read sequence into the reference sequence.
# 12: MD ; String for mismatching positions: This field describes the mismatching positions.
# 13: MC ; CIGAR string for mate/next read: CIGAR string for the mate/next read (MC tag). This field is set as ‘*’ when the information is unavailable.
# 14: AS ; Alignment Score: AS is the alignment score generated by the aligner. This can be very similar for clearly different alignments given I look at one nt mutations.
# 15: XS ; Suboptimal Alignment Score: XS is the suboptimal alignment score generated by the aligner.

def get_alignment_score(alignment):
    # Extracting AS tag from alignment
    return alignment.get_tag('AS')

def process_alignments(alignments, lib_name):
    # Process alignments
    aligned, discard_edit_distance, discarded_not_both, discarded_multiple_as = False, False, False, False
    valid_alignments = []
    edit_distance = 500
    
    # Get top and bottom strand alignments
    alignments_top = [alignment for alignment in alignments if not alignment.is_reverse]
    alignments_bottom = [alignment for alignment in alignments if alignment.is_reverse]
    
    # Get top and bottom strand alignment scores
    alignment_scores_top = [get_alignment_score(alignment) for alignment in alignments_top]
    alignment_scores_bottom = [get_alignment_score(alignment) for alignment in alignments_bottom]

    # Find max alignment scores and their indices
    max_score_top = max(alignment_scores_top, default=0)
    max_score_bottom = max(alignment_scores_bottom, default=0)
    
    # find the index of the alignments with the highest score (there may be multiple)
    max_indices_top = [i for i, score in enumerate(alignment_scores_top) if score == max_score_top]
    max_indices_bottom = [i for i, score in enumerate(alignment_scores_bottom) if score == max_score_bottom]
    
    # check if there are multiple alignments with the same score for both the top and bottom strand
    if (len(max_indices_top) > 1) and (len(max_indices_bottom) > 1):
        discarded_multiple_as = True
        max_index_top = -1
        max_index_bottom = -1
    elif (len(max_indices_top) > 1) and (len(max_indices_bottom) == 1):
        max_index_bottom = max_indices_bottom[0]
        # get the rname of the bottom strand alignment
        rname_bottom = alignments_bottom[max_index_bottom].reference_name
        rnames_top = [alignments_top[i].reference_name for i in max_indices_top]
        try:
            max_index_top = rnames_top.index(rname_bottom)
        except ValueError:
            max_index_top = -1
            discarded_not_both = True
    elif (len(max_indices_top) == 1) and (len(max_indices_bottom) > 1):
        max_index_top = max_indices_top[0]
        # get the rname of the top strand alignment
        rname_top = alignments_top[max_index_top].reference_name
        rnames_bottom = [alignments_bottom[i].reference_name for i in max_indices_bottom]
        try:
            max_index_bottom = rnames_bottom.index(rname_top)
        except ValueError:
            max_index_bottom = -1
            discarded_not_both = True
    elif (len(max_indices_top) < 1) or (len(max_indices_bottom) < 1):
        max_index_bottom = -1
        max_index_top = -1
        discarded_not_both = True
    else:
        max_index_top = max_indices_top[0]
        max_index_bottom = max_indices_bottom[0]
        
        # do these align to the same reference?
        if alignments_top[max_index_top].reference_name != alignments_bottom[max_index_bottom].reference_name:
            discarded_not_both = True
            max_index_top = -1
            max_index_bottom = -1
        
    if (max_index_top != -1) and (max_index_bottom != -1):        
        # get the read sequences and quality scores
        try:
            for alignment_top in alignments_top:
                read1_seq = str(alignment_top.query_sequence)
                read1_qual = alignment_top.query_qualities
                if read1_qual:
                    read1_qual = [q for q in read1_qual]
                    break
            for alignment_bottom in alignments_bottom:
                read2_seq = str(alignment_bottom.query_sequence)
                read2_qual = alignment_bottom.query_qualities
                if read2_qual:
                    read2_qual = [q for q in read2_qual]
                    break
        except TypeError:
            print_alignment_details(alignments_top, "alignments_top")
            print_alignment_details(alignments_bottom, "alignments_bottom")
        
        reference = ref_seqs_3UTR[alignments_top[max_index_top].reference_name]

        edit_distance_top = compute_edit_distance(read1_seq, read1_qual, reference[:150], inverse=False)
        edit_distance_bottom = compute_edit_distance(read2_seq, read2_qual, reference[-140:], inverse=True)
        
        # get the sum of the edit distance for the top and bottom strand alignments
        edit_distance = edit_distance_top + edit_distance_bottom
        
        # ---------------------------------------------------------
        # get the CIGAR tag for the top and bottom strand alignments
        CIGAR_top = alignments_top[max_index_top].cigarstring
        CIGAR_bottom = alignments_bottom[max_index_bottom].cigarstring
        
        insertions_top, deletions_top = get_insertions_and_deletions(CIGAR_top)
        insertions_bottom, deletions_bottom = get_insertions_and_deletions(CIGAR_bottom)
        
        # get the sum of the insertions and deletions for the top and bottom strand alignments
        insertions = insertions_top + insertions_bottom
        deletions = deletions_top + deletions_bottom
        
        # if it is below edit_dist_threshold, then it is a valid alignment
        if edit_distance < dist_threshold and insertions < insertion_threshold and deletions < deletion_threshold:
            aligned = True
            # append the top strand alignment
            valid_alignments.append(alignments_top[max_index_top])
            # append the bottom strand alignments
            valid_alignments.append(alignments_bottom[max_index_bottom])
        else:
            discard_edit_distance = True
            
    return aligned, valid_alignments, edit_distance, discard_edit_distance, discarded_not_both, discarded_multiple_as

def process_bam_file(input_filepath, output_filepath, lib_name):
    # Open the input BAM file for reading
    with pysam.AlignmentFile(input_filepath, "rb") as bam_in, \
            pysam.AlignmentFile(output_filepath, "wb", template=bam_in) as bam_out:
                
        # Statistics
        n_seqs = 0
        n_aligned = 0
        n_discarded_multiple_as = 0
        n_discarded_not_both = 0
        n_discarded_edit_distance = 0
        alignment_dict = {}

        current_qname = None
        alignments = []

        print(f"Starting library {lib_name} at: {datetime.datetime.now()}")
        # Iterate over each alignment
        for read in bam_in:
            # Check if this is a new query name
            if read.query_name != current_qname and alignments:
                n_seqs += 1
                
                if n_seqs % 1000000 == 0:
                    print(f"Processed {n_seqs} reads of {lib_name}...")
                
                # Process the collected alignments
                aligned, valid_alignments, edit_distance, discard_edit_distance, \
                    discarded_not_both, discarded_multiple_as = process_alignments(alignments, lib_name)
                    
                if aligned:
                    # Write valid alignments to output BAM file
                    for valid_alignment in valid_alignments:
                        bam_out.write(valid_alignment)
                    n_aligned += 1
                    
                if discarded_multiple_as:
                    n_discarded_multiple_as += 1
                if discard_edit_distance:
                    n_discarded_edit_distance += 1
                if discarded_not_both:
                    n_discarded_not_both += 1
                    
                alignment_dict[current_qname] = {
                    'reference_name': alignments[0].reference_name,
                    'aligned': aligned,
                    'edit_distance': edit_distance,
                    'discarded_not_both': discarded_not_both,
                    'discard_edit_distance': discard_edit_distance,
                    'discarded_multiple_as': discarded_multiple_as,
                }

                # Reset for the next set of alignments
                alignments = []

            # Add the read to the current alignments and update current_qname
            alignments.append(read)
            current_qname = read.query_name
        
        if create_alignment_overview_df:
            alignment_df = pandas.DataFrame.from_dict(alignment_dict, orient='index')
            if create_alignment_df:
                alignment_df.to_csv(os.path.join(align_output_dir, f"{lib_name}_alignment_stats.csv"))
            
            alignment_df_overview = pandas.DataFrame(index = alignment_df['reference_name'].unique(),
                                                 columns = ['aligned', 'discard_edit_distance', 'discarded_not_both', 'discarded_multiple_as'])
            alignment_df_overview['aligned'] = alignment_df.groupby('reference_name')['aligned'].sum()
            alignment_df_overview['discard_edit_distance'] = alignment_df.groupby('reference_name')['discard_edit_distance'].sum()
            alignment_df_overview['discarded_not_both'] = alignment_df.groupby('reference_name')['discarded_not_both'].sum()
            alignment_df_overview['discarded_multiple_as'] = alignment_df.groupby('reference_name')['discarded_multiple_as'].sum()
            alignment_df_overview['total'] = alignment_df_overview['aligned'] + alignment_df_overview['discard_edit_distance'] + alignment_df_overview['discarded_not_both'] + alignment_df_overview['discarded_multiple_as']
            alignment_df_overview['per_dis_edit'] = alignment_df_overview['discard_edit_distance'] / ( alignment_df_overview['total'] + 1 )
            alignment_df_overview['per_dis_not_both'] = alignment_df_overview['discarded_not_both'] / ( alignment_df_overview['total'] + 1 )
            alignment_df_overview['per_dis_multiple_as'] = alignment_df_overview['discarded_multiple_as'] / ( alignment_df_overview['total'] + 1 )
            alignment_df_overview.to_csv(os.path.join(align_output_dir, f"{lib_name}_alignment_stats_overview.csv"))
        
        print(f"Finished library {lib_name} at: {datetime.datetime.now()}")
        print(f"Total reads: {n_seqs}")
        print(f"Properly aligned reads: {n_aligned}")
        print(f"Discarded (not both aligned): {n_discarded_not_both}")
        print(f"Discarded (multiple alignment scores): {n_discarded_multiple_as}")
        print(f"Discarded (edit distance): {n_discarded_edit_distance}")
        
        # Output statistics
        with open(os.path.join(align_output_dir, f"{lib_name}_alignment_stats.txt"), "w") as stats_file:
            stats_file.write(f"Total reads: {n_seqs}\n")
            stats_file.write(f"Properly aligned reads: {n_aligned}\n")
            stats_file.write(f"Discarded (not both aligned): {n_discarded_not_both}\n")
            stats_file.write(f"Discarded (multiple alignments): {n_discarded_multiple_as}\n")
            stats_file.write(f"Discarded (edit distance): {n_discarded_edit_distance}\n")

def task_wrapper(func, *args, **kwargs):
    try:
        return func(*args, **kwargs)
    except Exception as e:
        # Capturing the traceback
        error_msg = f"An error occurred: {e}\n"
        error_msg += "".join(traceback.format_exception(None, e, e.__traceback__))
        return error_msg

# Function to prepare arguments for process_bam_file
def prepare_args(lib_name):
    align_input_filename = f'{lib_name}.bam'
    align_input_filepath = os.path.join(align_input_dir, align_input_filename)
    align_output_filename = f'{lib_name}_filter.bam'
    align_output_filepath = os.path.join(align_output_dir, align_output_filename)
    return (align_input_filepath, align_output_filepath, lib_name)

# Execute functions concurrently, modified to use the wrapper
with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
    # Prepare arguments for process_bam_file
    args = [prepare_args(lib_name) for lib_name in lib_names]
    
    # Schedule the wrapped process_bam_file function to run with the arguments
    futures = [executor.submit(task_wrapper, process_bam_file, *arg) for arg in args]
    
    # Use as_completed to handle the results as they are completed
    for future in concurrent.futures.as_completed(futures):
        result = future.result()
        if isinstance(result, str) and result.startswith("An error occurred:"):
            # If the result is a string starting with "An error occurred:", it's an error message with a traceback
            print(result)