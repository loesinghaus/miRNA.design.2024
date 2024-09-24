import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
import pandas as pd
import datetime

fixed_suffix = "CKDL240004006-1A_22GVGLLT3_L2"
lib_names = [
    'HEK293T_r1',
]

input_dir = '1_fastq'
output_dir = '2_fastq_split_rem_umi'

# Create output dir if necessary
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# this is for TruSeq Read2 in this design
internal_index_start = 26
internal_index_end = 32

# set the index sequences
threeUTRindex = "TCTATG"

# set up a dataframe to store the statistics
stats_df = pd.DataFrame(columns=["Library", "3UTR", "not found"])

def count_mismatching_characters(str1, str2):
    count = 0
    for char1, char2 in zip(str1, str2):
        if char1 != char2:
            count += 1
    return count

umi_len = 10
umi_pos_r2 = 0
# set this to None if you don't want to remove the UMI
# set it to -1 if the UMI is at the end of the read
umi_pos_r1 = None

# Iterate over read files
for lib_name in lib_names:
    no_reads_3UTR = 0
    no_reads_neither = 0
    total_reads = 0
    
    # Generate relevant names
    read1_filename = f'{lib_name}/{lib_name}_{fixed_suffix}_1.fq.gz'
    read2_filename = f'{lib_name}/{lib_name}_{fixed_suffix}_2.fq.gz'
    split_lib_name = lib_name.split('_')
    output1_filename_3UTR = f"{split_lib_name[0]}_3UTR_{split_lib_name[1]}_Read1.fq.gz"
    output2_filename_3UTR = f"{split_lib_name[0]}_3UTR_{split_lib_name[1]}_Read2.fq.gz"
    output1_filename_5UTR = f"{split_lib_name[0]}_5UTR_{split_lib_name[1]}_Read1.fq.gz"
    output2_filename_5UTR = f"{split_lib_name[0]}_5UTR_{split_lib_name[1]}_Read2.fq.gz"
    
    print(f"Starting library {lib_name} at: {datetime.datetime.now()}")
    
    with gzip.open(os.path.join(input_dir, read1_filename), 'rt') as input_read1_handle, \
        gzip.open(os.path.join(input_dir, read2_filename), 'rt') as input_read2_handle:

        with gzip.open(os.path.join(output_dir, output1_filename_3UTR), 'wt', compresslevel=5) as output1_3UTR_handle, \
            gzip.open(os.path.join(output_dir, output2_filename_3UTR), 'wt', compresslevel=5) as output2_3UTR_handle, \
            gzip.open(os.path.join(output_dir, output1_filename_5UTR), 'wt', compresslevel=5) as output1_5UTR_handle, \
            gzip.open(os.path.join(output_dir, output2_filename_5UTR), 'wt', compresslevel=5) as output2_5UTR_handle:
            for (title1, seq1, qual1), (title2, seq2, qual2) in zip(FastqGeneralIterator(input_read1_handle),\
                FastqGeneralIterator(input_read2_handle)):
                # Sanity check: titles must match
                assert(title1.split()[0] == title2.split()[0])
                # get the index on read2
                index = seq2[internal_index_start:internal_index_end]
                
                # ------------ umi processing ----------------------
                # get the UMI
                umi = seq2[umi_pos_r2:umi_pos_r2+umi_len]
                
                # add the UMI to the first portion of the title
                title1 = f"{title1.split()[0]}_{umi} {title1.split()[1]}"
                title2 = f"{title2.split()[0]}_{umi} {title2.split()[1]}"
                
                # remove the umi from the sequence
                seq2 = seq2[:umi_pos_r2] + seq2[umi_pos_r2+umi_len:]
                qual2 = qual2[:umi_pos_r2] + qual2[umi_pos_r2+umi_len:]
                
                # check if the read1 umi is not None
                if umi_pos_r1 is not None:
                    if umi_pos_r1 == -1:
                        seqlen = len(seq1)
                        seq1 = seq1[:seqlen-umi_len]
                        qual1 = qual1[:seqlen-umi_len]
                    else:
                        seq1 = seq1[:umi_pos_r1] + seq1[umi_pos_r1+umi_len:]
                        qual1 = qual1[:umi_pos_r1] + qual1[umi_pos_r1+umi_len:]
                
                # ------------ index processing ----------------------
                
                # check if it is a 3UTR index
                if count_mismatching_characters(threeUTRindex, index) < 2:
                    output1_3UTR_handle.write(f"@{title1}\n{seq1}\n+\n{qual1}\n")
                    output2_3UTR_handle.write(f"@{title2}\n{seq2}\n+\n{qual2}\n")
                    no_reads_3UTR += 1
                else:
                    no_reads_neither += 1
                total_reads += 1
                
                if total_reads % 1000000 == 0:
                    print(f"Processed {total_reads} reads")
                    
    print(f"Number of 3UTR reads: {no_reads_3UTR}")
    print(f"Number of reads w/o index: {no_reads_neither}")
    print(f"Total number of reads: {total_reads}")
    print(f"Finished library {lib_name} at: {datetime.datetime.now()}")
    
    new_row_df = pd.DataFrame({"Library": lib_name, "3UTR": no_reads_3UTR, "not found": no_reads_neither}, index=[0])
    stats_df = pd.concat([stats_df, new_row_df], ignore_index=True)
    
stats_df.to_csv(f"{output_dir}/stats.csv", index=False)