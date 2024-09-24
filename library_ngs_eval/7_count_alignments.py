import pysam
import sys
import pandas as pd
import os

def count_alignments(bam_file):
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Initialize a dictionary to store counts with all templates from the header
        counts = {ref: 0 for ref in bam.references}

        # Iterate over each alignment
        for read in bam:
            # Get reference name (template)
            # Ensure the read is aligned
            if read.reference_id != -1:  
                ref_name = bam.get_reference_name(read.reference_id)

                # Increment count for this reference name
                counts[ref_name] += 1

        return counts

def main():
    lib_names = [
        'HEK293T_3UTR_r1',
    ]
    output_path = '6_counts'

    # create folder if it doesn't exist
    if not os.path.exists('6_counts'):
        os.makedirs('6_counts')

    for lib_name in lib_names:
        print(f"Processing {lib_name}...")
        filename = f'5_deduplicated/{lib_name}_deduplicated.bam'
        counts = count_alignments(filename)

        # Convert the counts dictionary to a DataFrame
        df = pd.DataFrame(list(counts.items()), columns=['Reference', 'Count'])

        # Save the DataFrame to the specified path
        df.to_csv(os.path.join(output_path, lib_name+".csv"), index=False)

    # get all filenames in the output folder
    filenames = os.listdir(output_path)
    all_counts = []
    for filename in filenames:
        df = pd.read_csv(os.path.join(output_path, filename), index_col=0)
        df.rename(columns={'Count': filename.split('.')[0]}, inplace=True)
        all_counts.append(df)
    all_counts = pd.concat(all_counts, axis=1)
    all_counts.to_csv(os.path.join(output_path, 'all_counts.csv'))

if __name__ == "__main__":
    main()