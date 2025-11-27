import concurrent.futures
import os
import subprocess

lib_names = [
    'HEK293T_3UTR_r1',
]

align_output_dir = '4_filtered_alignments'

def discard_sort_and_index(bam_file):
    print(f"Discard the second read of {bam_file}...")
    command = f'samtools view -b -f 64 {bam_file} > {bam_file.replace(".bam", "_discard.bam")}'
    proc = subprocess.run(command, shell=True, stderr=subprocess.PIPE)
    print(proc.stderr.decode())

    # sort
    print(f"Sorting {bam_file}...")
    command = f'samtools sort {bam_file.replace(".bam", "_discard.bam")} -o {bam_file.replace(".bam", "_sorted.bam")}'
    proc = subprocess.run(command, shell=True, stderr=subprocess.PIPE)
    print(proc.stderr.decode())
    
    # index
    print(f"Indexing {bam_file}...")
    command = f'samtools index {bam_file.replace(".bam", "_sorted.bam")}'
    proc = subprocess.run(command, shell=True, stderr=subprocess.PIPE)
    print(proc.stderr.decode())

    # remove the unsorted bam file
    print(f"Removing {bam_file}...")
    os.remove(bam_file)
    os.remove(bam_file.replace(".bam", "_discard.bam"))
    
with concurrent.futures.ProcessPoolExecutor() as executor:
    futures = [executor.submit(discard_sort_and_index,
             os.path.join(align_output_dir, lib_name + '_filter.bam')) for lib_name in lib_names]
    