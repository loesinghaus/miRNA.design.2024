import concurrent.futures
import os
import subprocess
import datetime

lib_names = [
    'HEK293T_3UTR_r1',
]

input_dir = '4_filtered_alignments'
output_dir = '5_deduplicated'

# create output directory if it does not exist
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# Function to deduplicate a library
def deduplicate_lib(lib_name):
    print(f"Starting library {lib_name} at: {datetime.datetime.now()}")
    command = f'umi_tools dedup --method unique -I {input_dir}/{lib_name}_filter_sorted.bam ' + \
              f'--output-stats={output_dir}/{lib_name}_deduplicated -S {output_dir}/{lib_name}_deduplicated.bam > {output_dir}/{lib_name}_umi_tools.log'
    proc = subprocess.run(command, shell=True, stderr=subprocess.PIPE)
    if proc.stderr:
        print("Error:", proc.stderr.decode())
    print(f"Finished library {lib_name} at: {datetime.datetime.now()}")

# Use ProcessPoolExecutor to run deduplications in parallel
with concurrent.futures.ProcessPoolExecutor() as executor:
    # Map the deduplicate_lib function to each lib_name
    results = executor.map(deduplicate_lib, lib_names)

    # Process results if needed (here we just pass)
    for result in results:
        pass  