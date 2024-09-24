import datetime
import os
import pandas
import multiprocessing
import shlex
import subprocess

import Bio
import Bio.SeqIO
import Bio.SeqRecord
import gzip

align_output_dir = '3_alignments'
ref_filepath = '2_references/references/'
data_dir = './2_fastq_split_rem_umi/'
temp_dir = './temp/'

lib_names = [
    'HEK293T_3UTR_r1',
]

if not os.path.exists(align_output_dir):
    os.mkdir(align_output_dir)
if not os.path.exists(temp_dir):
    os.mkdir(temp_dir)
    
# align
for lib_name in lib_names:
    print(f"Starting library {lib_name} at: {datetime.datetime.now()}")
    
    align_output_filename = lib_name + '.sam'
    align_output_filepath = os.path.join(align_output_dir, align_output_filename)
    align_output_logname = lib_name + '_align.log'
    align_output_logpath = os.path.join(align_output_dir, align_output_logname)

    file_path_R1 = [os.path.join(data_dir, lib_name + '_Read1.fq.gz')]
    file_path_R2 = [os.path.join(data_dir, lib_name + '_Read2.fq.gz')]
    
    # create temporary unpacked files
    temp_R1 = os.path.join(temp_dir, lib_name + '_Read1.fastq')
    temp_R2 = os.path.join(temp_dir, lib_name + '_Read2.fastq')
    
    print(f"Unpacking {file_path_R1[0]} and {file_path_R2[0]}")
    # unpack using the gzip library
    # chunk it to avoid memory issues
    chunk_size = 1024 * 1024 * 68
    with gzip.open(file_path_R1[0], 'rt') as f_in, open(temp_R1, 'w') as f_out:
        while True:
            chunk = f_in.read(chunk_size)
            if not chunk:
                break
            f_out.write(chunk)
    with gzip.open(file_path_R2[0], 'rt') as f_in, open(temp_R2, 'w') as f_out:
        while True:
            chunk = f_in.read(chunk_size)
            if not chunk:
                break
            f_out.write(chunk)
    
    print(f"Aligning {temp_R1} and {temp_R2} to {ref_filepath}")
    # -a: output all alignments for SE or unpaired PE
    # -M: mark shorter split hits as secondary
    # -t: number of threads
    # -o: output file
    command = f'bwa mem -a -M -t {multiprocessing.cpu_count()} -o {align_output_filepath} ' + \
            f'{ref_filepath} {temp_R1} {temp_R2}'
            
    command_split = shlex.split(command)
    log_str = ''
    
    with subprocess.Popen(command_split, stderr=subprocess.PIPE) as proc:
        myline = proc.stderr.readline().decode("utf-8")
        while myline:
            log_str += myline
            print(f"[{datetime.datetime.now().isoformat(sep=' ')}][{lib_name}] {myline.strip()}")
            myline = proc.stderr.readline().decode("utf-8")

    # Save log
    with open(align_output_logpath, 'w') as f:
        f.write(log_str)
        
    print("Converting to bam...")
    command = f'samtools view -b {align_output_filepath} > {align_output_filepath.replace(".sam", ".bam")}'
    proc = subprocess.run(command, shell=True, stderr=subprocess.PIPE)
        
    # remove temporary files
    os.remove(temp_R1)
    os.remove(temp_R2)
    
    # remove the sam file
    os.remove(align_output_filepath)
    
    print(f"Finished library {lib_name} at: {datetime.datetime.now()}")

print('Done')