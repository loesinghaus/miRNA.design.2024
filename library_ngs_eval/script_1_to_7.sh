#!/bin/bash
echo "Start time: $(date)" > time.txt

echo "Starting 1_split_UTR_rem_umis.py: $(date)" >> time.txt
python -u 1_split_UTR_rem_umis.py > log/1_split_rem_umis.txt
echo "Finished 1_split_UTR_rem_umis.py: $(date)" >> time.txt

echo "Starting 3_align_to_references.py: $(date)" >> time.txt
python -u 3_align_to_references.py > log/3_align_to_references.txt
echo "Finished 3_align_to_references.py: $(date)" >> time.txt

echo "Starting 4_filter_alignments.py: $(date)" >> time.txt
python -u 4_filter_alignments.py > log/4_filter_alignments.txt
echo "Finished 4_filter_alignments.py: $(date)" >> time.txt

echo "Starting 5_discard_sort_and_index.py: $(date)" >> time.txt
python -u 5_discard_sort_and_index.py > log/5_sort_and_index.txt
echo "Finished 5_discard_sort_and_index.py: $(date)" >> time.txt

echo "Starting 6_deduplicate_umis.py: $(date)" >> time.txt
python -u 6_deduplicate_umis.py > log/6_deduplicate_umis.txt
echo "Finished 6_deduplicate_umis.py: $(date)" >> time.txt

echo "Starting 7_count_alignments.py: $(date)" >> time.txt
python -u 7_count_alignments.py > log/7_count_alignments.txt
echo "Finished 7_count_alignments.py: $(date)" >> time.txt

echo "End time: $(date)" >> time.txt