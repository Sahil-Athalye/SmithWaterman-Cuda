#!/usr/bin/env python3
"""
run_pairwise_alignments.py

Workflow:
1. Parse command-line arguments:
   - input_root: root directory containing subfolders of single-sequence FASTA files (e.g., pairwise_fasta)
   - cuda_binary: path to your CUDA alignment executable (smithWaterman)
   - output_root: directory where MSF outputs will be written
   - strip_prefix: whether to strip the folder name prefix from sequence names in MSF files

2. For each subdirectory under input_root:
   a. Create a corresponding subdirectory under output_root with the same name.
   b. List all FASTA files (*.fa, *.fasta) in the subdirectory.
   c. For each unordered pair of FASTA files:
      i. Construct an output filename: <subdir>_<seq1>__<seq2>.msf
     ii. Invoke the CUDA binary: [cuda_binary, seq1_path, seq2_path], capture stdout.
    iii. Process the output to strip prefixes from sequence names if requested
     iv. Write the processed output to the output file in the output subdirectory.

3. Print a summary of total alignments processed.

Usage:
    python run_pairwise_alignments.py \
        /path/to/pairwise_fasta \
        /path/to/smithWaterman \
        /path/to/output_msf \
        [--strip-prefix]
"""
import os
import argparse
import itertools
import subprocess
import sys
import re

def strip_prefixes_from_msf(msf_content, prefix):
    """
    Remove the folder name prefix from sequence names in MSF content
    E.g., change "BB11001_1aab_" to "1aab_" throughout the MSF file
    """
    # Pattern to match: Name lines and sequence data lines
    pattern = rf'(^|\s)({prefix}_)([^\s]+)'
    
    # Process line by line
    lines = msf_content.split('\n')
    processed_lines = []
    
    for line in lines:
        # Replace the prefix in the line
        processed_line = re.sub(pattern, r'\1\3', line)
        processed_lines.append(processed_line)
    
    return '\n'.join(processed_lines)

def main():
    parser = argparse.ArgumentParser(description="Run CUDA Smith-Waterman on all FASTA pairs and collect MSF outputs.")
    parser.add_argument('input_root', help='Root directory of FASTA subfolders')
    parser.add_argument('cuda_binary', help='Path to CUDA alignment executable (smithWaterman)')
    parser.add_argument('output_root', help='Directory to store MSF output files')
    parser.add_argument('--strip-prefix', action='store_true', 
                        help='Strip folder name prefix from sequence names in MSF files')
    args = parser.parse_args()

    input_root = os.path.abspath(args.input_root)
    cuda_bin = os.path.abspath(args.cuda_binary)
    output_root = os.path.abspath(args.output_root)
    os.makedirs(output_root, exist_ok=True)

    total = 0
    # Iterate over subdirectories
    for sub in sorted(os.listdir(input_root)):
        sub_in = os.path.join(input_root, sub)
        if not os.path.isdir(sub_in):
            continue
        # Create corresponding output subdir
        sub_out = os.path.join(output_root, sub)
        os.makedirs(sub_out, exist_ok=True)
        # Gather FASTA files
        fastas = [f for f in sorted(os.listdir(sub_in)) if f.lower().endswith(('.fa', '.fasta'))]
        # For each unordered pair
        for f1, f2 in itertools.combinations(fastas, 2):
            path1 = os.path.join(sub_in, f1)
            path2 = os.path.join(sub_in, f2)
            # Derive sequence IDs (filename without extension)
            id1 = os.path.splitext(f1)[0]
            id2 = os.path.splitext(f2)[0]
            out_name = f"{sub}_{id1}__{id2}.msf"
            out_path = os.path.join(sub_out, out_name)
            # Run CUDA binary
            try:
                result = subprocess.run([cuda_bin, path1, path2], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
            except subprocess.CalledProcessError as e:
                sys.stderr.write(f"Error running {cuda_bin} on {path1}, {path2}: {e}\n")
                continue
            
            # Process the output if needed
            output_content = result.stdout
            if args.strip_prefix:
                output_content = strip_prefixes_from_msf(output_content, sub)
            
            # Write MSF output
            with open(out_path, 'w') as wf:
                wf.write(output_content)
            total += 1
            print(f"Wrote {out_path}")

    print(f"Processed {total} alignments across all folders.")

if __name__ == '__main__':
    main()