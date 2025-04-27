#!/usr/bin/env python3
"""
split_msf_to_pairwise.py

Split multi-sequence MSF files into pairwise two-sequence MSF files.
For each MSF, creates a subdirectory under 'pairwise_msf' named after the file (without extension),
then writes one MSF per unordered pair, with headers only for those two sequences.

Usage:
    python split_msf_to_pairwise.py /path/to/msf_folder

Example:
  Input MSF: BB11001.msf (contains sequences: 1aab_, 1j46_A, 1k99_A, 2lef_A)

  Output folder structure:
    /path/to/msf_folder/pairwise_msf/BB11001/
        BB11001_1aab___1j46_A.msf
        BB11001_1aab___1k99_A.msf
        BB11001_1aab___2lef_A.msf
        BB11001_1j46_A__1k99_A.msf
        BB11001_1j46_A__2lef_A.msf
        BB11001_1k99_A__2lef_A.msf

This results in 6 files (choose(4,2)) in the BB11001 folder.
"""
import os
import argparse
import itertools


def read_msf(path):
    """
    Read an MSF file. Return header_lines (up to and including '//') and list of blocks.
    Each block is a list of sequence lines (strings without newline).
    """
    header_lines = []
    blocks = []
    with open(path) as f:
        # Read header
        for line in f:
            header_lines.append(line.rstrip('\n'))
            if line.strip() == '//':
                break
        # Read blocks
        current = []
        for line in f:
            if line.strip() == '':
                if current:
                    blocks.append(current)
                    current = []
            else:
                current.append(line.rstrip('\n'))
        if current:
            blocks.append(current)
    return header_lines, blocks


def get_ids_from_block(block):
    """
    Extract sequence IDs from the first token of each line in a block.
    """
    ids = []
    for line in block:
        parts = line.split()
        if parts:
            ids.append(parts[0])
    return ids


def filter_header(header_lines, seq1, seq2):
    """
    Filter header lines to include only Name: entries for seq1 and seq2.
    Keep all non-Name lines intact.
    """
    filtered = []
    for line in header_lines:
        stripped = line.lstrip()
        if stripped.startswith('Name:'):
            parts = stripped.split()
            if len(parts) >= 2 and parts[1] in (seq1, seq2):
                # preserve original indentation
                indent = len(line) - len(stripped)
                filtered.append(' ' * indent + stripped)
        else:
            filtered.append(line)
    return filtered


def write_pair_msf(out_path, header_lines, blocks, seq1, seq2):
    """
    Write a two-sequence MSF file containing only seq1 and seq2.
    """
    filtered_header = filter_header(header_lines, seq1, seq2)
    with open(out_path, 'w') as w:
        # Write filtered header
        for line in filtered_header:
            w.write(line + '\n')
        w.write('\n')
        # Write each block, only lines for seq1 and seq2
        for block in blocks:
            for line in block:
                name = line.split()[0]
                if name == seq1 or name == seq2:
                    w.write(line + '\n')
            w.write('\n')


def main():
    parser = argparse.ArgumentParser(
        description="Split MSF files into pairwise MSF for scoring."
    )
    parser.add_argument('input_dir', help="Directory containing .msf files")
    parser.add_argument('--outname', default='pairwise_msf',
                        help="Name of output subfolder")
    args = parser.parse_args()

    in_dir = os.path.abspath(args.input_dir)
    out_root = os.path.join(in_dir, args.outname)
    os.makedirs(out_root, exist_ok=True)

    for fname in os.listdir(in_dir):
        if not fname.lower().endswith('.msf'):
            continue
        base = os.path.splitext(fname)[0]
        msf_path = os.path.join(in_dir, fname)
        header, blocks = read_msf(msf_path)
        if not blocks:
            continue
        # infer sequence IDs from first block
        ids = get_ids_from_block(blocks[0])
        if len(ids) < 2:
            continue
        # create output folder for this MSF
        out_dir = os.path.join(out_root, base)
        os.makedirs(out_dir, exist_ok=True)
        # write each pair
        for seq1, seq2 in itertools.combinations(ids, 2):
            out_fname = f"{base}_{seq1}__{seq2}.msf"
            out_path = os.path.join(out_dir, out_fname)
            write_pair_msf(out_path, header, blocks, seq1, seq2)

    print("For each MSF, created a folder under {}/ with pairwise MSFs.".format(out_root))

if __name__ == '__main__':
    main()
