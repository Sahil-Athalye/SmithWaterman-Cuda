#!/usr/bin/env python3
"""
split_tfa_to_single_fasta.py

This script takes a directory of multi-sequence BAliBASE .tfa files and
extracts each sequence into its own FASTA file, organizing outputs into
subdirectories named after the original .tfa file.

Example input folder:
    /home/warehouse/a.sahil/cse4059/SmithWaterman/Sequences/

Example output structure:
    /home/warehouse/a.sahil/cse4059/SmithWaterman/Sequences/
        pairwise_fasta/
            BB11001/
                BB11001_1aab_.fa
                BB11001_1j46_A.fa
                BB11001_1k99_A.fa
                BB11001_2lef_A.fa
            BB11002/
                ...

Usage:
    python split_tfa_to_single_fasta.py /path/to/Sequences

"""
import os
import argparse

def parse_fasta(path):
    """
    Parse a multi-FASTA (.tfa) file.
    Returns a list of (seq_id, sequence_string).
    """
    records = []
    header = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_lines)))
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records.append((header, ''.join(seq_lines)))
    return records


def write_fasta(path, entries, line_width=60):
    """
    Write entries—a list of (seq_id, seq_str)—to FASTA at path.
    """
    with open(path, 'w') as w:
        for seq_id, seq in entries:
            w.write(f'>{seq_id}\n')
            for i in range(0, len(seq), line_width):
                w.write(seq[i:i+line_width] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Extract each sequence from BAliBASE .tfa into its own FASTA file"
    )
    parser.add_argument('input_dir',
                        help="Path to folder containing BAliBASE .tfa files")
    parser.add_argument('--outname', default='pairwise_fasta',
                        help="Name of the output subfolder (default: pairwise_fasta)")
    args = parser.parse_args()

    in_dir = os.path.abspath(args.input_dir)
    out_root = os.path.join(in_dir, args.outname)
    os.makedirs(out_root, exist_ok=True)

    total = 0
    for fname in os.listdir(in_dir):
        if not fname.lower().endswith('.tfa'):
            continue
        base = os.path.splitext(fname)[0]
        tfa_path = os.path.join(in_dir, fname)
        records = parse_fasta(tfa_path)

        # create subdirectory for this .tfa file
        out_dir = os.path.join(out_root, base)
        os.makedirs(out_dir, exist_ok=True)

        # write one FASTA per sequence
        for seq_id, seq in records:
            out_fname = f"{base}_{seq_id}.fa"
            out_path = os.path.join(out_dir, out_fname)
            write_fasta(out_path, [(f"{base}_{seq_id}", seq)])
            total += 1

    print(f"Extracted {total} sequences into separate FASTA files under:\n  {out_root}")

if __name__ == '__main__':
    main()
