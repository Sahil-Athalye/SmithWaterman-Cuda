import os
import sys
import subprocess
import argparse
import re
from itertools import combinations
import matplotlib.pyplot as plt

def find_subdirs(parent_dir):
    """Return full paths of immediate subdirectories of parent_dir."""
    return [
        os.path.join(parent_dir, d)
        for d in os.listdir(parent_dir)
        if os.path.isdir(os.path.join(parent_dir, d))
    ]

def find_files(dir_path):
    """Return full paths of all regular files in dir_path."""
    return [
        os.path.join(dir_path, f)
        for f in os.listdir(dir_path)
        if os.path.isfile(os.path.join(dir_path, f))
    ]

def main():
    parser = argparse.ArgumentParser(
        description="Run pairwise alignments and plot scores by folder."
    )
    parser.add_argument("parent_dir", help="Parent folder containing subfolders of sequences")
    parser.add_argument("executable", help="Path to smithWaterman executable")
    args = parser.parse_args()

    parent = os.path.abspath(args.parent_dir)
    exe = os.path.abspath(args.executable)

    if not os.path.isdir(parent):
        print(f"Error: '{parent}' is not a directory", file=sys.stderr)
        sys.exit(1)
    if not os.access(exe, os.X_OK):
        print(f"Error: '{exe}' is not executable", file=sys.stderr)
        sys.exit(1)

    all_scores = {}  # folder_name -> list of (label, score)

    for subdir in find_subdirs(parent):
        folder = os.path.basename(subdir)
        files = sorted(find_files(subdir))
        scores = []
        for f1, f2 in combinations(files, 2):
            cmd = [exe, f1, f2]
            print(f"Running: {' '.join(cmd)}")
            proc = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            # parse "Alignment score: N"
            m = re.search(r"Alignment score:\s*(\d+)", proc.stdout)
            score = int(m.group(1)) if m else None
            label = f"{os.path.basename(f1)} vs {os.path.basename(f2)}"
            scores.append((label, score))
        all_scores[folder] = scores

   # make an Images/ directory next to where the user ran the script
    images_dir = os.path.join(os.getcwd(), "Images")
    os.makedirs(images_dir, exist_ok=True)

    # for each folder, generate and save exactly one plot
    for folder, data in all_scores.items():
        if not data:
            continue
        labels, vals = zip(*data)
        x = range(len(vals))

        plt.figure()
        plt.plot(x, vals, marker='o')
        plt.title(f"Alignment scores in {folder}")
        plt.xlabel("Pair index")
        plt.ylabel("Alignment score")
        plt.xticks(x, labels, rotation=90)
        plt.tight_layout()

        # save it
        out_path = os.path.join(images_dir, f"{folder}.png")
        plt.savefig(out_path)
        plt.close()

if __name__ == "__main__":
    main()

