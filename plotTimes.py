import os
import sys
import subprocess
import argparse
import re
from itertools import combinations
import matplotlib.pyplot as plt
import time
from collections import defaultdict

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
        if os.path.isfile(os.path.join(dir_path, f)) and f.endswith('.fa')
    ]

def extract_timing(stderr_output):
    """Extract timing information from stderr output."""
    # Pattern for microsecond format: "CPU/GPU Execution time: 4743 μs (4743667 ns)"
    micro_match = re.search(r"Execution time: (\d+) μs", stderr_output)
    if micro_match:
        return float(micro_match.group(1)) / 1000.0  # Convert to milliseconds
    
    # Pattern for millisecond format: "CPU/GPU Execution time: 359.260 ms"
    ms_match = re.search(r"Execution time: (\d+\.\d+) ms", stderr_output)
    if ms_match:
        return float(ms_match.group(1))
    
    # If no match found, return None
    return None

def get_sequence_length(file_path):
    """Extract the length of a sequence from a FASTA file."""
    length = 0
    with open(file_path, 'r') as f:
        # Skip header line
        line = f.readline()
        if not line.startswith('>'):
            return 0
        
        # Count characters in sequence lines (skip whitespace)
        for line in f:
            if line.startswith('>'):
                break
            length += sum(1 for c in line if not c.isspace())
    
    return length

def run_alignment(cpu_exe, gpu_exe, file1, file2, verbose=True):
    """Run both CPU and GPU alignments, return timing and score."""
    result = {}
    
    # Get sequence lengths
    seq1_len = get_sequence_length(file1)
    seq2_len = get_sequence_length(file2)
    pair_size = seq1_len + seq2_len
    
    # Run CPU version
    if verbose:
        print(f"Running CPU: {os.path.basename(file1)} vs {os.path.basename(file2)}")
    
    cpu_proc = subprocess.run(
        [cpu_exe, file1, file2],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )
    
    # Extract CPU timing and score
    cpu_time = extract_timing(cpu_proc.stderr)
    score_match = re.search(r"Alignment score:\s*(\d+)", cpu_proc.stdout)
    score = int(score_match.group(1)) if score_match else 0
    
    # Run GPU version
    if verbose:
        print(f"Running GPU: {os.path.basename(file1)} vs {os.path.basename(file2)}")
    
    gpu_proc = subprocess.run(
        [gpu_exe, file1, file2],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )
    
    # Extract GPU timing
    gpu_time = extract_timing(gpu_proc.stderr)
    
    # Verify both implementations got the same score
    gpu_score_match = re.search(r"Alignment score:\s*(\d+)", gpu_proc.stdout)
    gpu_score = int(gpu_score_match.group(1)) if gpu_score_match else 0
    
    if score != gpu_score:
        print(f"Warning: Different scores for {file1} vs {file2} - CPU: {score}, GPU: {gpu_score}")
    
    return {
        'file1': os.path.basename(file1),
        'file2': os.path.basename(file2),
        'seq1_len': seq1_len,
        'seq2_len': seq2_len,
        'pair_size': pair_size,
        'score': score,
        'cpu_time': cpu_time,
        'gpu_time': gpu_time
    }

def sort_by_key(data_list, key):
    """Sort a list of dictionaries by a specific key."""
    return sorted(data_list, key=lambda x: x[key])

def plot_timing_comparison(results, images_dir, folder_name=None):
    """Plot CPU vs GPU timing comparison graphs."""
    # Sort by pair size
    sorted_results = sort_by_key(results, 'pair_size')
    
    # Extract data for plotting
    pair_sizes = [r['pair_size'] for r in sorted_results]
    cpu_times = [r['cpu_time'] for r in sorted_results]
    gpu_times = [r['gpu_time'] for r in sorted_results]
    
    # Create figure directory if needed
    os.makedirs(images_dir, exist_ok=True)
    
    # Plot 1: CPU vs GPU times by pair size
    plt.figure(figsize=(12, 8))
    plt.plot(pair_sizes, cpu_times, 'b-', marker='o', label='CPU')
    plt.plot(pair_sizes, gpu_times, 'r-', marker='x', label='GPU')
    
    title = "CPU vs GPU Execution Time by Sequence Pair Size"
    if folder_name:
        title += f" - {folder_name}"
    
    plt.title(title, fontsize=14)
    plt.xlabel("Sequence Pair Size (Total Characters)", fontsize=12)
    plt.ylabel("Execution Time (ms)", fontsize=12)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    file_prefix = folder_name if folder_name else "all"
    filename = f"{file_prefix}_timing_comparison.png"
    plt.savefig(os.path.join(images_dir, filename), dpi=150)
    
    # Plot 2: Log scale for better visualization of both small and large sequences
    plt.figure(figsize=(12, 8))
    plt.plot(pair_sizes, cpu_times, 'b-', marker='o', label='CPU')
    plt.plot(pair_sizes, gpu_times, 'r-', marker='x', label='GPU')
    
    title = "CPU vs GPU Execution Time by Sequence Pair Size (Log Scale)"
    if folder_name:
        title += f" - {folder_name}"
    
    plt.title(title, fontsize=14)
    plt.xlabel("Sequence Pair Size (Total Characters)", fontsize=12)
    plt.ylabel("Execution Time (ms)", fontsize=12)
    plt.yscale('log')
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save log-scale plot
    filename = f"{file_prefix}_timing_comparison_log.png"
    plt.savefig(os.path.join(images_dir, filename), dpi=150)
    
    # Plot 3: Speedup ratio (CPU time / GPU time)
    plt.figure(figsize=(12, 8))
    speedup = [c/g if g > 0 else 1 for c, g in zip(cpu_times, gpu_times)]
    plt.plot(pair_sizes, speedup, 'g-', marker='o')
    
    # Add horizontal line at y=1 (where CPU and GPU have equal performance)
    plt.axhline(y=1, color='k', linestyle='--', alpha=0.7, label='Equal Performance')
    
    title = "CPU/GPU Speedup Ratio by Sequence Pair Size"
    if folder_name:
        title += f" - {folder_name}"
    
    plt.title(title, fontsize=14)
    plt.xlabel("Sequence Pair Size (Total Characters)", fontsize=12)
    plt.ylabel("Speedup Ratio (CPU time / GPU time)", fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    plt.tight_layout()
    
    # Save speedup plot
    filename = f"{file_prefix}_speedup_ratio.png"
    plt.savefig(os.path.join(images_dir, filename), dpi=150)
    
    plt.close('all')
    
    # Plot 4: Scatterplot with annotations for outliers
    plt.figure(figsize=(12, 8))
    
    # Create scatterplot
    plt.scatter(cpu_times, gpu_times, alpha=0.7)
    
    # Add diagonal line for equal performance
    max_time = max(max(cpu_times), max(gpu_times))
    plt.plot([0, max_time], [0, max_time], 'k--', alpha=0.7, label='Equal Performance')
    
    # Add annotations for significant outliers
    threshold = 5  # Annotate points where one implementation is 5x faster than the other
    for i, result in enumerate(sorted_results):
        if gpu_times[i] == 0:
            continue
        ratio = cpu_times[i] / gpu_times[i]
        if ratio > threshold or ratio < 1/threshold:
            plt.annotate(
                f"{result['file1']} vs {result['file2']}",
                (cpu_times[i], gpu_times[i]),
                fontsize=8,
                alpha=0.8,
                textcoords="offset points",
                xytext=(5, 5),
                ha='left'
            )
    
    title = "CPU vs GPU Execution Time Comparison"
    if folder_name:
        title += f" - {folder_name}"
    
    plt.title(title, fontsize=14)
    plt.xlabel("CPU Execution Time (ms)", fontsize=12)
    plt.ylabel("GPU Execution Time (ms)", fontsize=12)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    plt.tight_layout()
    
    # Save comparison plot
    filename = f"{file_prefix}_time_scatter.png"
    plt.savefig(os.path.join(images_dir, filename), dpi=150)
    plt.close()

def write_csv(data, filepath):
    """Write benchmark data to a CSV file without pandas."""
    # Extract all possible keys from data
    keys = set()
    for item in data:
        keys.update(item.keys())
    keys = sorted(list(keys))
    
    with open(filepath, 'w') as f:
        # Write header
        f.write(','.join(keys) + '\n')
        
        # Write data
        for item in data:
            row = [str(item.get(key, '')) for key in keys]
            f.write(','.join(row) + '\n')
    
def main():
    parser = argparse.ArgumentParser(
        description="Benchmark and compare CPU vs GPU Smith-Waterman implementations."
    )
    parser.add_argument("parent_dir", help="Parent folder containing subfolders of FASTA sequences")
    parser.add_argument("cpu_exe", help="Path to CPU Smith-Waterman executable")
    parser.add_argument("gpu_exe", help="Path to GPU Smith-Waterman executable")
    parser.add_argument("--output", "-o", default="results", help="Output directory for plots and data")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print progress information")
    args = parser.parse_args()

    parent = os.path.abspath(args.parent_dir)
    cpu_exe = os.path.abspath(args.cpu_exe)
    gpu_exe = os.path.abspath(args.gpu_exe)
    output_dir = os.path.abspath(args.output)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    images_dir = os.path.join(output_dir, "images")
    os.makedirs(images_dir, exist_ok=True)
    
    # Validate inputs
    if not os.path.isdir(parent):
        print(f"Error: '{parent}' is not a directory", file=sys.stderr)
        sys.exit(1)
    if not os.access(cpu_exe, os.X_OK):
        print(f"Error: '{cpu_exe}' is not executable", file=sys.stderr)
        sys.exit(1)
    if not os.access(gpu_exe, os.X_OK):
        print(f"Error: '{gpu_exe}' is not executable", file=sys.stderr)
        sys.exit(1)

    # List to store all results
    all_results = []
    
    # Dictionary to store results by folder
    folder_results = defaultdict(list)

    # Process all subdirectories
    for subdir in find_subdirs(parent):
        folder = os.path.basename(subdir)
        if args.verbose:
            print(f"\nProcessing folder: {folder}")
        
        # Find all FASTA files in current subdirectory
        files = sorted(find_files(subdir))
        
        # Process all combinations of files
        total_pairs = len(list(combinations(files, 2)))
        for i, (f1, f2) in enumerate(combinations(files, 2)):
            if args.verbose:
                print(f"Pair {i+1}/{total_pairs}")
            
            # Run both implementations and collect results
            result = run_alignment(cpu_exe, gpu_exe, f1, f2, args.verbose)
            
            # Add folder name to result
            result['folder'] = folder
            
            # Store result
            all_results.append(result)
            folder_results[folder].append(result)
    
    # Save raw data to CSV
    csv_path = os.path.join(output_dir, "sw_benchmark_results.csv")
    write_csv(all_results, csv_path)
    print(f"\nSaved raw benchmark data to {csv_path}")
    
    # Generate overall comparison plot
    plot_timing_comparison(all_results, images_dir)
    
    # Generate per-folder comparison plots
    for folder, results in folder_results.items():
        plot_timing_comparison(results, images_dir, folder)
    
    print(f"\nBenchmark completed successfully!")
    print(f"- Found {len(all_results)} sequence pairs across {len(folder_results)} folders")
    print(f"- Generated plots in {images_dir}")
    print(f"- Raw data saved to {csv_path}")

if __name__ == "__main__":
    main()