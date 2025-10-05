#!/usr/bin/env python3
import re
import csv
import glob
import sys
from pathlib import Path

def parse_output_file(filepath):
    """Parse a SLURM output file and extract timing information."""
    with open(filepath, 'r') as f:
        content = f.read()

    # Extract job information
    job_id = re.search(r'Job ID: (\d+)', content)
    job_name = re.search(r'Job Name: (.+)', content)
    num_threads = re.search(r'CPUs per task \(OMP_NUM_THREADS\): (\d+)', content)

    # Extract timing information
    elapsed_time = re.search(r'Elapsed time: ([\d.]+) seconds', content)
    comp_time = re.search(r'Computation time: ([\d.]+) seconds', content)
    comm_time = re.search(r'Communication time: ([\d.e+-]+) seconds', content)

    if not all([job_id, num_threads, elapsed_time, comp_time, comm_time]):
        return None

    threads = int(num_threads.group(1))
    total_time = float(elapsed_time.group(1))
    computation_time = float(comp_time.group(1))
    communication_time = float(comm_time.group(1))

    # Calculate speedup (you'll need to specify the baseline)
    # For now, using threads as a placeholder - adjust as needed
    speedup = None  # Will be calculated after collecting all data

    return {
        'job_id': job_id.group(1),
        'job_name': job_name.group(1) if job_name else '',
        'threads': threads,
        'total_time': total_time,
        'computation_time': computation_time,
        'communication_time': communication_time,
        'speedup': speedup
    }

def main():
    # Find all output files
    output_pattern = sys.argv[1] if len(sys.argv) > 1 else 'output/threads_scaling/*.out'
    output_files = glob.glob(output_pattern)

    if not output_files:
        print(f"No output files found matching pattern: {output_pattern}")
        sys.exit(1)

    # Parse all files
    results = []
    for filepath in sorted(output_files):
        data = parse_output_file(filepath)
        if data:
            results.append(data)
        else:
            print(f"Warning: Could not parse {filepath}", file=sys.stderr)

    if not results:
        print("No valid results found")
        sys.exit(1)

    # Calculate speedup relative to single thread (or minimum threads)
    # Find baseline (typically single thread or minimum thread count)
    baseline = min(results, key=lambda x: x['threads'])
    baseline_time = baseline['total_time']

    for result in results:
        result['speedup'] = baseline_time / result['total_time']

    # Write to CSV
    output_csv = 'results/threads_scaling/threads_results.csv'
    fieldnames = ['job_id', 'job_name', 'threads', 'total_time',
                  'computation_time', 'communication_time', 'speedup']

    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"Results written to {output_csv}")
    print(f"Processed {len(results)} output files")

if __name__ == '__main__':
    main()
