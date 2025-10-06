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
    num_nodes = re.search(r'Running on (\d+) nodes', content)
    mpi_tasks = re.search(r'Total MPI tasks: (\d+)', content)
    tasks_per_node = re.search(r'Tasks per node: (\d+)', content)
    num_threads = re.search(r'CPUs per task \(OMP_NUM_THREADS\): (\d+)', content)
    

    # Extract timing information
    elapsed_time = re.search(r'Elapsed time: ([\d.]+) seconds', content)
    comp_time = re.search(r'Computation time: ([\d.]+) seconds', content)
    comm_time = re.search(r'Communication time: ([\d.e+-]+) seconds', content)

    if not all([job_id, num_nodes, mpi_tasks, tasks_per_node, num_threads, elapsed_time, comp_time, comm_time]):
        return None
    
    nodes = int(num_nodes.group(1))
    mpi_t = int(mpi_tasks.group(1))
    t_per_node = int(tasks_per_node.group(1))
    threads = int(num_threads.group(1))
    total_time = float(elapsed_time.group(1))
    computation_time = float(comp_time.group(1))
    communication_time = float(comm_time.group(1))

    # Calculate speedup and efficiency (you'll need to specify the baseline)
    # For now, using threads as a placeholder - adjust as needed
    speedup = None  # Will be calculated after collecting all data
    efficiency = None # Will be calculated after collecting all data
    return {
        'job_id': job_id.group(1),
        'job_name': job_name.group(1) if job_name else '',
        'nodes': nodes,
        'mpi_tasks': mpi_t,
        'tasks/node': t_per_node,
        'threads': threads,
        'total_time': total_time,
        'computation_time': computation_time,
        'communication_time': communication_time,
        'speedup': speedup,
        'efficiency': efficiency
    }

def main():
    # Find all output files
    threads_scaling_output_files = glob.glob('output/threads_scaling/*.out')
    weak_scaling_output_files = glob.glob('output/weak_scaling/*.out')
    strong_scaling_output_files = glob.glob('output/strong_scaling/*.out')

    # Parse all files
    threads_results = []
    weak_results = []
    strong_results = []
    for filepath in sorted(threads_scaling_output_files):
        threads_data = parse_output_file(filepath)
        if threads_data:
            threads_results.append(threads_data)
        else:
            print(f"Warning: Could not parse {filepath}", file=sys.stderr)

    for filepath in sorted(weak_scaling_output_files):
        weak_data = parse_output_file(filepath)
        if weak_data:
            weak_results.append(weak_data)
        else:
            print(f"Warning: Could not parse {filepath}", file=sys.stderr)

    for filepath in sorted(strong_scaling_output_files):
        strong_data = parse_output_file(filepath)
        if strong_data:
            strong_results.append(strong_data)
        else:
            print(f"Warning: Could not parse {filepath}", file=sys.stderr)
    
    # Calculate speedup relative to single thread (or minimum threads)
    # Find baseline (typically single thread or minimum thread count)
    threads_baseline = min(threads_results, key=lambda x: x['threads'])
    threads_baseline_time = threads_baseline['total_time']
    weak_baseline = min(weak_results, key=lambda x: x['nodes'])
    weak_baseline_time = weak_baseline['total_time']
    strong_baseline = min(strong_results, key=lambda x: x['nodes'])
    strong_baseline_time = strong_baseline['total_time']

    for result in threads_results:
        result['speedup'] = threads_baseline_time / result['total_time']
    for result in weak_results:
        result['speedup'] = weak_baseline_time / result['total_time']
    for result in strong_results:
        result['speedup'] = strong_baseline_time / result['total_time']

    # Calculate efficiency -> efficiency = speedup/num_threads for threads, speedup for weak (should be ~1), speedup/num_nodes for strong
    for result in threads_results:
        result['efficiency'] = result['speedup'] / result['threads']
    for result in weak_results:
        result['efficiency'] = result['speedup']  # For weak scaling, efficiency is just the speedup (ideally ~1.0)
    for result in strong_results:
        result['efficiency'] = result['speedup'] / result['nodes']

    # Write to CSV
    threads_output_csv = 'results/threads_results.csv'
    weak_output_csv = 'results/weak_results.csv'
    strong_output_csv = 'results/strong_results.csv'
    fieldnames = ['job_id', 'job_name', 'nodes', 'mpi_tasks', 'tasks/node', 'threads', 'total_time',
                  'computation_time', 'communication_time', 'speedup', 'efficiency']

    with open(threads_output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(threads_results)

    with open(weak_output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(weak_results)

    with open(strong_output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(strong_results)

    print(f"Results written to results/ folder")

if __name__ == '__main__':
    main()
