#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import glob
import re

threads_results = pd.read_csv('results/threads_results.csv')
weak_results = pd.read_csv('results/weak_results.csv')
strong_results = pd.read_csv('results/strong_results.csv')

# -----------------------------------------------------------

# Threads Scaling

# -----------------------------------------------------------

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 5))

# Plot 1: Total time vs threads
ax1.plot(threads_results['threads'], threads_results['total_time'], marker='o', label='Total Time')
ax1.plot(threads_results['threads'], threads_results['computation_time'], marker='s', label='Computation Time')
ax1.plot(threads_results['threads'], threads_results['communication_time'], marker='^', label='Communication Time')
ax1.set_xlabel('Number of Threads')
ax1.set_ylabel('Time (s)')
ax1.legend()
ax1.grid(True)

# Plot 2: Speedup vs threads
ax2.plot(threads_results['threads'], threads_results['speedup'], marker='o', label='Actual Speedup')
ax2.plot(threads_results['threads'], threads_results['threads'], "--", alpha=0.5, label='Ideal Speedup')
ax2.set_xlabel('Number of Threads')
ax2.set_ylabel('Speedup')
ax2.legend()
ax2.grid(True)

# Plot 3: Efficiency vs threads
ax3.plot(threads_results['threads'], threads_results['efficiency'], marker='o', label='Actual Efficiency')
ax3.axhline(y=1.0, color='b', linestyle='--', label='Ideal Efficiency')
ax3.set_xlabel('Number of threads')
ax3.set_ylabel('Efficiency')
ax3.set_ylim(0, 1.1)
ax3.legend()
ax3.grid(True)

plt.tight_layout()
plt.savefig("plots/threads_scaling.png", dpi=300, bbox_inches='tight')
print(f"Plot saved to plots/threads_scaling.png")

# -----------------------------------------------------------

# Weak Scaling

# -----------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

# Plot 1: Total time vs nodes
ax1.plot(weak_results['nodes'], weak_results['total_time'], marker='o', label='Total Time')
ax1.plot(weak_results['nodes'], weak_results['computation_time'], marker='s', label='Computation Time')
ax1.plot(weak_results['nodes'], weak_results['communication_time'], marker='^', label='Communication Time')
ax1.set_xlabel('Number of Nodes')
ax1.set_ylabel('Time (s)')
ax1.legend()
ax1.grid(True)

# Plot 2: Efficiency vs nodes
ax2.plot(weak_results['nodes'], weak_results['efficiency'], marker='o', label='Actual Efficiency')
ax2.axhline(y=1.0, color='b', linestyle='--', label='Ideal Efficiency')
ax2.set_xlabel('Number of Nodes')
ax2.set_ylabel('Efficiency')
ax2.set_ylim(0, 1.1)
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig("plots/weak_scaling.png", dpi=300, bbox_inches='tight')
print(f"Plot saved to plots/weak_scaling.png")

# -----------------------------------------------------------

# Strong Scaling

# -----------------------------------------------------------

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 5))

# Plot 1: Total time vs nodes
ax1.plot(strong_results['nodes'], strong_results['total_time'], marker='o', label='Total Time')
ax1.plot(strong_results['nodes'], strong_results['computation_time'], marker='s', label='Computation Time')
ax1.plot(strong_results['nodes'], strong_results['communication_time'], marker='^', label='Communication Time')
ax1.set_xlabel('Number of Nodes')
ax1.set_ylabel('Time (s)')
ax1.legend()
ax1.grid(True)

# Plot 2: Speedup vs nodes
ax2.plot(strong_results['nodes'], strong_results['speedup'], marker='o', label='Actual Speedup')
ax2.plot(strong_results['nodes'], strong_results['nodes'], "--", alpha=0.5, label='Ideal Speedup')
ax2.set_xlabel('Number of Nodes')
ax2.set_ylabel('Speedup')
ax2.legend()
ax2.grid(True)

# Plot 3: Efficiency vs nodes
ax3.plot(strong_results['nodes'], strong_results['efficiency'], marker='o', label='Actual Efficiency')
ax3.axhline(y=1.0, color='b', linestyle='--', label='Ideal Efficiency')
ax3.set_xlabel('Number of Nodes')
ax3.set_ylabel('Efficiency')
ax3.set_ylim(0, 1.1)
ax3.legend()
ax3.grid(True)

plt.tight_layout()
plt.savefig("plots/strong_scaling.png", dpi=300, bbox_inches='tight')
print(f"Plot saved to plots/strong_scaling.png")
