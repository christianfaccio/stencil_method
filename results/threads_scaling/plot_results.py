#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import glob
import re

# Pattern to match output files
files = glob.glob("output/threads_scaling/*.out")

data = []

for f in files:
    with open(f, 'r') as file:
        content = file.read()

        # Extract number of threads
        threads_match = re.search(r'CPUs per task \(OMP_NUM_THREADS\): (\d+)', content)
        # Extract total time
        time_match = re.search(r'Elapsed time: ([\d.]+) seconds', content)

        if threads_match and time_match:
            threads = int(threads_match.group(1))
            total_time = float(time_match.group(1))
            data.append((threads, total_time))

# Convert to DataFrame
results = pd.DataFrame(data, columns=["Threads", "Total_time"])
results.sort_values("Threads", inplace=True)

# Calculate speedup (relative to minimum threads)
baseline_time = results[results["Threads"] == results["Threads"].min()]["Total_time"].iloc[0]
results["Speedup"] = baseline_time / results["Total_time"]

# Create two plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Total time vs threads
ax1.plot(results["Threads"], results["Total_time"], marker="o")
ax1.set_xlabel("Number of Threads")
ax1.set_ylabel("Total Time (s)")
ax1.set_title("Total Time vs Number of Threads")
ax1.grid(True)

# Plot 2: Speedup vs threads
ax2.plot(results["Threads"], results["Speedup"], marker="o", label="Actual Speedup")
# Ideal speedup line
ax2.plot(results["Threads"], results["Threads"] / results["Threads"].min(), "--", alpha=0.5, label="Ideal Speedup")
ax2.set_xlabel("Number of Threads")
ax2.set_ylabel("Speedup")
ax2.set_title("Speedup vs Number of Threads")
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig("plots/threads_scaling.png", dpi=300, bbox_inches='tight')
print(f"Plot saved to plots/threads_scaling.png")
print(f"\nProcessed {len(results)} output files")
print(results)
plt.show()
