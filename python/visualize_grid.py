#!/usr/bin/env python3
"""
Script to visualize energy grid evolution from binary files
Usage: python visualize_grid.py <grid_size_x> <grid_size_y> [--animate]
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import argparse
import os

# GitHub background colors
GITHUB_DARK = '#0d1117'
GITHUB_LIGHT = '#ffffff'

def read_grid(filename, sizex, sizey):
    """Read binary grid file in float32 format"""
    with open(filename, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32, count=sizex*sizey)
        return data.reshape(sizey, sizex)

def plot_single_grid(filename, sizex, sizey, output_file=None, github_bg=True):
    """Plot a single grid snapshot"""
    grid = read_grid(filename, sizex, sizey)

    fig = plt.figure(figsize=(10, 8))
    ax = plt.gca()

    if github_bg:
        fig.patch.set_facecolor(GITHUB_DARK)
        ax.set_facecolor(GITHUB_DARK)
        # Set text colors for dark background
        ax.tick_params(colors='white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')
        ax.spines['bottom'].set_color('white')
        ax.spines['top'].set_color('white')
        ax.spines['left'].set_color('white')
        ax.spines['right'].set_color('white')

    im = plt.imshow(grid, cmap='hot', interpolation='nearest', origin='lower')
    cbar = plt.colorbar(im, label='Energy')

    if github_bg:
        cbar.ax.yaxis.label.set_color('white')
        cbar.ax.tick_params(colors='white')
        cbar.outline.set_edgecolor('white')

    plt.title(f'Energy Grid - {os.path.basename(filename)}')
    plt.xlabel('X')
    plt.ylabel('Y')

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight',
                   facecolor=fig.get_facecolor())
        print(f"Saved plot to {output_file}")
    else:
        plt.show()
    plt.close()

def create_animation(sizex, sizey, pattern='plane_*.bin', output_file='energy_evolution.mp4', github_bg=True):
    """Create animation from all grid files"""
    files = sorted(glob.glob(pattern))

    if not files:
        print(f"No files found matching pattern: {pattern}")
        return

    print(f"Found {len(files)} grid files")

    # Read all grids
    grids = [read_grid(f, sizex, sizey) for f in files]

    # Find global min/max for consistent colorbar
    vmin = min(g.min() for g in grids)
    vmax = max(g.max() for g in grids)

    fig, ax = plt.subplots(figsize=(10, 8))

    if github_bg:
        fig.patch.set_facecolor(GITHUB_DARK)
        ax.set_facecolor(GITHUB_DARK)
        ax.tick_params(colors='white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')
        ax.spines['bottom'].set_color('white')
        ax.spines['top'].set_color('white')
        ax.spines['left'].set_color('white')
        ax.spines['right'].set_color('white')

    # Initial plot
    im = ax.imshow(grids[0], cmap='hot', interpolation='nearest',
                   origin='lower', vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(im, ax=ax, label='Energy')

    if github_bg:
        cbar.ax.yaxis.label.set_color('white')
        cbar.ax.tick_params(colors='white')
        cbar.outline.set_edgecolor('white')

    title = ax.set_title(f'Energy Grid - Step 0')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    def update(frame):
        im.set_array(grids[frame])
        title.set_text(f'Energy Grid - Step {frame}')
        return [im, title]

    anim = animation.FuncAnimation(fig, update, frames=len(grids),
                                   interval=100, blit=True, repeat=True)

    # Save animation
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=10, metadata=dict(artist='Stencil Method'), bitrate=1800)
    anim.save(output_file, writer=writer, savefig_kwargs={'facecolor': fig.get_facecolor()})
    print(f"Animation saved to {output_file}")
    plt.close()

def create_subplot_grid(sizex, sizey, pattern='plane_*.bin', output_file='energy_snapshots.png', github_bg=True):
    """Create a grid of snapshots at different time steps"""
    files = sorted(glob.glob(pattern))

    if not files:
        print(f"No files found matching pattern: {pattern}")
        return

    print(f"Found {len(files)} grid files")

    # Select evenly spaced snapshots (max 9)
    n_snapshots = min(9, len(files))
    indices = np.linspace(0, len(files)-1, n_snapshots, dtype=int)
    selected_files = [files[i] for i in indices]

    # Read selected grids
    grids = [read_grid(f, sizex, sizey) for f in selected_files]

    # Find global min/max for consistent colorbar
    vmin = min(g.min() for g in grids)
    vmax = max(g.max() for g in grids)

    # Create subplot grid
    n_rows = int(np.ceil(np.sqrt(n_snapshots)))
    n_cols = int(np.ceil(n_snapshots / n_rows))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows))

    if github_bg:
        fig.patch.set_facecolor(GITHUB_DARK)

    if n_snapshots == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for i, (grid, idx) in enumerate(zip(grids, indices)):
        ax = axes[i]

        if github_bg:
            ax.set_facecolor(GITHUB_DARK)
            ax.tick_params(colors='white')
            ax.xaxis.label.set_color('white')
            ax.yaxis.label.set_color('white')
            ax.title.set_color('white')
            ax.spines['bottom'].set_color('white')
            ax.spines['top'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.spines['right'].set_color('white')

        im = ax.imshow(grid, cmap='hot', interpolation='nearest',
                      origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(f'Step {idx}')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')

    # Hide unused subplots
    for i in range(n_snapshots, len(axes)):
        axes[i].axis('off')

    # Add colorbar
    cbar = fig.colorbar(im, ax=axes, label='Energy', fraction=0.046, pad=0.04)

    if github_bg:
        cbar.ax.yaxis.label.set_color('white')
        cbar.ax.tick_params(colors='white')
        cbar.outline.set_edgecolor('white')

    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight', facecolor=fig.get_facecolor())
    print(f"Snapshot grid saved to {output_file}")
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize energy grid evolution')
    parser.add_argument('sizex', type=int, help='Grid size in X direction')
    parser.add_argument('sizey', type=int, help='Grid size in Y direction')
    parser.add_argument('--pattern', default='plane_*.bin', help='File pattern to match (default: plane_*.bin)')
    parser.add_argument('--animate', action='store_true', help='Create animation (requires ffmpeg)')
    parser.add_argument('--snapshots', action='store_true', help='Create snapshot grid')
    parser.add_argument('--single', type=str, help='Plot a single file')
    parser.add_argument('--output', type=str, help='Output filename')
    parser.add_argument('--no-github-bg', action='store_true', help='Use white background instead of GitHub dark')

    args = parser.parse_args()

    github_bg = not args.no_github_bg

    if args.single:
        output = args.output or args.single.replace('.bin', '.png')
        plot_single_grid(args.single, args.sizex, args.sizey, output, github_bg)
    elif args.animate:
        output = args.output or 'energy_evolution.mp4'
        create_animation(args.sizex, args.sizey, args.pattern, output, github_bg)
    elif args.snapshots:
        output = args.output or 'energy_snapshots.png'
        create_subplot_grid(args.sizex, args.sizey, args.pattern, output, github_bg)
    else:
        # Default: create snapshot grid
        output = args.output or 'energy_snapshots.png'
        create_subplot_grid(args.sizex, args.sizey, args.pattern, output, github_bg)
