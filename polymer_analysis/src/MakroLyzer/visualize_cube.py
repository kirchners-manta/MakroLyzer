import argparse

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np


def read_cube(path):
    with open(path, "r") as fh:
        lines = fh.readlines()

    natoms = int(lines[2].split()[0])
    nx = int(lines[3].split()[0])
    ny = int(lines[4].split()[0])
    nz = int(lines[5].split()[0])
    data_start = 6 + abs(natoms)
    values = np.fromstring(" ".join(lines[data_start:]), sep=" ")
    expected = nx * ny * nz
    if values.size < expected:
        raise ValueError("Cube file does not contain enough volumetric values.")
    return values[:expected].reshape((nx, ny, nz), order="C")


def show_slices(grid):
    cx, cy, cz = (dim // 2 for dim in grid.shape)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].imshow(grid[:, :, cz].T, origin="lower", aspect="auto", cmap="viridis")
    axes[0].set_title(f"XY slice (z={cz})")
    axes[1].imshow(grid[:, cy, :].T, origin="lower", aspect="auto", cmap="viridis")
    axes[1].set_title(f"XZ slice (y={cy})")
    axes[2].imshow(grid[cx, :, :].T, origin="lower", aspect="auto", cmap="viridis")
    axes[2].set_title(f"YZ slice (x={cx})")
    for ax in axes:
        ax.set_xlabel("grid index")
        ax.set_ylabel("grid index")
    plt.tight_layout()
    plt.show()


def show_projections(grid):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].imshow(np.sum(grid, axis=2).T, origin="lower", aspect="auto", cmap="magma")
    axes[0].set_title("XY projection")
    axes[1].imshow(np.sum(grid, axis=1).T, origin="lower", aspect="auto", cmap="magma")
    axes[1].set_title("XZ projection")
    axes[2].imshow(np.sum(grid, axis=0).T, origin="lower", aspect="auto", cmap="magma")
    axes[2].set_title("YZ projection")
    for ax in axes:
        ax.set_xlabel("grid index")
        ax.set_ylabel("grid index")
    plt.tight_layout()
    plt.show()


def show_3d(grid, threshold=None, max_points=50000):
    nonzero = grid[grid > 0]
    if nonzero.size == 0:
        raise ValueError("Cube field contains only zeros; nothing to visualize in 3D.")

    if threshold is None:
        threshold = float(np.percentile(nonzero, 90))

    mask = grid >= threshold
    points = np.argwhere(mask)
    values = grid[mask]

    if points.size == 0:
        raise ValueError("No voxels are above the selected threshold.")

    if points.shape[0] > max_points:
        idx = np.linspace(0, points.shape[0] - 1, max_points, dtype=int)
        points = points[idx]
        values = values[idx]

    min_count = int(np.floor(values.min()))
    max_count = int(np.ceil(values.max()))
    boundaries = np.arange(min_count - 0.5, max_count + 1.5, 1.0)
    norm = mcolors.BoundaryNorm(boundaries, ncolors=plt.get_cmap("inferno").N, clip=True)

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(
        points[:, 0],
        points[:, 1],
        points[:, 2],
        c=values,
        cmap="inferno",
        s=100,
        alpha=0.55,
        edgecolors="none",
        norm=norm,
    )
    ax.set_title(f"3D HBond Density (threshold >= {threshold:.3g})")
    ax.set_xlabel("x index")
    ax.set_ylabel("y index")
    ax.set_zlabel("z index")
    cbar = fig.colorbar(sc, ax=ax, pad=0.1, label="H-bond count")
    cbar.set_ticks(np.arange(min_count, max_count + 1, 1))
    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Visualize MakroLyzer .cube scalar fields.")
    parser.add_argument("cube_file", help="Path to .cube file.")
    parser.add_argument(
        "--mode",
        choices=["slices", "projection", "3d"],
        default="slices",
        help="Visualization mode (default: slices).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Absolute threshold for mode=3d. If omitted, uses 90th percentile of non-zero values.",
    )
    parser.add_argument(
        "--max-points",
        type=int,
        default=50000,
        help="Maximum number of points shown in mode=3d (default: 50000).",
    )
    args = parser.parse_args()

    grid = read_cube(args.cube_file)
    if args.mode == "projection":
        show_projections(grid)
    elif args.mode == "3d":
        show_3d(grid, threshold=args.threshold, max_points=args.max_points)
    else:
        show_slices(grid)


if __name__ == "__main__":
    main()
