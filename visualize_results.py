import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_gc_content(metrics_path: str | Path, output_path: str | Path) -> None:
    df = pd.read_csv(metrics_path)
    required = {"organism", "gc_percentage"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in metrics file: {', '.join(sorted(missing))}")
    fig, ax = plt.subplots(figsize=(11, 6))
    bars = ax.bar(df["organism"], df["gc_percentage"], color="#4C72B0")
    ax.set(xlabel="Bacterial sequence", ylabel="GC content (%)", ylim=(0, 100),
           title="GC Content of 16S rRNA Sequences")
    ax.tick_params(axis="x", rotation=35)
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f"{bar.get_height():.1f}%", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_identity_heatmap(matrix_path: str | Path, output_path: str | Path) -> None:
    matrix = pd.read_csv(matrix_path, index_col=0)
    matrix = matrix.apply(pd.to_numeric)
    if matrix.empty or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Identity matrix must be a non-empty square matrix.")
    fig_size = max(7, min(16, 4 + len(matrix) * 0.65))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    image = ax.imshow(matrix.values, cmap="viridis", vmin=0, vmax=100)
    ax.set_xticks(range(len(matrix.columns)), matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(matrix.index)), matrix.index)
    ax.set_title("Pairwise 16S rRNA Sequence Identity (%)")
    fig.colorbar(image, ax=ax, label="Identity (%)", shrink=0.8)
    if len(matrix) <= 20:
        for row in range(len(matrix.index)):
            for col in range(len(matrix.columns)):
                ax.text(col, row, f"{matrix.iloc[row, col]:.1f}", ha="center", va="center", color="white", fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_all(results_dir: str | Path = "results") -> None:
    results_dir = Path(results_dir)
    plot_gc_content(results_dir / "sequence_metrics.csv", results_dir / "gc_comparison.png")
    plot_identity_heatmap(results_dir / "sequence_identity_matrix.csv", results_dir / "identity_heatmap.png")
    print(f"Saved visualizations to {results_dir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot pipeline results.")
    parser.add_argument("--results-dir", default="results")
    args = parser.parse_args()
    plot_all(args.results_dir)
