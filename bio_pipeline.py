import argparse
import csv
from pathlib import Path
from typing import Dict, Mapping, Sequence

import matplotlib.pyplot as plt
from Bio import Phylo, SeqIO
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix, DistanceTreeConstructor


def _normalise_sequence(sequence: str) -> str:
    """Return an uppercase sequence with whitespace removed."""
    if not isinstance(sequence, str):
        raise TypeError("Sequences must be strings.")
    normalised = "".join(sequence.split()).upper()
    if not normalised:
        raise ValueError("Sequences must not be empty.")
    return normalised


def _make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    return aligner


def pairwise_sequence_alignment(sequence_a: str, sequence_b: str) -> dict:
    """Return the best global alignment and its percentage identity."""
    sequence_a = _normalise_sequence(sequence_a)
    sequence_b = _normalise_sequence(sequence_b)
    alignment = _make_aligner().align(sequence_a, sequence_b)[0]
    aligned_a, aligned_b = str(alignment[0]), str(alignment[1])
    alignment_length = len(aligned_a)
    matches = sum(a == b for a, b in zip(aligned_a, aligned_b))
    return {
        "sequence_a_aligned": aligned_a,
        "sequence_b_aligned": aligned_b,
        "matches": matches,
        "alignment_length": alignment_length,
        "identity_percent": (matches / alignment_length * 100) if alignment_length else 0.0,
        "score": float(alignment.score),
    }


def pairwise_sequence_identity(sequence_a: str, sequence_b: str) -> float:
    """Return global alignment identity as a percentage."""
    return pairwise_sequence_alignment(sequence_a, sequence_b)["identity_percent"]


def calculate_identity_matrix(sequences: Mapping[str, str]) -> Dict[str, Dict[str, float]]:
    """Compare every sequence to every other sequence with PairwiseAligner."""
    if len(sequences) < 2:
        raise ValueError("At least two sequences are required to build an identity matrix.")
    names = list(sequences)
    cleaned = {name: _normalise_sequence(seq) for name, seq in sequences.items()}
    matrix = {row: {} for row in names}
    for i, name_a in enumerate(names):
        for j in range(i, len(names)):
            name_b = names[j]
            identity = 100.0 if i == j else pairwise_sequence_identity(cleaned[name_a], cleaned[name_b])
            matrix[name_a][name_b] = identity
            matrix[name_b][name_a] = identity
    return matrix


def load_fasta_sequences(fasta_path: str | Path) -> Dict[str, str]:
    """Load named nucleotide sequences from a FASTA file."""
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}.")
    sequences = {record.id: _normalise_sequence(str(record.seq)) for record in records}
    if len(sequences) != len(records):
        raise ValueError("FASTA record identifiers must be unique.")
    return sequences


def write_identity_matrix(matrix: Mapping[str, Mapping[str, float]], output_path: str | Path) -> None:
    names = list(matrix)
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence"] + names)
        for name in names:
            writer.writerow([name] + [f"{matrix[name][other]:.6f}" for other in names])


def build_neighbor_joining_tree(matrix: Mapping[str, Mapping[str, float]]):
    """Build a Neighbor-Joining tree from percentage identities."""
    names = list(matrix)
    distances = [[0.0 if i == j else (100.0 - matrix[names[i]][names[j]]) / 100.0
                  for j in range(i + 1)] for i in range(len(names))]
    distance_matrix = DistanceMatrix(names, distances)
    return DistanceTreeConstructor().nj(distance_matrix)


def write_neighbor_joining_tree(matrix: Mapping[str, Mapping[str, float]], output_path: str | Path) -> None:
    tree = build_neighbor_joining_tree(matrix)
    Phylo.write(tree, str(output_path), "newick")


def analyze_sequence(seq: str) -> dict:
    seq = _normalise_sequence(seq)
    counts = {base: seq.count(base) for base in "ATGC"}
    return {
        "length": len(seq),
        "A": counts["A"], "T": counts["T"], "G": counts["G"], "C": counts["C"],
        "GC_content": (counts["G"] + counts["C"]) / len(seq) * 100,
    }


def analyze_fasta(fasta_path: str | Path, output_dir: str | Path = "results") -> dict:
    """Run metrics, all-vs-all alignment, and NJ tree generation for a FASTA file."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sequences = load_fasta_sequences(fasta_path)
    if len(sequences) < 2:
        raise ValueError("The FASTA file must contain at least two sequences.")

    matrix = calculate_identity_matrix(sequences)
    write_identity_matrix(matrix, output_dir / "sequence_identity_matrix.csv")
    write_neighbor_joining_tree(matrix, output_dir / "neighbor_joining_tree.nwk")

    with open(output_dir / "sequence_metrics.csv", "w", newline="") as handle:
        fieldnames = ["organism", "length", "gc_percentage"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for name, sequence in sequences.items():
            metrics = analyze_sequence(sequence)
            writer.writerow({"organism": name, "length": metrics["length"],
                             "gc_percentage": f"{metrics['GC_content']:.6f}"})
    return matrix


def fetch_ncbi_gene_data() -> str:
    """Return the legacy demonstration sequence used by the original pipeline."""
    return "ATCGATCGATCGATCGATCGGCGCGCATATCGATCGATCGATCGATCGGCGCGCAT"


def plot_nucleotide_counts(results_dict: dict, output_path: str | Path = "genomic_analysis_chart.png") -> None:
    bases = ["Adenine (A)", "Thymine (T)", "Guanine (G)", "Cytosine (C)"]
    counts = [results_dict["A"], results_dict["T"], results_dict["G"], results_dict["C"]]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(bases, counts, color=["#ff9999", "#66b3ff", "#99ff99", "#ffcc99"])
    ax.set(xlabel="Nucleotide Base", ylabel="Count", title="Nucleotide Distribution")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare bacterial 16S rRNA sequences.")
    parser.add_argument("--fasta", required=True, help="Input FASTA containing 16S rRNA sequences")
    parser.add_argument("--output-dir", default="results", help="Directory for CSV, Newick, and plots")
    args = parser.parse_args()
    matrix = analyze_fasta(args.fasta, args.output_dir)
    print(f"Compared {len(matrix)} sequences.")
    print(f"Wrote identity matrix and Neighbor-Joining tree to {args.output_dir}/")


if __name__ == "__main__":
    main()
