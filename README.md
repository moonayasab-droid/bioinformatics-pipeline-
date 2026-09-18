# Comparative 16S rRNA Bioinformatics Pipeline

A reproducible Python workflow for comparing bacterial 16S rRNA sequences. The pipeline calculates nucleotide metrics, performs an all-vs-all global alignment with Biopython's `PairwiseAligner`, creates a percentage identity matrix, generates a publication-ready heatmap, and builds a Neighbor-Joining tree from alignment distances.

## Workflow

```text
FASTA sequences
      |
      v
Biopython PairwiseAligner ---> identity matrix (CSV) ---> identity heatmap (PNG)
      |
      +--> sequence metrics (CSV) ---------------------> GC-content chart (PNG)
      |
      +--> 1 - identity distance -----------------------> Neighbor-Joining tree (Newick)
```

## Installation

Python 3.10+ is recommended.

```bash
python -m venv .venv
# macOS/Linux: source .venv/bin/activate
# Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
```

## Sample input

A sample FASTA is included in `data/16s_sequences.fasta` so the pipeline can run immediately.

```fasta
>Escherichia_coli
ATCGATCGATCGATCGATCGGCGCGCATATCGATCGATCGATCGATCGGCGCGCAT
>Salmonella_enterica
ATCGATCGATCGATCGATCGGCGCGCATATCGATCGATCGATCGATCGGCGCGCAT
>Bacillus_subtilis
ATCGATCGATCGATCCGATCGGCGCGCATATCGATCGATCGATCGGCGCGCAT
>Pseudomonas_aeruginosa
ATCGATCGATCGATCGATCGGCGCGCAGATCGATCGATCGATCGATCGGCGCGCAT
```

## Run the analysis

```bash
python bio_pipeline.py --fasta data/16s_sequences.fasta --output-dir results
```

This creates:

- `results/sequence_metrics.csv` — sequence length and GC content
- `results/sequence_identity_matrix.csv` — all-vs-all identity percentages
- `results/neighbor_joining_tree.nwk` — a Neighbor-Joining tree in Newick format

Generate the visualizations:

```bash
python visualize_results.py --results-dir results
```

This creates:

- `results/gc_comparison.png` — GC-content comparison chart
- `results/identity_heatmap.png` — pairwise identity heatmap

The Newick tree can be opened in FigTree, iTOL, or rendered with Biopython.

## Method notes

Identity is calculated from a global alignment as matching alignment columns divided by total alignment columns, including gaps. Distances for the Neighbor-Joining tree are `1 - identity` represented as a fraction. This is a simple exploratory phylogenetic workflow, not a substitute for curated multiple-sequence alignment, model selection, or bootstrap support.

## Repository layout

```text
.
├── bio_pipeline.py       # Metrics, pairwise alignment, matrix, and tree generation
├── visualize_results.py  # GC-content and identity heatmap plots
├── requirements.txt
├── accession_metadate.csv
├── data/
│   └── 16s_sequences.fasta
├── results/              # Generated outputs
├── README.md
└── genomic_analysis_chart.png
```

## Reproducibility and quality checks

- Sequence names must be unique in the FASTA file.
- At least two sequences are required.
- Empty sequences are rejected.
- Matrix diagonal values are 100%.
- The matrix is symmetric by construction.
- Generated files are written to the selected output directory.

## License

No license has been declared yet. Add a license before redistributing the project.
