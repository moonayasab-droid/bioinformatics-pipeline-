16S rRNA Reliability and Methodological Sensitivity Analysis
Overview & Research Question
This repository houses an independent computational research study investigating a central question in molecular microbial ecology: How reliably does 16S rRNA sequence similarity recover bacterial taxonomy, and how sensitive are those conclusions to methodological choices?
Rather than treating bioinformatics pipelines as black boxes, this project establishes a rigorous baseline using curated sequence data to evaluate clustering accuracy, distance thresholds, and the impact of quality control parameters on downstream phylogenetic reconstruction.
1. Research Design & Protocol
The investigation follows a controlled multi-phase computational design:
 * Core Topic: Assessment of taxonomic resolution, method comparison, and robustness against sequence anomalies.
 * Dataset Protocol (Dataset V1.0): Strict inclusion and exclusion criteria were established prior to analysis:
   * Minimum Length: \ge 1400 bp to ensure near-full-length 16S rRNA coverage.
   * Ambiguity Threshold: Maximum allowable fraction of ambiguous bases (N \le 1.0\%).
   * Normalization: Strict handling of valid nucleotide bases (ACGTN) with systematic sequence cleanup.
2. Engineering & Testing Foundation
The pipeline is built using a modular, test-driven engineering approach to ensure reproducibility and code integrity:
 * Environment & Libraries: Python-based architecture utilizing Biopython for sequence parsing, alignment, and distance matrix computations.
 * Modular Structure:
   * src/pipeline.py: Orchestrates end-to-end execution.
   * src/qc.py: Houses sequence length verification and ambiguity filtering logic (passes_qc).
  * src/utils.py: Manages sequence normalization routines.
 * Automated Testing: Comprehensive unit test suites (tests/test_qc.py, tests/test_pipeline.py) verified via pytest to guarantee robust validation logic.
3. Materials & Methodology
 * Data Ingestion & Normalization: Raw candidate sequences (data/raw/candidate_sequences.fasta) are parsed via Biopython, cleaned, and standardized.
 * Quality Control Execution: The quality control script evaluates sequences against the Dataset V1.0 protocol, generating an official audit trail (results/dataset_qc_report.csv) that categorizes accessions into INCLUDE or EXCLUDE status.
 * Distance Calculation & Visualization: Filtered sequences undergo pairwise alignment to generate identity matrices (sequence_identity_matrix.csv), identity heatmaps (results/identity_heatmap.png), and Neighbor-Joining phylogenetic trees (results/neighbor_joining_tree.nwk).
4. Results & Baseline Findings
 * QC Yield: Evaluation of 87 candidate sequences demonstrated high data fidelity, with 86 sequences successfully passing all quality control filters and 1 sequence flagged for exclusion due to formatting or ambiguity anomalies.
 * Clustering & Distance Distribution: Pairwise identity heatmaps revealed distinct boundaries between expected taxonomic groupings, confirming that full-length similarity metrics reliably delineate high-level clades.
 * Phylogenetic Recovery: Neighbor-Joining tree construction successfully clustered major lineages into monophyletic groups under standard alignment parameters, establishing a stable control baseline for ongoing sensitivity and failure-mode analyses.
Repository Structure
bioinformatics-pipeline/
├── data/
│   └── raw/               # Candidate FASTA sequences
├── src/
│   ├── pipeline.py        # Core pipeline orchestrator
│   ├── qc.py              # Quality control and validation rules
│   └── utils.py           # Sequence normalisation utilities
├── tests/
│   ├── test_qc.py         # Unit tests for quality control
│   └── test_pipeline.py   # Unit tests for pipeline execution
├── results/               # Generated heatmaps, trees, and audit reports
├── README.md              # Project documentation
└── requirements.txt       # Project dependencies
