Comparative 16S rRNA Bioinformatics Pipeline
A reproducible Python-based computational biology workflow for analyzing bacterial 16S rRNA sequences and investigating the taxonomic resolution and methodological robustness of 16S rRNA phylogenetic inference.

The repository currently contains a baseline sequence-analysis pipeline and is being developed into a controlled computational research study examining how reliably 16S rRNA sequence data recover established bacterial taxonomic relationships and how sensitive those inferences are to methodological choices.

Research Question
How reliably can full-length 16S rRNA sequence analysis recover bacterial taxonomic relationships, and how sensitive are those inferences to methodological choices and taxonomic resolution?

The study investigates three related aspects of 16S-based inference:

Taxonomic resolution — How well does 16S rRNA distinguish organisms at family, genus, and species levels?
Methodological sensitivity — Do reasonable choices of alignment, distance calculation, and phylogenetic inference alter the resulting relationships?
Robustness and limitations — Which relationships remain stable, and where does 16S-based inference produce ambiguous or conflicting results?
The complete research methodology is documented in:

research/protocol_v1.1.md

Project Overview
16S rRNA sequencing is widely used for bacterial identification and phylogenetic analysis because the gene contains conserved and variable regions that provide information about evolutionary relationships.

However, 16S rRNA is a single genetic marker and does not necessarily provide sufficient resolution for every bacterial lineage or taxonomic level. Closely related organisms may have highly similar 16S sequences, while relationships inferred from a single marker can differ from those supported by broader genomic evidence.

Rather than assuming that 16S-based inference is either universally reliable or universally unreliable, this project investigates where the method performs well, where its resolution decreases, and how methodological decisions affect its conclusions.

The computational pipeline serves as the experimental framework for this investigation.

Research Design
The study is organized around a baseline analysis followed by controlled experiments.

                    Curated 16S Dataset
                            │
                            ▼
                     Sequence QC
                            │
                            ▼
                    Baseline Analysis
                            │
             ┌──────────────┼──────────────┐
             ▼              ▼              ▼
       Taxonomic        Methodological   Robustness
       Resolution        Sensitivity      Analysis
             │              │              │
             └──────────────┼──────────────┘
                            ▼
                     Failure Cases
                            │
                            ▼
                  Genomic Validation
                            │
                            ▼
                    Statistical Analysis
                            │
                            ▼
                       Conclusions
The research protocol is designed so that reference taxonomy is used for evaluation, rather than for selecting parameters or optimizing the pipeline.

Current Baseline Pipeline
The current implementation provides the foundation for the experimental study.

It performs:

FASTA sequence loading
Sequence normalization and validation
All-vs-all pairwise global alignment
Pairwise sequence identity calculation
Identity matrix generation
Distance-matrix construction
Neighbor-Joining phylogenetic inference
Sequence-length calculation
GC-content calculation
CSV result generation
Heatmap and comparative visualization
Baseline workflow
FASTA sequences
      │
      ▼
Sequence normalization
      │
      ▼
Pairwise global alignment
      │
      ▼
Sequence identity matrix
      │
      ├──────────────► Identity heatmap
      │
      ▼
Distance matrix
      │
      ▼
Neighbor-Joining tree
      │
      └──────────────► Newick output

Sequence data
      │
      ├──────────────► Length metrics
      │
      └──────────────► GC-content metrics
                              │
                              ▼
                       Comparative chart
Scientific Hypotheses
The research protocol currently evaluates two primary hypotheses.

H1 — Hierarchical Taxonomic Resolution
Agreement between 16S-derived relationships and reference taxonomy will decrease as the required taxonomic resolution becomes finer, with stronger agreement expected at broader taxonomic levels than at the species level.

H2 — Methodological Sensitivity
Phylogenetic tree topologies and distance-based clustering will be sensitive to alternative methodological choices, particularly among closely related taxa within the same genus.

These hypotheses are treated as testable propositions, not predetermined conclusions. Results that contradict either hypothesis will be reported and investigated.

Planned Experiments
Experiment 1 — Taxonomic Resolution
The analysis will evaluate correspondence between inferred sequence relationships and reference taxonomy at multiple taxonomic levels:

Family
Genus
Species
Multiple complementary measures will be considered, including:

Monophyly
Pairwise sequence-identity distributions
Cluster purity
Tree/reference concordance
The purpose is to determine whether the resolution of 16S-based inference changes systematically across taxonomic levels.

Experiment 2 — Methodological Sensitivity
The same curated dataset will be analyzed under controlled methodological variations.

Only one methodological factor will be changed at a time wherever possible.

Distance models
Potential comparisons include:

Identity-derived distance
Jukes-Cantor
Kimura 2-parameter
Tree inference
Potential comparisons include:

Neighbor-Joining
UPGMA
Maximum Likelihood
The assumptions and input requirements of each approach will be documented rather than treating the algorithms as interchangeable.

Alignment sensitivity
The study will investigate the effect of alignment methodology and relevant scoring parameters, including:

Match/mismatch scores
Gap penalties
Alignment algorithms
Experiment 3 — Robustness and Statistical Support
Phylogenetic robustness will be investigated using resampling and tree-comparison approaches.

Bootstrap analysis
Bootstrap resampling will be used to determine whether important inferred relationships remain supported across resampled datasets.

Tree topology comparison
Trees generated using different methodological approaches will be compared using normalized Robinson-Foulds distance and complementary branch-support information where appropriate.

Experiment 4 — Failure-Case Analysis
The study will explicitly investigate cases in which sequence-based inference conflicts with reference taxonomy.

Candidate cases may include:

Very high 16S similarity between differently annotated taxa
Unexpected cross-genus clustering
Poorly supported branches
Relationships that change substantially across methodological approaches
These cases will be investigated using sequence characteristics, alignment behavior, metadata, and relevant scientific literature.

The purpose is not to remove inconvenient observations, but to understand why and where the method may lose resolution or robustness.

Experiment 5 — Independent Genomic Validation
For a predefined subset of organisms with suitable genomic data, 16S-derived relationships will be compared with broader genomic evidence.

Potential reference approaches include:

Core-genome phylogenies
Multi-locus genomic comparisons
Average Nucleotide Identity (ANI)
This analysis will investigate whether relationships inferred from a single conserved marker agree with relationships supported by broader genomic information.

The genomic validation subset will be selected using predefined criteria rather than based on whether the resulting comparison supports the 16S analysis.

Dataset Design
The evaluation dataset will be constructed using predefined selection criteria.

Inclusion Criteria
Sequences should:

Represent bacterial 16S rRNA
Contain at least 1,400 aligned nucleotide positions corresponding to the bacterial 16S rRNA locus
Have reliable accession records
Have sufficiently curated taxonomic metadata
Exclusion Criteria
Sequences may be excluded when they:

Fail the predefined completeness requirement
Contain excessive ambiguous nucleotide positions
Are recognized as chimeric
Represent redundant duplicate records of the same isolate
Sampling Strategy
The dataset will use a hierarchical sampling design containing:

Multiple bacterial families
Multiple genera within selected families
Multiple species within selected genera where available
The objective is to prevent a single dominant lineage from disproportionately determining global results.

The final dataset will be frozen according to the predefined selection criteria before hypothesis-testing analyses are performed.

Metadata
Each sequence will be associated with standardized metadata.

Expected fields include:

sequence_id
accession
organism
species
genus
family
order
class
phylum
sequence_length
GC_percentage
source
Additional metadata may be included when necessary for particular analyses.

Analysis Independence
Reference taxonomy will be used for evaluation rather than optimization.

Methodological parameters will be:

specified before analysis where possible, or
selected using criteria independent of reference-taxonomy agreement.
Sequences will not be removed simply because they produce unexpected or inconvenient results.

Any post hoc methodological changes will be documented and distinguished from the predefined analysis.

Statistical Analysis
The study will account for the non-independence of pairwise sequence comparisons.

Pairwise observations will not automatically be treated as independent biological observations because multiple comparisons can involve the same sequences and taxonomic groups.

Depending on the final dataset structure, analyses may use:

Group-level summaries
Permutation procedures
Bootstrap confidence intervals
Hierarchical statistical approaches
Appropriate parametric or non-parametric tests
Effect sizes and confidence intervals will be reported where appropriate.

Multiple comparisons will be addressed using appropriate correction procedures when required.

Null Model
An appropriate randomized null model will be used to establish the amount of taxonomic clustering expected by chance.

The planned framework will preserve the observed number and size distribution of taxonomic groups while randomly permuting taxonomic labels relative to sequence-derived relationships.

The exact permutation procedure will be specified before hypothesis testing.

Reproducibility
Reproducibility is treated as part of the research methodology rather than an optional software feature.

The repository will record:

Python version
Dependency versions
Dataset accession information
Dataset-selection criteria
Preprocessing procedures
Analysis parameters
Random seeds where applicable
Code used to generate tables and figures
Commands required to reproduce the analysis
The goal is for an independent researcher to be able to reconstruct the computational analysis from the documented inputs and procedures.

Repository Structure
The repository is organized around both the software and the research workflow.

comparative-16s-pipeline/
│
├── README.md
├── requirements.txt
├── LICENSE
│
├── research/
│   ├── protocol_v1.1.md
│   ├── dataset_design.md
│   ├── analysis_plan.md
│   └── references.md
│
├── data/
│   ├── raw/
│   └── processed/
│
├── src/
│   └── ...
│
├── tests/
│   └── ...
│
├── results/
│   └── ...
│
└── figures/
    └── ...
The exact structure may evolve as the research pipeline becomes more modular.

Baseline Implementation
The current baseline implementation contains:

bio_pipeline.py
Core sequence-analysis pipeline responsible for:

FASTA parsing
Sequence normalization
Pairwise alignment
Sequence identity calculation
Distance-matrix construction
Neighbor-Joining tree generation
Sequence metrics
CSV output
visualize_results.py
Visualization workflow responsible for generating:

Sequence identity heatmaps
GC-content comparative plots
Additional result visualizations as the research analysis expands
requirements.txt
Contains the Python dependencies required for the current computational workflow.

Installation
Requirements
Python 3.10 or newer
pip
Git
Clone the repository and enter the project directory:

git clone <REPOSITORY_URL>
cd comparative-16s-pipeline
Create a virtual environment:

macOS / Linux
python -m venv .venv
source .venv/bin/activate
Windows
python -m venv .venv
.venv\Scripts\activate
Install dependencies:

python -m pip install -r requirements.txt
Running the Baseline Pipeline
The baseline pipeline accepts a FASTA file and an output directory.

Example:

python bio_pipeline.py \
    --fasta data/example.fasta \
    --output-dir results
On Windows PowerShell:

python bio_pipeline.py `
    --fasta data/example.fasta `
    --output-dir results
The exact command-line interface may evolve as the research pipeline becomes modular.

Visualization
After generating the baseline results:

python visualize_results.py --results-dir results
The visualization workflow can generate outputs such as:

results/
├── sequence_identity_matrix.csv
├── sequence_metrics.csv
├── neighbor_joining_tree.nwk
├── identity_heatmap.png
└── gc_comparison.png
The exact filenames may change as the experimental pipeline develops.

Example Input
A minimal FASTA file follows the standard format:

>Escherichia_coli
ATG...
>Salmonella_enterica
ATG...
>Bacillus_subtilis
ATG...
>Pseudomonas_aeruginosa
ATG...
For the research dataset, sequences will be accompanied by standardized metadata and accession information rather than relying solely on FASTA descriptions.

Methodological Notes
Sequence Identity
The baseline implementation calculates sequence identity from global pairwise alignments.

Identity is calculated as the proportion of matching alignment columns relative to the total alignment length under the baseline definition.

Because gap treatment and alignment strategy can affect identity values, the baseline identity metric is treated as an analysis measure, not as a universal evolutionary distance.

This distinction is important for the methodological-sensitivity experiments.

Baseline Quality Checks
The current implementation includes validation checks such as:

Sequence names must be unique
At least two sequences are required
Empty sequences are rejected
Sequences are normalized before analysis
Identity-matrix diagonal values should equal 100%
The identity matrix should be symmetric
Output files should be written to the designated directory
Additional QC checks will be introduced as the research dataset is constructed.

Planned Research Outputs
The completed study is expected to produce structured tables and figures including:

Tables
Dataset composition
QC filtering summary
Taxonomic agreement statistics
Methodological comparison results
Genomic-validation comparisons
Figures
Sequence identity distributions by taxonomic tier
Baseline Neighbor-Joining tree
Taxonomically annotated tree with bootstrap support
Methodological comparison plots
Tree-topology distance matrices
16S versus genome-scale phylogenetic comparisons
Figures will be generated programmatically where possible to support reproducibility.

Limitations
The study explicitly recognizes several limitations.

Biological limitations
Multiple 16S rRNA operons may occur within a genome
Closely related organisms may have insufficient 16S variation for reliable discrimination
Recombination or horizontal transfer can complicate phylogenetic interpretation
Bacterial taxonomy is continually revised
Dataset limitations
Public databases may contain sampling biases
Some lineages may have substantially more available sequences than others
Metadata quality may vary between records
Methodological limitations
Alignment methods involve assumptions and parameter choices
Distance models represent evolutionary processes imperfectly
Phylogenetic algorithms have different assumptions
Robinson-Foulds distance captures certain aspects of topology but not every dimension of phylogenetic similarity
Evaluation limitations
Reference taxonomy is not an infallible ground truth.

A monophyletic group does not automatically establish an absolute species boundary, and a phylogenetic tree based on a single gene should not automatically be interpreted as the complete species or genome evolutionary history.

Scientific Decision Criteria
The project is designed to allow the empirical results to determine the conclusions.

If the analysis demonstrates that 16S rRNA:

resolves some taxonomic levels more reliably than others,
produces unstable relationships under particular methodological choices,
agrees with genome-scale evidence in some lineages but not others, or
produces identifiable failure cases,
those observations will be reported as scientific findings rather than treated as failures of the project.

The purpose of the study is therefore not to demonstrate that one method is universally correct.

The purpose is to characterize the resolution, robustness, and limitations of 16S rRNA-based phylogenetic inference under a controlled computational framework.

Research Documentation
The research methodology is documented separately from the software implementation.

Key documents:

research/
├── protocol_v1.1.md
├── dataset_design.md
├── analysis_plan.md
└── references.md
The protocol defines the scientific question, hypotheses, dataset criteria, experiments, evaluation strategy, statistical framework, and reproducibility requirements.

Project Status
Current stage: Baseline pipeline completed; research study in dataset-design and experimental-design phase.

The baseline software infrastructure is functional.

The next stage is to construct and freeze the curated evaluation dataset according to the research protocol before implementing the experimental extensions.

Planned development stages:

[✓] Baseline computational pipeline
[✓] Research protocol
[ ] Dataset design
[ ] Dataset acquisition and QC
[ ] Dataset freeze
[ ] Baseline benchmark
[ ] Taxonomic-resolution analysis
[ ] Methodological-sensitivity experiments
[ ] Robustness analysis
[ ] Failure-case investigation
[ ] Genomic validation
[ ] Statistical analysis
[ ] Final figures
[ ] Research report
Reproducibility Statement
This repository is intended to develop into a fully reproducible computational research workflow.

All major analytical decisions will be documented, and changes to the research protocol or computational methodology will be recorded through version control.

The final release will include the code, documented dataset sources/accessions, analysis configuration, generated results, and instructions necessary to reproduce the reported computational findings, subject to the redistribution terms of the underlying datasets.

License
This project is currently distributed under the license specified in the repository’s LICENSE file.

If no license has yet been selected, the project remains unlicensed until an explicit license is added.

Author
sara ahmed

Computational Biology / Bioinformatics Research Project

Citation
If this repository develops into a published research project, citation information will be added here.

For the current research protocol, see:

research/protocol_v1.1.md

