# bioinformatics-pipeline-

# Comparative Genomic Pipeline: 16S rRNA GC-Content Analysis

## Overview
This repository contains a reproducible Python bioinformatics pipeline developed to retrieve, parse, analyze, and visualize nucleotide sequences across multiple bacterial species. By focusing on the **16S rRNA gene**, this project investigates how genomic metrics—specifically Guanine-Cytosine (GC) content and sequence length—vary across diverse prokaryotic lineages.

## Research Question
> *How does the GC content and nucleotide composition of the 16S rRNA gene vary across representative bacterial species (*E. coli, B. subtilis, P. aeruginosa, S. aureus, S. coelicolor*)?*

## Project Structure
```text
bioinformatics-pipeline/
├── results/
│   ├── sequence_metrics.csv       # Extracted quantitative sequence metrics
│   └── gc_comparison.png          # Comparative GC-content visualization chart
├── src/ (or root scripts)
│   ├── bio_pipeline.py            # NCBI Entrez data fetcher & metric extractor
│   └── visualize_results.py       # Matplotlib comparative visualization generator
├── requirements.txt               # Project dependencies
└── README.md                      # Project documentation
