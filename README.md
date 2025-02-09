# GC Content and ENC Analysis

This repository contains Python scripts for analyzing GC content and codon usage in genetic sequences.

## Scripts

### 1. `gc_content_calculator.py`
- Reads a FASTA file and calculates GC content at different codon positions (GC1, GC2, GC3, GC3s, GC12).
- Saves the results to `gc_content_results.csv`.

### 2. `enc_analysis.py`
- Reads `codonw_results.csv` and calculates expected ENC (`ENC_exp`).
- Computes the relative difference between observed and expected ENC values.
- Saves the results to `ENC_exp_ENC_obs_result.csv`.

## Usage
Run the scripts with:
```sh
python gc_content_calculator.py
python enc_analysis.py
```

## Requirements

- Python 3.x
- Biopython (`pip install biopython`)
- Pandas (`pip install pandas`)

## Contact

For questions, contact [bio_67] at [bio_67@163.com].