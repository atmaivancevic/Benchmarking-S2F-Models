# Benchmarking-S2F-Models

## Overview

Sequence-to-function (S2F) models like AlphaGenome, Borzoi and ChromBPNet predict regulatory activity directly from DNA sequence. This workflow performs **virtual CRISPR deletions**, mimicking experimental knockouts *in silico*, and compares the predicted chromatin and RNA-seq changes against published CRISPR KO data for TE-derived enhancers, promoters, and exons.

---

# 1. AlphaGenome

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![AlphaGenome](https://img.shields.io/badge/AlphaGenome-v0.5.1-green.svg)](https://pypi.org/project/alphagenome/)

Automated Python workflow for benchmarking [AlphaGenome](https://github.com/google-deepmind/alphagenome) predictions against CRISPR-validated TE-derived regulatory elements.

---

## Requirements

- Python 3.10+
- An [AlphaGenome API key](https://www.alphagenomedocs.com/installation.html) (free for non-commercial use)

### Installation
```bash
python3 -m venv alphagenome_venv
source alphagenome_venv/bin/activate
pip install notebook alphagenome numpy pandas matplotlib
```
---

## Files

| File | Description |
|------|-------------|
| `alphagenome_automated_analysis.ipynb` | Main interactive notebook for exploring individual elements |
| `helper_functions.py` | Core functions used by the notebook (plot_baseline, plot_deletion, plot_difference, score_deletion, etc.) |
| `LTR10_variants.tab` | Input variant file: e.g. CRISPR-validated LTR10 enhancers & control deletions |
| `my_api_key.txt` | Your personal AlphaGenome API key (**not included**) |

---

## Input Format

The variant file is tab-delimited with the following columns:

| Column | Description |
|--------|-------------|
| `ID` | Element identifier (e.g. `LTR10.ATG12`) |
| `CHROM` | Chromosome (e.g. `chr5`) |
| `POS` | Start position (0-based) |
| `REF` | Reference sequence to be deleted (the full element sequence) |
| `ALT` | Alternate allele (`.` = deletion, i.e. empty string) |
| `Output` | 1 = experimentally validated enhancer, 0 = negative control |
| `Study_ID` | Source publication |
| `Study_Variant_ID` | Element ID in the original study |

---

## Usage

### Interactive notebook

```bash
cd AlphaGenome
source alphagenome_venv/bin/activate
jupyter notebook
```

Open `alphagenome_automated_analysis.ipynb` and edit the **Configuration Section**:

```python
API_KEY_FILE  = 'my_api_key.txt'
VARIANT_FILE  = 'LTR10_variants.tab'
PRIMARY_ASSAYS       = ['RNA_SEQ', 'ATAC', 'CHIP_HISTONE']
PRIMARY_CELL_LINE    = 'HCT116'
PRIMARY_CELL_LINE_ID = 'EFO:0002824'
```

### Helper functions

```python
plot_baseline('LTR10.ATG12')                    # baseline genome state
plot_deletion('LTR10.ATG12')                    # predicted state after virtual KO
plot_difference('LTR10.ATG12')                  # difference: KO minus baseline
plot_difference_corrected_rna('LTR10.ATG12')    # RNA difference with frameshift correction
plot_overlaid_corrected_rna('LTR10.ATG12')      # reference vs deletion RNA overlay
score_deletion('LTR10.ATG12', gene_name='ATG12', assay='RNA_SEQ')  # quantile scores
```

---

## License

MIT License - see [LICENSE](LICENSE) file.

## Authors

Atma Ivancevic, PhD  
BioFrontiers Institute, University of Colorado Boulder
