# LTR10_variants

Data files supporting the analysis of LTR10 transposable elements as regulatory variants in the hg38 human genome, with a focus on LTR10A and LTR10F subfamilies in HCT116 colorectal cancer cells.

---

## Files

| File | n | Description |
|------|---|-------------|
| `2331_LTR10_elements_annotated_in_hg38.bed` | 2,331 | All LTR10 elements (subfamilies A–G) annotated in hg38 |
| `143_LTR10_elements_ranked_by_peak_signal.txt` | 143 | LTR10A/F elements ranked by FOSL1/H3K27ac signal in HCT116 cells |
| `71_LTR10_enhancers_predicted_by_ABC_model.txt` | 71 | LTR10A/F elements predicted as enhancers by FOSL1/MAPK silencing + ABC model |
| `6_CRISPRd_LTR10_sequences.txt` | 6 | Coordinates and sequences of CRISPR-deleted LTR10 enhancers |
| `15_LTR10_known_promoters.txt` | 15 | Known LTR10A/F promoters in hg38 |
| `DESeq2_tables/` | — | DESeq2 results from all 6 CRISPR deletion experiments (`.xlsx`) |

---

## File details

**`143_LTR10_elements_ranked_by_peak_signal.txt`**  
Ranked by FOSL1 ChIP-seq peak score. Based on in-house H3K27ac CUT&RUN and public FOSL1 ChIP-seq (Cistrome ID 46212). Restricted to LTR10A and LTR10F subfamilies.

**`71_LTR10_enhancers_predicted_by_ABC_model.txt`**  
Candidate enhancers identified through FOSL1 and MAPK silencing experiments combined with the [Activity-By-Contact (ABC) model](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction). 8 of these elements were selected for CRISPR deletion (6 total sequences after merging).

**`6_CRISPRd_LTR10_sequences.txt`**  
The 6 LTR10 sequences selected for CRISPR deletion. Note that some entries (e.g., `LTR10.ATG12`, `LTR10.XRCC4`) span two adjacent LTR10 elements.

**`15_LTR10_known_promoters.txt`**  
Curated set of LTR10A/F elements with evidence of promoter activity in hg38.

**`DESeq2_tables/`**  
One table per CRISPR experiment (6 total). Each `.xlsx` file contains differential expression results comparing CRISPR deletion to control.

---

## Bigwigs

Raw fastq files, bigwigs and count tables are available on GEO dataset [GSE186618](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5657882). Bigwigs have also been downloaded to the shared cluster dir (chuong-layer).
