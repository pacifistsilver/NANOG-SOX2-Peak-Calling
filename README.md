# Extract NANOG/SOX2 Sites vs Bursting Kinetics

This repository contains the analysis scripts used to cross-reference TF binding data (NANOG and SOX2) with transcriptional bursting kinetics.

## Pipeline Overview

1. **Peak Calling**: `run_macs2.pbs` pulls the original Chen *et al.* (2008) ChIP-Seq read alignments from GEO and runs MACS2 peak calling on the Imperial HPC.
2. **Site Counting & Liftover**: `count_sites_chen2008.py` takes the MACS2 `.narrowPeak` output (mm8), lifts the coordinates over to `mm9`, identifies active promoters (using H3K4me3 marks), and counts the discrete number of TF binding sites at each active promoter (`TSS +/- 5kb`).
3. **Correlation Plotting**: `plot_binding_vs_bursting.py` merges the site counts with the Ochiai *et al.* (2020) single-cell transcriptomics bursting kinetic data (Serum+LIF) and generates publication-ready plots.

## Results

We demonstrate a strong positive correlation between the number of pluripotency factor binding sites (SOX2/NANOG) and **Burst Frequency**, with negligible impact on **Burst Size**, confirming enhancer-promoter communication models.
