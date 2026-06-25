# EMT, proteoglycan core proteins and GAG biosynthetic enzymes
This repository contains R scripts "Epithelial-mesenchymal transition in embryonic development has a major impact on proteoglycan core proteins expression in adult human tissue".

The scripts document the workflow used to generate the processed analyses, summary tables, and figure panels reported in the manuscript.

## Repository structure
1. Adult tissue expression analyses
  - 1a_Proteoglycan_heatmap_related.R
    - Generates heatmaps and summary analyses for proteoglycan core protein expression across adult human cell types.
  - 1b_GAG_heatmap_related.R
    - Generates heatmaps and summary analyses for GAG biosynthetic enzyme expression across adult human cell types.
  - 1c_Expression_of_Subgroups.R
    - Summarizes expression patterns of predefined proteoglycan and GAG-related gene subgroups.
2. Embryonic development analyses
  - 2a_Embryonic_Data_Processing.R
    - Processes embryonic single-cell RNA-seq data for downstream analyses.
  - 2b_Embryonic_Lineage_Processing.R
    - Prepares lineage-specific embryonic data and expression summaries.
  - 2c_Association_towards_germlayers.R
    - Tests associations between proteoglycan/GAG gene expression and embryonic germ-layer identity.
  - 2d_Embryo_adult_sankey.R
    - Generates summaries linking embryonic expression patterns to adult cell-type expression patterns.
3. EMT transcription-factor analyses
  - 3a_Adult_Tissue_Correlation_0.4.R
    - Performs correlation analyses between EMT transcription factors and proteoglycan/GAG genes in adult expression data.
  - 3b_Bubble_Plot_for_Correlation.R
    - Generates bubble plots summarizing EMT transcription-factor correlations.
  - 3c_Representative_TF_perturbation_loop.R
    - Summarizes EMT transcription-factor perturbation datasets.
  - 3c_1_Representative_RNAseq_DE.R
    - Provides a representative RNA-seq differential-expression workflow.
  - 3c_2_Representative_Affymetrix_CEL_oligo_RMA_DE.R
    - Provides a representative Affymetrix microarray differential-expression workflow.
  - 3c_3_Representative_Illumina_array_DE.R
    - Provides a representative Illumina microarray differential-expression workflow.
  - 3c_4_Representative_Agilent_single_channel_array_DE.R
    - Provides a representative Agilent single-channel microarray differential-expression workflow.
  - 3d_TF_Making_Excel.R
    - Compiles EMT transcription-factor perturbation results into summary Excel files.
  - 3e_Bubble_Plots_for_TF.R
    - Generates bubble plots of EMT transcription-factor perturbation effects.
  - 3f_Box_Plots_for_TF.R
    - Generates box plots of EMT transcription-factor perturbation effects.
4. ChIP-seq and regulatory-overlap analyses
  - 4a_Finding_Overlap.R
    - Identifies overlaps between EMT transcription-factor ChIP-seq peaks and GeneHancer regulatory regions linked to proteoglycan/GAG genes.
  - 4b_Chipseq_Getting_summary_files.R
    - Generates summary files from ChIP-seq/GeneHancer overlap results.
  - 4c_Chipseq_Combined_summary_files.R
    - Combines ChIP-seq/GeneHancer overlap summaries into final tables.

**Input data**\
The scripts use processed and publicly available datasets, including adult single-cell expression data, human gastrulation single-cell RNA-seq data, EMT transcription-factor perturbation datasets, ChIP-seq peak files, and GeneHancer regulatory annotations.

Some input files are not included in this repository because they originate from public databases or external repositories. Processed summary files associated with the manuscript will be provided through the accompanying Zenodo record (*will be updated when available*).

**The scripts are provided to document the analyses performed for the manuscript and to support transparency and reproducibility.**

Citation
If you use or refer to these scripts, please cite the associated manuscript and Zenodo release (*will be updated when available*).
