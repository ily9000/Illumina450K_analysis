Methylation, the addition of a methyl group to a DNA nucleotide, plays a role in psychiatric disorders. To detect which nucleotides (sites) in our genome play a role in methylation, we use the Illumina 450K array to quantify methylation status of 480,000 sites across the genome. In a given cell, each site can assume 3 states: unmethylated,  hemi-methylated, or methylated. Since the array measures methylation on a collection of cells, it reports a continous value between 0 (completely unmethylated) and 1 (completely methylated) for each site.

# Illumina450K_analysis
Analysis for epigenetic data from Illumina 450K 

1. runbeads.R is used to transform the raw IDAT image files from the methylation array into a R object using the RnBeads package.
2. chk_qc/pca_qc.Rnw runs a PCA on the 600 negative controls in the green and red channel (Cy3 & Cy5)
3. chk_qc/chk_cntrl/chk_cntrol.Rnw
In addition to the negative controls, there are 6 categories of quality control probes in each sample to ensure each step of the methylation assay was completed. The probes are categorized: Staining, Specificity I, Specificity II, Bisulfite converation I, Bisulfite converation II, Extension. Samples which did not have the expected values for all the probes should be removed prior to analysis and are output to rmsamps.csv and are also listed in chk_cntrol.pdf
