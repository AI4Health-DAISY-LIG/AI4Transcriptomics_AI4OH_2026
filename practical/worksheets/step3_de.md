# Step 3: differential expression with PyDESeq2

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

## Design

Input matrix (which file, which filter): 
Design formula: 
Covariates included and why, or why not: 
Covariates you wanted and could not include, and why: 
Number of samples per group: 

## Results

| Comparison | Genes tested | padj < 0.05 | of which up | of which down |
|------------|--------------|-------------|-------------|---------------|
| FSHD vs Control | | | | |
| Group 1 vs Control | | | | |
| Group 2 vs Control | | | | |
| Group 3 vs Control | | | | |
| Group 4 vs Control | | | | |

Figures: MA plot, volcano plot, p-value histogram (in `results/step3_de/`).

Top 10 genes by adjusted p-value for one comparison of your choice, with log2
fold change:

## Checks

- Describe the p-value histogram: a flat part with a spike near 0, or another shape?
- Run the analysis with and without batch in the design. Report how many genes
  are significant in both, and in one only.
- Group 1 biopsies show almost no expression of the biomarkers. Count the genes
  that differ between group 1 and the controls, and list the explanations you can
  think of other than the disease.

## Decisions I made

-
