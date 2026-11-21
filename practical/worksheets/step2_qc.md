# Step 2: quality control and exploration

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

## Filtering

| Filter | Threshold I chose | Why | Samples and genes before | after |
|--------|-------------------|-----|--------------------------|-------|
| | | | | |

Does the paper state its thresholds? Quote them if it does, say so if it does not.
Do the four biomarker genes of the paper pass your filter? 

## Exploration

Figures (in `results/step2_qc/`):

- PCA coloured by condition: 
- PCA coloured by each technical variable: 
- PCA coloured by the authors' group (`rnaseq.rank`): 

Name the variable that separates the samples on PC1, then on PC2: the condition,
a technical variable, the authors' group, or none you can name. Give the variance
PC1 and PC2 explain.

Compare your PCA with Figure 3 of the paper. The same or not, and why:

## Samples I excluded, and why

-

## Decisions I made
1.- I matched count columns to GEO metadata using `title` instead of `char:sample_id`, because all 43 samples matched by title, while one sample (`01-0022-1`) did not match by sample_id. This allowed us to retain all 43 samples without dropping any.

-
