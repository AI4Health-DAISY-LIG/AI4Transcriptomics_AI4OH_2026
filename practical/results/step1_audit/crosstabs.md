# Cross-tabulations and statistical analysis

## Statistical methodology

- Categorical variables: Fisher's exact test for sparse 2x2 tables; otherwise Pearson's chi-square test.
- Continuous variables: Mann-Whitney U for two conditions; Kruskal-Wallis for more than two.
- Categorical effect size: bias-corrected Cramér's V.
- Missingness: separately tested against condition.
- Multiple testing: Benjamini-Hochberg FDR correction.
- Statistical significance: adjusted p-value < 0.05.

Missing observations are excluded from association tests but retained as a category in descriptive cross-tabs.

## char:batch

### Counts

| condition |   I |  II | III |  IV | All |
| :-------- | --: | --: | --: | --: | --: |
| Control   |   0 |   0 |   6 |   3 |   9 |
| FSHD      |   5 |   5 |  16 |   8 |  34 |
| All       |   5 |   5 |  22 |  11 |  43 |

### Row percentages (%)

| condition |    I |   II |  III |   IV |
| :-------- | ---: | ---: | ---: | ---: |
| Control   |    0 |    0 | 66.7 | 33.3 |
| FSHD      | 14.7 | 14.7 | 47.1 | 23.5 |

## instrument_model

### Counts

| condition | Illumina HiSeq 2500 | All |
| :-------- | ------------------: | --: |
| Control   |                   9 |   9 |
| FSHD      |                  34 |  34 |
| All       |                  43 |  43 |

### Row percentages (%)

| condition | Illumina HiSeq 2500 |
| :-------- | ------------------: |
| Control   |                 100 |
| FSHD      |                 100 |

## char:tissue

### Counts

| condition | muscle biopsy | All |
| :-------- | ------------: | --: |
| Control   |             9 |   9 |
| FSHD      |            34 |  34 |
| All       |            43 |  43 |

### Row percentages (%)

| condition | muscle biopsy |
| :-------- | ------------: |
| Control   |           100 |
| FSHD      |           100 |

## char:gender

### Counts

| condition | female | male | missing | All |
| :-------- | -----: | ---: | ------: | --: |
| Control   |      0 |    0 |       9 |   9 |
| FSHD      |     12 |   22 |       0 |  34 |
| All       |     12 |   22 |       9 |  43 |

### Row percentages (%)

| condition | female | male | missing |
| :-------- | -----: | ---: | ------: |
| Control   |      0 |    0 |     100 |
| FSHD      |   35.3 | 64.7 |       0 |

### Missingness by condition

| condition | missing | present |
| :-------- | ------: | ------: |
| Control   |       9 |       0 |
| FSHD      |       0 |      34 |

## char:age

### Counts

| condition | 30-44 | 45-59 | 60+ | <30 | missing | All |
| :-------- | ----: | ----: | --: | --: | ------: | --: |
| Control   |     0 |     0 |   0 |   0 |       9 |   9 |
| FSHD      |     5 |    14 |  13 |   2 |       0 |  34 |
| All       |     5 |    14 |  13 |   2 |       9 |  43 |

### Row percentages (%)

| condition | 30-44 | 45-59 |  60+ | <30 | missing |
| :-------- | ----: | ----: | ---: | --: | ------: |
| Control   |     0 |     0 |    0 |   0 |     100 |
| FSHD      |  14.7 |  41.2 | 38.2 | 5.9 |       0 |

### Missingness by condition

| condition | missing | present |
| :-------- | ------: | ------: |
| Control   |       9 |       0 |
| FSHD      |       0 |      34 |

## char:infer_gender

### Counts

| condition | Female | Male | All |
| :-------- | -----: | ---: | --: |
| Control   |      6 |    3 |   9 |
| FSHD      |     11 |   23 |  34 |
| All       |     17 |   26 |  43 |

### Row percentages (%)

| condition | Female | Male |
| :-------- | -----: | ---: |
| Control   |   66.7 | 33.3 |
| FSHD      |   32.4 | 67.6 |

## char:site

### Counts

| condition |  01 |  32 | All |
| :-------- | --: | --: | --: |
| Control   |   9 |   0 |   9 |
| FSHD      |  16 |  18 |  34 |
| All       |  25 |  18 |  43 |

### Row percentages (%)

| condition |   01 |   32 |
| :-------- | ---: | ---: |
| Control   |  100 |    0 |
| FSHD      | 47.1 | 52.9 |

# Statistical test results

| test           | statistic | p_value   | cramers_v | min_expected | warning                                                                  | variable          | analysis                         | effect_size | effect_measure | p_adjusted | significant |
| :------------- | :-------- | :-------- | --------: | :----------- | :----------------------------------------------------------------------- | :---------------- | :------------------------------- | :---------- | :------------- | :--------- | :---------- |
| Chi-square     | 3.449     | 0.3274    | 0.0948663 | 1.05         | Sparse expected frequencies: chi-square approximation may be unreliable. | char:batch        | Condition vs observed categories | 0.0949      | Cramér's V     | 0.3274     | False       |
| Not applicable | NA        | NA        |       nan | NA           | Insufficient categories                                                  | instrument_model  | Condition vs observed categories | NA          | Cramér's V     | NA         | False       |
| Not applicable | NA        | NA        |       nan | NA           | Insufficient categories                                                  | char:tissue       | Condition vs observed categories | NA          | Cramér's V     | NA         | False       |
| Not applicable | NA        | NA        |       nan | NA           | Insufficient categories                                                  | char:gender       | Condition vs observed categories | NA          | Cramér's V     | NA         | False       |
| Fisher exact   | inf       | 1.773e-09 |         1 | 1.88         |                                                                          | char:gender       | Condition vs missingness         | 1.0000      | Cramér's V     | 4.433e-09  | True        |
| Not applicable | NA        | NA        |       nan | NA           | Insufficient observations                                                | char:age          | Condition vs continuous variable | NA          | Not calculated | NA         | False       |
| Fisher exact   | inf       | 1.773e-09 |         1 | 1.88         |                                                                          | char:age          | Condition vs missingness         | 1.0000      | Cramér's V     | 4.433e-09  | True        |
| Fisher exact   | 4.182     | 0.1221    |   0.24313 | 3.56         |                                                                          | char:infer_gender | Condition vs observed categories | 0.2431      | Cramér's V     | 0.1526     | False       |
| Fisher exact   | inf       | 0.005649  |  0.413336 | 3.77         |                                                                          | char:site         | Condition vs observed categories | 0.4133      | Cramér's V     | 0.009415   | True        |

## Interpretation

Adjusted p-values below 0.05 indicate statistically significant associations after FDR correction.

Cramér's V measures the strength of categorical associations. Statistical significance should be interpreted alongside effect size and sample size.
