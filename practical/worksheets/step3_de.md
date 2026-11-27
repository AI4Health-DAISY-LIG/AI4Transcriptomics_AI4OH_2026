# Step 3: differential expression with PyDESeq2

_Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence._

## Design

Input matrix (which file, which filter):
`counts_filtered.csv` (from Step 2) with outliers `32-0016` and `01-0037` removed (41 samples: 9 Controls, 32 FSHD across G1–G4; 19,258 genes).

Design formula:
Primary model (condition): `~ batch + sex + site + hemoglobin_pct + neutrophil_pct + condition`
Primary model (subgroups): `~ batch + sex + site + hemoglobin_pct + neutrophil_pct + group`
Sensitivity models: `no_batch`, `no_site`, `no_blood` (omitting `hemoglobin_pct` and `neutrophil_pct`), `site01`, `site01_no_blood` (benchmarked against `site01`), `batch_III_IV`, `without_32-0014`, and `with_outliers`.

Covariates included and why, or why not:

- **`sex` (`infer_gender`):** Top driver of PC2 in Step 2 ($\eta^2 = 0.51, p = 1.35 \times 10^{-7}$), orthogonal to disease status (`PC1` $\eta^2 = 0.058$), cleanly removing sex-linked expression variance.
- **`batch`:** Major source of technical variance (`PC1` $\eta^2 = 0.39$, `PC2` $\eta^2 = 0.16$), included to account for batch/library-depth effects across the cohort.
- **`site`:** Present across all 41 samples (`PC1` $\eta^2 = 0.21, p = 0.0027$); included to adjust for collection-site variation. Sensitivity runs (`no_site`, `site01`) confirmed that `site` is estimable without distorting effect directions (`0` sign flips) and improves statistical power in mild subgroups (17 vs. 8 DEGs in G1).
- **`hemoglobin_pct` & `neutrophil_pct`:** Sample-level compositional proxies available for all 41 samples, included to control for peripheral blood and acute neutrophil contamination introduced during biopsy extraction without regressing out muscle-intrinsic pathology.

Covariates you wanted and could not include, and why:

- **`age`:** Missing for all 9 healthy control samples, making it non-estimable in the full cohort (and within FSHD cases alone, `age` explained virtually zero global variance: `PC1` $R^2 = 0.0035, p = 0.75$; `PC2` $R^2 = 0.0077, p = 0.63$).
- **Disease severity & imaging metrics (`fat.fraction`, `inflam`, `path.score`, `rnaseq.score`, `rnaseq.rank`, `stir`, `t1`, `active`):** Missing for healthy controls (`n = 31–32`) and strongly collinear with `condition` and subgroup severity (`PC1` $R^2 / \eta^2 = 0.21\text{–}0.66$), meaning inclusion would either drop all controls or regress out the primary biological signal of interest.

## Results

run: primary samples: 41 genes: 19258\
design: `~ batch + sex + site + hemoglobin_pct + neutrophil_pct + group`

Number of samples per group:\
Control=9, FSHD=32 (Total 41).\
Breakdown by Group for secondary analysis: Control=9, G1=8, G2=5, G3=13, G4=6.

| Comparison         | Genes tested | padj < 0.05 | of which up | of which down |
| ------------------ | ------------ | ----------- | ----------- | ------------- |
| FSHD vs Control    | 19258        | 901         | 811         | 90            |
| Group 1 vs Control | 19258        | 17          | 16          | 1             |
| Group 2 vs Control | 19258        | 497         | 444         | 53            |
| Group 3 vs Control | 19258        | 2233        | 1789        | 444           |
| Group 4 vs Control | 19258        | 6392        | 3603        | 2789          |

Figures: MA plot, volcano plot, p-value histogram (in `results/step3_de/`).

Top 10 genes by adjusted p-value for FSHD vs Control comparison, with log2 fold change:

| Rank | Gene Symbol | ENSG ID            | Log2FC | Padj     | Note                                                                                   |
| ---- | ----------- | ------------------ | ------ | -------- | -------------------------------------------------------------------------------------- |
| 1    | COMP        | ENSG00000105664.10 | 6.80   | 1.64e-07 | Cartilage Oligomeric Matrix Protein; ECM remodeling and fibrosis marker.               |
| 2    | SPP1        | ENSG00000118785.13 | 6.44   | 4.44e-06 | Osteopontin (Secreted Phosphoprotein 1); macrophage/inflammation marker.               |
| 3    | TNMD        | ENSG00000000005.5  | 5.53   | 4.44e-06 | Tenomodulin; tendon/ECM structural remodeling glycoprotein.                            |
| 4    | EGR1        | ENSG00000120738.7  | 3.47   | 1.66e-05 | Early Growth Response 1; immediate-early stress/injury transcription factor.           |
| 5    | TNC         | ENSG00000041982.15 | 3.50   | 1.66e-05 | Tenascin C; ECM glycoprotein upregulated during tissue injury/remodeling.              |
| 6    | CTGF        | ENSG00000118523.5  | 2.30   | 1.95e-05 | Connective Tissue Growth Factor (CCN2); profibrotic signaling mediator.                |
| 7    | EGR2        | ENSG00000122877.13 | 3.16   | 1.95e-05 | Early Growth Response 2; stress response and peripheral nerve/Schwann cell regulation. |
| 8    | PRG4        | ENSG00000116690.12 | 4.27   | 3.64e-05 | Proteoglycan 4 (Lubricin); ECM/joint marker induced by inflammation.                   |
| 9    | CDKN1A      | ENSG00000124762.13 | 3.21   | 6.59e-05 | p21/CIP1; cyclin-dependent kinase inhibitor, cell cycle arrest/stress.                 |
| 10   | AREG        | ENSG00000109321.10 | 6.27   | 8.57e-05 | Amphiregulin; EGFR ligand involved in tissue repair and inflammation.                  |

_Observation:_ The top hits are dominated by **ECM/fibrotic remodeling** (COMP, TNMD, TNC, CTGF, PRG4), **inflammatory/macrophage** (SPP1, AREG), and **cellular stress/cell-cycle arrest** (EGR1, EGR2, CDKN1A) markers rather than canonical DUX4 targets (e.g., KHDC1L). This suggests that while DUX4 activation is restricted to focal bursts or specific advanced subgroups (G3/G4), the global FSHD vs Control transcriptional signature is heavily driven by secondary fibrotic remodeling and chronic stress/repair responses.

## Check

### Describe the p-value histogram

Each p-value histogram shows a sharp peak near 0 (reflecting true differential expression signal) followed by a flat, uniform distribution toward 1. The height of the peak near 0 scales directly with disease progression/group severity : smallest in Group 1 vs Control and largest in Group 4 vs Control.

### Run the analysis with and without batch in the design

| Contrast            | Sig. (With Batch) | Sig. (Without Batch) | Up (With / Without) | Down (With / Without) | Change from Adding Batch |
| :------------------ | :---------------- | :------------------- | :------------------ | :-------------------- | :----------------------- |
| **FSHD vs Control** | 901               | 358                  | 811 / 324           | 90 / 34               | +543 (+151.7%)           |
| **G1 vs Control**   | 17                | 5                    | 16 / 4              | 1 / 1                 | +12 (+240.0%)            |
| **G2 vs Control**   | 497               | 553                  | 444 / 492           | 53 / 61               | -56 (-10.1%)             |
| **G3 vs Control**   | 2,233             | 2,703                | 1,789 / 2,174       | 444 / 529             | -470 (-17.4%)            |
| **G4 vs Control**   | 6,392             | 8,185                | 3,603 / 4,408       | 2,789 / 3,777         | -1,793 (-21.9%)          |

Including `batch` in the model shifts the number of significant genes in opposite directions depending on the comparison:

- **Global FSHD vs Control & G1:** Adding `batch` _increases_ statistical power, raising significant DEGs from 358 to 901 in FSHD vs Control (and from 5 to 17 in G1 vs Control) by absorbing unexplained technical variance across sequencing batches. Notably, top hits such as _COMP_ ($padj$: $2.78\times 10^{-6} \to 1.64\times 10^{-7}$), _EGR1_ ($6.86\times 10^{-4} \to 1.66\times 10^{-5}$), and _PRG4_ ($1.08\times 10^{-3} \to 3.64\times 10^{-5}$) become substantially stronger, and _SPP1_ emerges as the #2 hit ($padj = 4.44\times 10^{-6}$).
- **Subgroups G2–G4:** Adding `batch` _reduces_ significant DEG counts by ~10–22% (e.g., 8,185 down to 6,392 in G4 vs Control), removing false positives caused by partial confounding between sequencing batches and more severe disease subgroups.

### Group 1 vs control analysis

Under the primary subgroup design (`~ batch + sex + site + hemoglobin_pct + neutrophil_pct + group`, comparing `Control=9` vs `G1=8` across 19,258 tested genes):

- **Differentially expressed genes (`padj < 0.05`):** **17 genes** total (**16 upregulated**, **1 downregulated**), compared to 901 DEGs in the overall `FSHD_vs_Control` model (`~ batch + sex + site + hemoglobin_pct + neutrophil_pct + condition`, `Control=9` vs `FSHD=32`):
  - **16 upregulated:** _MTND6P4_ ($\log_2\text{FC} = 4.67, padj = 0.0068$), _COMP_ ($5.22, 0.0068$), _TNC_ ($2.69, 0.0107$), _PRG4_ ($3.97, 0.0107$), _EGR2_ ($2.63, 0.0107$), _EGR1_ ($2.83, 0.0107$), _KERA_ ($5.57, 0.0142$), _APLN_ ($1.73, 0.0152$), _TNMD_ ($4.43, 0.0152$), _PHYKPL_ ($0.54, 0.0172$), _SERP2_ ($1.48, 0.0172$), _NR4A2_ ($2.37, 0.0223$), _CTGF_ ($1.83, 0.0396$), _TCEAL7_ ($1.89, 0.0412$), _DNER_ ($3.23, 0.0440$), and _EPB41L3_ ($1.32, 0.0496$).
  - **1 downregulated:** _MTRNR2L1_ ($\log_2\text{FC} = -6.95, padj = 0.0068$).

- **Robustness across sensitivity analyses:**
  - Across all sensitivity runs (`no_batch`, `no_site`, `no_blood`, `site01`, `site01_no_blood` vs. `site01`, `batch_III_IV`, `without_32-0014`, `with_outliers`), **zero shared DEGs ever reverse direction** (`shared_opposite_direction = 0` across all 40 comparisons), and **100% of the top 100 reference-ranked genes maintain identical fold-change direction** (`top_same_direction = 100/100`).
  - Both genome-wide rank correlations (`lfc_spearman` $= 0.909\text{–}0.963$ across full-cohort models; $0.972$ for `site01_no_blood` vs. `site01`) and effect-size magnitudes (`median_effect_ratio` $= 0.865\text{–}1.055$) remain highly stable in `G1_vs_Control`.
  - Removing `32-0014` retains **17/17** primary DEGs while expanding G1 DEGs to 45 (`median_effect_ratio = 1.049`), confirming the signature is not driven by an outlier or covariate artifact.

- **Biological & Technical Interpretation:**
  - **Earliest subclinical FSHD remodeling (primary interpretation):** Rather than active _DUX4_ target bursts or overt inflammation (which emerge in G2–G4), G1 muscle exhibits a reproducible early signature of extracellular matrix remodeling (_COMP_, _TNC_, _CTGF_, _PHYKPL_), myogenic/angiogenic adaptation (_TCEAL7_, _APLN_), and mitochondrial stress (_MTRNR2L1_). This is the same results that were explained by the authors and are from their supplementary 4 table.
  - **Non-disease factors that may still modulate effect sizes:**
    1. **Anatomical micro-sampling & connective tissue proportion:** Even after adjusting for anatomical `site`, subtle needle-biopsy variation in depth or proximity to fascia/myotendinous junctions can enrich for resident tenogenic/cartilage ECM transcripts (_KERA_, _TNMD_, _COMP_, _PRG4_).
    2. **Ex vivo ischemia & biopsy excision stress:** Mechanical excision and brief delays before flash-freezing rapidly induce immediate-early transcription factors (_EGR1_, _EGR2_, _NR4A2_) within minutes of sampling.
    3. **Pre-biopsy physical activity / mechanical loading:** Routine exertion prior to clinic visits transiently elevates angiogenic (_APLN_), immediate-early (_EGR1_, _NR4A2_), and matrix-turnover (_TNC_, _CTGF_) transcripts.
    4. **Unmodeled clinical/demographic covariates:** Differences in age (not included in the model), habitual activity levels, or NSAID/medication use between G1 patients and healthy volunteers.

## Decisions I made

- **Adopted the full covariate-adjusted models (`~ batch + sex + site + hemoglobin_pct + neutrophil_pct + condition` and `+ group`) as the primary design:** Adjusting for `batch`, `site`, and sample-level blood/immune composition proxies (`hemoglobin_pct`, `neutrophil_pct`) alongside `sex` suppressed technical and batch-driven inflation in severe subgroups (reducing G4 DEGs from 8,185 in `no_batch` to 6,392) while increasing statistical power to detect consistent disease signal in the global `FSHD_vs_Control` comparison (increasing significant DEGs from 358 to 901).
- **Retained outlier exclusions (`32-0016` and `01-0037`) from Step 2:** Kept the filtered 41-sample dataset (`9` Controls, `32` FSHD; `19,258` expressed genes) to prevent extreme sample outliers from inflating gene-wise dispersion estimates in PyDESeq2, while verifying via the `with_outliers` run that excluding them does not distort fold-change directions (`lfc_spearman` $= 0.942\text{–}0.996$, `0` sign flips).
- **Evaluated both global (`condition`) and severity-stratified (`group`: G1–G4) contrasts:** Modeling subgroups separately captured the monotonic escalation of transcriptomic remodeling across disease stages (**G1: 17** $\rightarrow$ **G2: 497** $\rightarrow$ **G3: 2,233** $\rightarrow$ **G4: 6,392** DEGs) that is otherwise masked in a pooled FSHD comparison.
- **Anchored conclusions on cross-model concordance rather than a single specification:** Benchmarking across 9 sensitivity models (`no_batch`, `no_site`, `no_blood`, `site01`, `site01_no_blood`, `batch_III_IV`, `without_32-0014`, `with_outliers`) demonstrated **zero directional reversals** among shared DEGs (`shared_opposite_direction = 0`), **100% sign concordance** across the top 100 reference genes, and high genome-wide rank correlations (`lfc_spearman` $= 0.909\text{–}0.999$ across full-cohort models), confirming that the core biological findings are robust to model specification.
