# FSHD Translator reference lists

Snapshot prepared 24 November 2026 for the AI for Omics transcriptomics practical.

**Where the files live.** The lists are data, so they stay out of git, like everything in `data/`; the course asks that the pull request contain no data. Extract the full snapshot (all six files and this manifest) into `practical/data/raw/bonus_lists/`, which is where `scripts/bonus_select_genes.py` reads them. The repository tracks only this README, `sources.tsv` (source URL, version, retrieval date, method, record count and SHA-256 of each file) and `SHA256SUMS`.

**Getting them again, exactly.** The five files read by `scripts/bonus_select_genes.py` were committed once, at commit `23edfb7`, before data left the repository; the Yao workbook comes from the EBI BioStudies mirror of its PMC supplement (the PMC address in `sources.tsv` may answer with a browser challenge instead of the file). From `practical/`:

```bash
mkdir -p data/raw/bonus_lists
for f in EpiGenes_main.csv QuickGO_GO0043403_human.tsv TF_names_v_1.01.txt \
         Tullai2007_NHLBI_ESBL_ieg.tsv Yao2014_TableS1_DUX4_induced_genes.tsv; do
  git show 23edfb7:practical/bonus_lists/$f > data/raw/bonus_lists/$f
done
curl -fL -o data/raw/bonus_lists/supp_ddu251_ddu251supp_table1.xls \
  https://www.ebi.ac.uk/biostudies/files/S-EPMC4168822/supp_ddu251_ddu251supp_table1.xls
(cd data/raw/bonus_lists && sha256sum -c ../../../bonus_lists/SHA256SUMS)
```

**Getting them from the sources.** Several source pages change over time, so a rebuild may not match the checksums.
- Download the three original files from the URLs in `sources.tsv`: the Lambert TF list, EpiFactors and the Yao workbook.
- Rebuild the three extracts as described in the `acquisition_and_scope` column: the Tullai ESBL table, the QuickGO human annotations, and the two-column Yao TSV.
- Then check them against the tracked manifest. From `practical/`: `(cd data/raw/bonus_lists && sha256sum -c ../../../bonus_lists/SHA256SUMS)`.

All five reference datasets are included. Yao Supplementary Table S1 is provided as the original workbook and an extracted two-column TSV.

## Contents

- `Tullai2007_NHLBI_ESBL_ieg.tsv` contains the 46 populated rows of the NHLBI ESBL table, preserving source symbols and UniProt accessions. It represents the ESBL list of immediate-early genes from growth-factor-stimulated T98G cells. Use human gene symbols and documented aliases for mapping, because the table's accessions are not a human-only identifier set.
- `TF_names_v_1.01.txt` contains the 1,639 symbols in Lambert and colleagues' posted TF-only list, reconstructed from the complete plain-text source with normalized line endings.
- `EpiGenes_main.csv` is the original EpiFactors v2.1 all-proteins CSV with 801 records. It includes the `Function` and `Modification` fields. Choose and record the eligible EpiFactors classes before ranking.
- `QuickGO_GO0043403_human.tsv` contains the 45 annotations shown across both source pages. The query requested human taxon 9606 and GO:0043403 with descendants reached through `is_a` and `part_of`. All shown qualifiers were `involved_in`. The TSV captures the source gene product ID, symbol, qualifier, GO ID, evidence, reference, taxon and agency. Symbols are intentionally preserved as shown, including aliases and an accession-like symbol. Deduplicate after resolving human GENCODE IDs.

## Yao DUX4-induced gene exclusion

- `supp_ddu251_ddu251supp_table1.xls` is the original Supplementary Table S1 workbook, kept byte-for-byte.
- `Yao2014_TableS1_DUX4_induced_genes.tsv` contains every one of the 213 gene rows from worksheet `Dux4.target.csv`, in the published order. It preserves the Ensembl gene ID and the original symbol. The workbook defines which genes belong to the set. The displayed fold changes are rounded.
- Match on the Ensembl gene ID, dropping any version suffix from the analysis IDs before matching. Keep a record of unmatched identifiers when applying the list to GENCODE v24. The 213 IDs have 212 distinct source symbols because `KHDC1L` labels two IDs (ENSG00000243501, ENSG00000256980).

## Provenance

`sources.tsv` records source URLs, versions, acquisition method, record counts and SHA-256 checksums. From the extracted directory, `sha256sum -c SHA256SUMS` checks the six packaged data files. Record the exact composition exclusions and the human identifier mapping separately before running the bonus ranking.

The Track B output branch was verified at [fc37b89](https://github.com/pierrebaumierresearch/AI4Transcriptomics_AI4OH_2026/commit/fc37b89a388acc31012f7e64851230ba8a869478), the Rochester replacement check and DUX4-score models commit.

EpiFactors v2.1 is attributed to Yulia A. Medvedeva, Daria Marakulina, Ilya E. Vorontsov, Ivan V. Kulakovskiy, Andreas Lennartsson and Finn Drabløs under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The source URL in the manifest links to the unmodified data. Cite the underlying papers and resources when using these lists.
