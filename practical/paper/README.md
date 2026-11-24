# The paper

Wang LH, Friedman SD, Shaw D, Snider L, Wong CJ, Budech CB, Poliachik SL, Gove NE,
Lewis LM, Campbell AE, Lemmers RJLF, van der Maarel SM, Tapscott SJ, Tawil RN.
**MRI-informed muscle biopsies correlate MRI with pathology and DUX4 target gene
expression in FSHD.** Human Molecular Genetics, 2019, 28(3):476-486.

| Where | Link |
|-------|------|
| DOI | https://doi.org/10.1093/hmg/ddy364 |
| PubMed | https://pubmed.ncbi.nlm.nih.gov/30312408/ |
| PubMed Central (open access full text and supplementary data) | https://pmc.ncbi.nlm.nih.gov/articles/PMC6337697/ |
| GEO series | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115650 |

Save the PDF and the supplementary archive in `doc/`, which git ignores. Share
the citation, not the files.

The supplementary archive holds the tables the main text refers to, among them
the description of each subject and biopsy and the list of DUX4 target genes.
The Methods section also links the authors' R code. If you read it, say so in
your report and say what you took from it.

## paper_claims.csv

One row per result of the paper that you test against your own analysis. Fill
it as you read the paper.

| column | content |
|--------|---------|
| claim_id | short id, for example `C1` |
| figure_or_table | where the paper shows it, for example `Fig 4`, `Table S5`, `Results, last section` |
| comparison | the groups compared, in the paper's words |
| my_contrast | the result file of yours that tests it, for example `G4_vs_Control` for `results/step3_de/de_G4_vs_Control.csv`. Leave empty to test the row against all your contrasts |
| gene | gene symbol or Ensembl id, one per row |
| direction_in_paper | `up` or `down` in the first group of the comparison |
| value_in_paper | the number the paper gives, if any |
| notes | threshold, software, sample subset |

In this fork, the DE results sit in two tracks, so `my_contrast` takes a qualified key:
`track_a/<contrast>` for `results/step3_de/track_a/de_<contrast>.csv.gz` (for example
`track_a/S10_G4`), or `track_b/groups/<level>_vs_Control` for
`results/step3_de/track_b/groups/de_<level>_vs_Control.csv.gz` (for example
`track_b/groups/G4_vs_Control`). The script reports any key that matches no result file.

A claim that does not concern one gene, such as a count of biopsies or a
correlation, goes in the tables of `worksheets/step4_comparison.md`.
