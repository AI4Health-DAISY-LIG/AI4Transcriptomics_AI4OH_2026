# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|------|------|---------------|------------------------|
|GSE115650_all.soft.txt|13M|desciption of the study and all of its 43 samples | I need it to differentiate FSHD from control |
|GSE115650_biopsy_genecounts.csv.gz|1.7M| table of counts for each gene (rows) in each biopsy (columns)|input for all the analysis |
|gencode.v24.annotation.gtf.gz|37M| table of each gene and its corresponding ENSG code| I need it to know which gene they are talking about in the article |

## Sample table

Number of samples: 43
Condition groups and counts: 34 for muscular dystrophy (FSHD) and 9 for control, 
Variables available per sample (list the columns of `sample_table.csv`): subject status, age, gender, path.score, inflam, active, stir, t1, fat.fraction, tissue, rnaseq.score, rnaseq.rank, sample_id, batch.
Variables that have a value for one condition only (see `missing_values.csv`): age, gender, path.score, inflam, active, stir, t1, fat.fraction, rnaseq.score and rnaseq.rank. 
I've also noticed that for the variables path.score, inflam and active, there are 32 values for FSHD but we have 34 samples of FSHD.

## Do the counts match the sample table?

What `counts_header_check.txt` shows, and what you decided: This file shows that for all 43 samples, GEO and "sample_id" match. For 42 samples, GEO and "title" match. There's one sample, the one on line 13 in 'sample_table.csv', for which the GEO and title don't match. It's either 01-0022-1 for title or 01-0022 for sample_id and GEO.

Here's what 'counts_header_check.txt' says : first columns: 01-0021, 01-0022-1, 01-0023, 01-0024, 01-0025, 01-0026, 01-0027, 01-0029

This means that there is not a column named "01-0022" in the raw counts file, so this is supposedly the same sample. 
This sample doesn't seem to be very different from the others when looking at variables.
In the article, this sample is never mentioned like being different from the others. 

I think there is no point in removing a sample without any biological explanation, also I will not change directly the title from 01-0022-1 to 01-002 since it will be based on the hypothesis that the samples are the same, which I'm not 100% sure. 

So, I decided that I will use "title" as a key instead of "sample_id" as it is done in script_2. Otherwise, the program might not recognise this sample. Furthermore, in raw_data, in the gene count file, samples are already named after "title" and not "sample_ID" in the columns. 


## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

### char:batch

Technical variable, it is not balanced accross conditions, the control samples were only run with batch III and IV, it is a confounder.

| condition   |   I |   II |   III |   IV |   All |
|:------------|----:|-----:|------:|-----:|------:|
| Control     |   0 |    0 |     6 |    3 |     9 |
| FSHD        |   5 |    5 |    16 |    8 |    34 |
| All         |   5 |    5 |    22 |   11 |    43 |

### instrument_model

Technical variable, it is balanced accross conditions because there is only one instrument model.

| condition   |   Illumina HiSeq 2500 |   All |
|:------------|----------------------:|------:|
| Control     |                     9 |     9 |
| FSHD        |                    34 |    34 |
| All         |                    43 |    43 |

### char:tissue

Biological variable, it is balanced accross conditions because there is only one tissue type.

| condition   |   muscle biopsy |   All |
|:------------|----------------:|------:|
| Control     |               9 |     9 |
| FSHD        |              34 |    34 |
| All         |              43 |    43 |

### char:gender

Biological variable, since I have no info for the controls, I can't tell if it is balanced.

| condition   |   female |   male |   missing |   All |
|:------------|---------:|-------:|----------:|------:|
| Control     |        0 |      0 |         9 |     9 |
| FSHD        |       12 |     22 |         0 |    34 |
| All         |       12 |     22 |         9 |    43 |

### char:age

Biological variable, since I have no info for the controls, I can't tell if it is balanced.

| condition   |   [0, 30) |   [30, 45) |   [45, 60) |   [60, 120) |   missing |   All |
|:------------|----------:|-----------:|-----------:|------------:|----------:|------:|
| Control     |         0 |          0 |          0 |           0 |         9 |     9 |
| FSHD        |         2 |          5 |         14 |          13 |         0 |    34 |
| All         |         2 |          5 |         14 |          13 |         9 |    43 |



## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

- The paper describes from which muscles the samples were collected and on GEO you don't have the information of which sample came from what muscle : "Each subject underwent a lower extremity muscle biopsy in a muscle safely accessible by needle biopsy and included tibialis anterior, gastrocnemius muscles, vastus lateralis and one hamstring (biceps femoris) muscle. T2-STIR hyperintense muscles were preferentially biopsied when possible."

-All procedures, including muscle biopsies, handling and processing of samples and muscle MRI protocols, were conducted accross two sites, but we don't know which one comes from where on GEO. My Hypothesis is that for geo-accession variable, some start with -32 and others start with -1, that might represent two sites.

-In the paper, it is said that "Recruitment for the study included 36 subjects (34% female) with genetically confirmed FSHD (34 FSHD1, 2 FSHD2)", but on GEO we don't know which sample is FSHD1 and which is FSHD2.

-"In three subjects, the histopathologic samples were inadequate for histopathologic scoring, and in two subjects, a muscle sample for RNA extraction and sequencing was not successfully obtained.", on GEO you don't have info about the two subjects without rna sequencing and only two subjects don't have histopatgologic scoring, what about the third one? 



## FAIR score

Total: /30 (see `fair_score.md`). The two weakest criteria and why:

## Decisions I made

-In script_2, changing the key that is supposed to link the gene counts in raw data and the sample_table _from "sample_ID" to "Title"