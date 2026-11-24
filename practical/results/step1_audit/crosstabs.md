# Condition against sample variables

## char:batch

| condition   |   I |   II |   III |   IV |   All |
|:------------|----:|-----:|------:|-----:|------:|
| Control     |   0 |    0 |     6 |    3 |     9 |
| FSHD        |   5 |    5 |    16 |    8 |    34 |
| All         |   5 |    5 |    22 |   11 |    43 |

## instrument_model

| condition   |   Illumina HiSeq 2500 |   All |
|:------------|----------------------:|------:|
| Control     |                     9 |     9 |
| FSHD        |                    34 |    34 |
| All         |                    43 |    43 |

## char:tissue

| condition   |   muscle biopsy |   All |
|:------------|----------------:|------:|
| Control     |               9 |     9 |
| FSHD        |              34 |    34 |
| All         |              43 |    43 |

## char:gender

| condition   |   female |   male |   missing |   All |
|:------------|---------:|-------:|----------:|------:|
| Control     |        0 |      0 |         9 |     9 |
| FSHD        |       12 |     22 |         0 |    34 |
| All         |       12 |     22 |         9 |    43 |

## char:age

| condition   |   [0, 30) |   [30, 45) |   [45, 60) |   [60, 120) |   missing |   All |
|:------------|----------:|-----------:|-----------:|------------:|----------:|------:|
| Control     |         0 |          0 |          0 |           0 |         9 |     9 |
| FSHD        |         2 |          5 |         14 |          13 |         0 |    34 |
| All         |         2 |          5 |         14 |          13 |         9 |    43 |
