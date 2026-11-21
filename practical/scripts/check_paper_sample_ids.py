"""Check the paper's Table S1 against the GEO sample table: subject IDs, then per-subject values.

The paper (Wang et al. 2019, Supplementary Table S1) lists one row per FSHD subject.
This script
  1. compares the Table S1 subject IDs with the GEO sample identifiers, using both
     candidate keys (`char:sample_id` and `title`), in both directions;
  2. summarises the cohort as Table S1 describes it (genotype, sex, age, missing samples),
     for all subjects and for those with RNA-seq, so it can be compared with the paper's text;
  3. compares, subject by subject, the values that Table S1 and GEO both report.

Input:
  results/step1_audit/sample_table.csv   written by step 1
  Supplemental Table S1.docx             from the paper's supplementary archive, kept in doc/ (git-ignored)

Output:
  results/step1_audit/paper_table_s1_id_check.txt

Run:  uv run python scripts/check_paper_sample_ids.py "doc/<folder>/Supplemental Table S1.docx"
"""

from __future__ import annotations

import argparse
import re
import zipfile
import xml.etree.ElementTree as ET

import pandas as pd

import config

OUT = config.RESULTS / "step1_audit" / "paper_table_s1_id_check.txt"
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
SUBJECT_ID = re.compile(r"^\d{2}-\d{4}$")  # e.g. 01-0026, 32-0002

# Table S1 columns, in order. The four muscle columns hold the side (R/L) of the biopsied muscle.
S1_COLUMNS = ["id", "age", "sex", "css", "short4q", "Q", "G", "TA", "H",
              "path", "inflam", "active", "stir", "t1", "fat_fraction", "score", "rank"]
MUSCLES = ["Q", "G", "TA", "H"]

# Values reported by both sources: Table S1 column -> GEO column.
SHARED = {"age": "char:age", "sex": "char:gender", "path": "char:path.score",
          "inflam": "char:inflam", "active": "char:active", "stir": "char:stir",
          "t1": "char:t1", "fat_fraction": "char:fat.fraction",
          "score": "char:rnaseq.score", "rank": "char:rnaseq.rank"}


def cell_text(cell) -> str:
    return "".join(t.text or "" for t in cell.iter(f"{W}t")).strip()


def read_table_s1(docx_path) -> tuple[pd.DataFrame, list[str]]:
    """Subject rows of Table S1 (header rows skipped) and the notes printed under the title."""
    with zipfile.ZipFile(docx_path) as z:
        root = ET.fromstring(z.read("word/document.xml"))
    rows = []
    for row in root.iter(f"{W}tr"):
        cells = [cell_text(c) for c in row.iter(f"{W}tc")]
        if cells and SUBJECT_ID.match(cells[0]):
            if len(cells) != len(S1_COLUMNS):
                raise ValueError(f"row {cells[0]}: {len(cells)} cells, expected {len(S1_COLUMNS)}")
            rows.append(cells)
    body = root.find(f"{W}body")
    notes = [cell_text(p) for p in body.findall(f"{W}p") if cell_text(p)]  # paragraphs outside the table
    return pd.DataFrame(rows, columns=S1_COLUMNS), notes


def genotype(short4q: str) -> str:
    """Table S1 marks FSHD2 with '*' and mosaic with '**' after the 4q allele size."""
    if short4q.endswith("**"):
        return "mosaic"
    if short4q.endswith("*"):
        return "FSHD2"
    return "FSHD1 (unmarked)"


def cohort_summary(df: pd.DataFrame, label: str) -> list[str]:
    age = pd.to_numeric(df["age"])
    n_female = int((df["sex"] == "F").sum())
    return [
        f"{label}: n = {len(df)}",
        f"  genotype (from Short4q marks): {df['genotype'].value_counts().to_dict()}",
        f"  sex: {df['sex'].value_counts().to_dict()}; female = {n_female}/{len(df)} = {n_female / len(df):.1%}",
        f"  age: mean {age.mean():.1f}, range {age.min()}-{age.max()}",
        f"  muscle biopsied: {df['muscle'].value_counts().to_dict()}",
    ]


def compare_values(merged: pd.DataFrame) -> list[str]:
    """Subject-by-subject comparison of the values both sources report (numeric when possible)."""
    lines = []
    for s1_col, geo_col in SHARED.items():
        s1 = merged[s1_col].fillna("").astype(str).str.strip()
        geo = merged[geo_col].fillna("").astype(str).str.strip()
        if s1_col == "sex":
            geo = geo.str[:1].str.upper()  # GEO writes male/female, Table S1 M/F
            differs = s1 != geo
        else:
            s1_num = pd.to_numeric(s1.replace("", None))
            geo_num = pd.to_numeric(geo.replace("", None))
            both_missing = s1_num.isna() & geo_num.isna()
            differs = ~(both_missing | ((s1_num - geo_num).abs() < 1e-6))
        bad = merged.loc[differs, "char:sample_id"].tolist()
        detail = "; ".join(f"{i}: S1={a!r} GEO={b!r}" for i, a, b in
                           zip(bad, s1[differs], geo[differs])) or "none"
        lines.append(f"  {s1_col:<13} vs {geo_col:<20} mismatches: {len(bad)}  {detail}")
    group = merged["source_name_ch1"].str.extract(r"Group(\d)")[0]
    bad = merged.loc[group != merged["rank"], "char:sample_id"].tolist()
    lines.append(f"  {'rank':<13} vs {'Group in source_name':<20} mismatches: {len(bad)}  {', '.join(bad) or 'none'}")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("table_s1", help="path to 'Supplemental Table S1.docx'")
    args = parser.parse_args()

    s1, notes = read_table_s1(args.table_s1)
    s1["genotype"] = s1["short4q"].map(genotype)
    s1["muscle"] = s1[MUSCLES].apply(lambda r: next((m for m in MUSCLES if r[m]), "none"), axis=1)
    s1["has_rnaseq"] = s1["score"] != ""

    samples = pd.read_csv(config.RESULTS / "step1_audit" / "sample_table.csv", dtype=str)
    fshd = samples[samples["condition"] == "FSHD"]
    controls = samples[samples["condition"] == "Control"]
    paper = set(s1["id"])

    lines = ["== 1. Subject IDs ==",
             f"Table S1 subject IDs: {len(s1)} (unique: {len(paper)})",
             f"GEO FSHD samples: {len(fshd)}; GEO control samples: {len(controls)}"]
    for key in ["char:sample_id", "title"]:
        geo = set(fshd[key])
        lines += [
            "",
            f"[{key}] GEO FSHD unique values: {len(geo)}",
            f"[{key}] in both: {len(paper & geo)}",
            f"[{key}] only in Table S1: {', '.join(sorted(paper - geo)) or 'none'}",
            f"[{key}] only in GEO FSHD: {', '.join(sorted(geo - paper)) or 'none'}",
        ]
    found = sorted(paper & set(controls["char:sample_id"]))
    lines += ["", f"GEO controls found in Table S1: {', '.join(found) or 'none'}"]

    lines += ["", "== 2. Cohort as described by Table S1 ==", "Table S1 notes:"]
    lines += [f"  {n}" for n in notes]
    lines += [""] + cohort_summary(s1, "All subjects")
    lines += [""] + cohort_summary(s1[s1["has_rnaseq"]], "Subjects with an RNA-seq score")
    lines += ["",
              f"No RNA-seq score: {', '.join(s1.loc[~s1['has_rnaseq'], 'id'])}",
              f"No pathology score: {', '.join(s1.loc[s1['path'] == '', 'id'])}"]

    merged = fshd.merge(s1, left_on="char:sample_id", right_on="id", how="inner")
    lines += ["", f"== 3. Values reported by both (joined on char:sample_id, {len(merged)} subjects) =="]
    lines += compare_values(merged)

    OUT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwritten: {OUT.relative_to(config.ROOT)}")


if __name__ == "__main__":
    main()
