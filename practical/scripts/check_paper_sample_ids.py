"""Check that the subject IDs of the paper's Table S1 match the GEO sample identifiers.

The paper (Wang et al. 2019, Supplementary Table S1) lists one row per FSHD subject.
This script compares those IDs with the GEO sample table written by step 1, using
both candidate keys (`char:sample_id` and `title`), in both directions.

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


def table_s1_ids(docx_path) -> list[str]:
    """First-column values of the Word table that look like subject IDs (header rows skipped)."""
    with zipfile.ZipFile(docx_path) as z:
        root = ET.fromstring(z.read("word/document.xml"))
    ids = []
    for row in root.iter(f"{W}tr"):
        first_cell = next(row.iter(f"{W}tc"), None)
        if first_cell is None:
            continue
        text = "".join(t.text or "" for t in first_cell.iter(f"{W}t")).strip()
        if SUBJECT_ID.match(text):
            ids.append(text)
    return ids


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("table_s1", help="path to 'Supplemental Table S1.docx'")
    args = parser.parse_args()

    paper_ids = table_s1_ids(args.table_s1)
    samples = pd.read_csv(config.RESULTS / "step1_audit" / "sample_table.csv", dtype=str)
    fshd = samples[samples["condition"] == "FSHD"]
    controls = samples[samples["condition"] == "Control"]
    paper = set(paper_ids)

    lines = [
        f"Table S1 subject IDs: {len(paper_ids)} (unique: {len(paper)})",
        f"GEO FSHD samples: {len(fshd)}; GEO control samples: {len(controls)}",
    ]
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

    OUT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwritten: {OUT.relative_to(config.ROOT)}")


if __name__ == "__main__":
    main()
