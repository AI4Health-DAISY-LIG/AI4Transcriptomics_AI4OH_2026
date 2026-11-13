"""Fetch a GEO series by code.

Every function here makes one HTTPS request per file and nothing else:
no crawling, no retries in loops, no third-party mirrors. Files that are
already on disk are not downloaded again.

Two GEO endpoints are used:

* the series record in SOFT text format
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSExxxx&targ=all&form=text&view=brief
  which holds the series description, the list of supplementary files and
  one block per sample (GSM) with its characteristics;
* the supplementary files on the GEO FTP server, served over HTTPS at
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSExxxnnn/GSExxxx/suppl/<file>.
"""

from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path

import pandas as pd

GEO_TEXT = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc={acc}&targ={targ}&form=text&view=brief"
)
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}nnn/{acc}/suppl/{name}"
USER_AGENT = "ai4transcriptomics-practical (python urllib; course material, one request per file)"


def suppl_url(acc: str, name: str) -> str:
    """HTTPS URL of a series-level supplementary file, e.g. GSE115650_biopsy_genecounts.csv.gz."""
    return GEO_FTP.format(stem=acc[:-3], acc=acc, name=name)


def download(url: str, dest: Path, overwrite: bool = False) -> Path:
    """Download one URL to dest. Skips the request when dest already exists."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not overwrite:
        print(f"[skip] {dest.name} already on disk")
        return dest
    print(f"[get ] {url}")
    part = dest.with_name(dest.name + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=300) as resp, open(part, "wb") as out:
        shutil.copyfileobj(resp, out)
    part.replace(dest)
    print(f"[done] {dest.name}: {dest.stat().st_size / 1e6:.1f} MB")
    return dest


def fetch_soft(acc: str, destdir: Path, targ: str = "all") -> Path:
    """Series and sample records in SOFT text format, in one request."""
    return download(GEO_TEXT.format(acc=acc, targ=targ), Path(destdir) / f"{acc}_{targ}.soft.txt")


def parse_soft(path: Path) -> tuple[dict[str, list[str]], pd.DataFrame]:
    """Parse a SOFT text file into the series record and a sample table.

    Returns (series, samples). `series` maps each `!Series_*` key to the list of
    its values. `samples` has one row per GSM; `!Sample_characteristics_ch1`
    lines of the form "key: value" become columns named `char:key`, every other
    `!Sample_*` key becomes a column (multi-valued keys are joined with " | ").
    """
    series: dict[str, list[str]] = {}
    samples: list[dict] = []
    current: dict | None = None
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("^SERIES"):
            current = series
            continue
        if line.startswith("^SAMPLE"):
            current = {"gsm": line.partition("=")[2].strip()}
            samples.append(current)
            continue
        if line.startswith("^"):
            current = None  # platform block and any other entity: not parsed
            continue
        if current is None or not line.startswith("!"):
            continue
        key, _, value = line[1:].partition(" = ")
        key = key.removeprefix("Sample_").removeprefix("Series_")
        if key == "characteristics_ch1":
            name, _, val = value.partition(": ")
            current[f"char:{name.strip()}"] = val.strip()
        else:
            current.setdefault(key, []).append(value)
    table = pd.DataFrame(samples)
    for col in table.columns:
        if table[col].map(lambda x: isinstance(x, list)).any():
            table[col] = table[col].map(
                lambda x: " | ".join(x) if isinstance(x, list) else x
            )
    return series, table


def series_supplementary_files(series: dict[str, list[str]]) -> list[str]:
    """File names listed under !Series_supplementary_file."""
    return [url.rsplit("/", 1)[-1] for url in series.get("supplementary_file", [])]
