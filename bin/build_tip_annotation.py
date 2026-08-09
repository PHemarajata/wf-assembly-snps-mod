#!/usr/bin/env python3
"""
Map isolate metadata onto the tips of a grafted phylogeny, and emit the input
files that phylogeographic and molecular-clock tools expect.

Why this exists
---------------
The tip labels in the grafted tree are not the metadata `sample_id`. They carry
assembly-version suffixes, assembler suffixes, FASTA extensions and free-text
geography appended by earlier steps (`GCF_000756925_1_Australia_Townsville_
Queensland`, `SRR33188703_SPAdes`, `IP-0009-1-R.fas`). `GCF_` and `GCA_` are used
interchangeably between files for the same assembly. A naive join drops tips
silently, and a tip dropped from a phylogeographic analysis is not an error you
notice -- it is a country that quietly loses a genome.

This script normalises both sides to a single join key, requires a 100% match by
default, and writes the mapping and its QC alongside the annotation so the join
is auditable rather than assumed.

Geography
---------
`--country-col` selects which column carries country of ACQUISITION (default
`Country_Final`). `--diagnosis-col` carries country of DIAGNOSIS. Travel-
associated isolates -- an infection acquired in country A and diagnosed in
country B -- must enter a phylogeographic model under A; attributing them to B
invents transmission that never happened. The two columns are kept side by side
and every disagreement is written to `country_change_log.tsv`.

Values that name more than one country ("Panama and Peru") or a continent
("Africa") are NOT valid states in a discrete-trait model. They are flagged and
excluded from `geo_state` rather than silently coerced.

Dates
-----
`final_collection_dates` is mixed precision (YYYY, YYYY-MM, YYYY-MM-DD,
"unknown"). A year-only date is not 1 January: it is an interval. The script
emits a midpoint plus explicit bounds so TreeTime and LSD2 can integrate over
the uncertainty instead of pretending to a precision the data do not have.

Usage
-----
    build_tip_annotation.py \
        --tree        results/graft/global_grafted.treefile \
        --metadata    metadata.tsv \
        --clusters    Summaries/clusters.tsv \
        --cluster-qc  Summaries/cluster_phylogeny_summary.csv \
        --outdir      Tip_annotation

`--tips-from-csv FILE:COLUMN` substitutes for `--tree` when the Newick file is
not to hand (e.g. re-annotating a previous run from its `annotated_tips.csv`).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------
# Country -> region. Ten regions, matching the scheme used in the first
# interpretation pass so figures and tables stay comparable across runs.
# --------------------------------------------------------------------------
REGION = {
    "Australia": "Australasia", "New Zealand": "Australasia",
    "Papua New Guinea": "Oceania", "Micronesia": "Oceania",
    "Thailand": "Mainland SE Asia", "Vietnam": "Mainland SE Asia",
    "Laos": "Mainland SE Asia", "Myanmar": "Mainland SE Asia",
    "Cambodia": "Mainland SE Asia",
    "Malaysia": "Maritime SE Asia", "Singapore": "Maritime SE Asia",
    "Indonesia": "Maritime SE Asia", "Philippines": "Maritime SE Asia",
    "Brunei": "Maritime SE Asia",
    "China": "East Asia", "Hong Kong": "East Asia", "Taiwan": "East Asia",
    "Japan": "East Asia", "South Korea": "East Asia",
    "India": "South Asia", "Bangladesh": "South Asia", "Pakistan": "South Asia",
    "Sri Lanka": "South Asia", "Nepal": "South Asia",
    "USA": "Americas", "Mexico": "Americas", "Guatemala": "Americas",
    "Costa Rica": "Americas", "Panama": "Americas", "Ecuador": "Americas",
    "Peru": "Americas", "Brazil": "Americas", "Puerto Rico": "Americas",
    "Aruba": "Americas", "Martinique": "Americas", "Guadeloupe": "Americas",
    "Trinidad and Tobago": "Americas", "Virgin Islands": "Americas",
    "Panama and Peru": "Americas",
    "Ghana": "Africa", "South Africa": "Africa", "Nigeria": "Africa",
    "Kenya": "Africa", "Madagascar": "Africa", "Africa": "Africa",
    "France": "Europe", "United Kingdom": "Europe", "Portugal": "Europe",
    "Czech Republic": "Europe", "Russia": "Europe", "Switzerland": "Europe",
    "Netherlands": "Europe", "Germany": "Europe", "Spain": "Europe",
    "Italy": "Europe", "Sweden": "Europe", "Denmark": "Europe",
    "Israel": "Middle East", "Saudi Arabia": "Middle East",
    "Iran": "Middle East", "Oman": "Middle East", "United Arab Emirates": "Middle East",
}

# Spelling variants collapsed to one label per country. Country of acquisition
# is a model STATE; two spellings of one country are two states, which silently
# splits the transition matrix.
COUNTRY_ALIASES = {
    "viet nam": "Vietnam",
    "vietnam": "Vietnam",
    "united states": "USA",
    "united states of america": "USA",
    "us": "USA",
    "u.s.a.": "USA",
    "new caledodia": "New Caledonia",
    "korea, south": "South Korea",
    "republic of korea": "South Korea",
    "lao pdr": "Laos",
    "laos": "Laos",
    "burma": "Myanmar",
    "czechia": "Czech Republic",
    "uk": "United Kingdom",
    "united kingdom of great britain and northern ireland": "United Kingdom",
}

# Labels that name more than one country, or no country at all. These cannot be
# a discrete state; they are recorded and excluded, not guessed.
AMBIGUOUS_COUNTRY = {"africa", "panama and peru", "unknown", "", "nan", "none"}

_UNKNOWN = {"", "nan", "none", "null", "unknown", "na", "n/a", "not collected",
            "missing", "not applicable"}


def is_blank(value) -> bool:
    return str(value).strip().lower() in _UNKNOWN


# --------------------------------------------------------------------------
# Newick
# --------------------------------------------------------------------------
def read_newick_tips(path: Path) -> list[str]:
    """Return leaf labels in file order.

    Deliberately dependency-free: a leaf label is the token that sits directly
    after '(' or ',', an internal label the token after ')'. Quoted labels are
    honoured, and [&...] comments (BEAST/figtree annotations) are stripped.
    """
    text = re.sub(r"\[[^\]]*\]", "", path.read_text())
    tips: list[str] = []
    label: list[str] = []
    # What preceded the label being read: a label after '(' or ',' is a leaf,
    # a label after ')' is an internal node.
    preceded_by = ""
    in_branch_length = False
    i = 0
    n = len(text)

    while i < n:
        ch = text[i]
        if ch in "'\"" and not in_branch_length:
            quote = ch
            i += 1
            while i < n and text[i] != quote:
                label.append(text[i])
                i += 1
            i += 1
            continue
        if ch == ":":
            in_branch_length = True
        elif ch in "(),;":
            name = "".join(label).strip()
            if name and preceded_by in ("(", ","):
                tips.append(name)
            label = []
            in_branch_length = False
            preceded_by = ch if ch in "(,)" else preceded_by
        elif not in_branch_length and not ch.isspace():
            label.append(ch)
        i += 1
    return tips


# --------------------------------------------------------------------------
# Identifier normalisation
# --------------------------------------------------------------------------
_FASTA_EXT = re.compile(r"\.(fa|fas|fasta|fna|fsa)$", re.IGNORECASE)
_ASSEMBLY = re.compile(r"^GC[AF][_-]?(\d+)")
_SRA = re.compile(r"^((?:SRR|ERR|DRR)\d+)")


def join_key(label) -> str:
    """Collapse a tip label or a sample_id to one comparable key.

    - `GCF_000756925_1_Australia_Townsville_Queensland` -> `ACC000756925`
      (drops the assembly version and any appended free text, and makes the
      GCF/GCA prefixes interchangeable -- they denote the same assembly and are
      used inconsistently between the tree and the metadata sheet).
    - `SRR33188703_SPAdes`, `SRR33188703_S` -> `SRR33188703` (assembler suffixes
      are truncated differently in different files).
    - `IP-0009-1-R.fas` -> `IP-0009-1-R`.
    """
    s = _FASTA_EXT.sub("", str(label).strip())
    m = _ASSEMBLY.match(s)
    if m:
        return "ACC" + m.group(1)
    m = _SRA.match(s)
    if m:
        return m.group(1)
    return s


def canonical_country(value) -> str:
    if is_blank(value):
        return ""
    text = str(value).strip()
    return COUNTRY_ALIASES.get(text.lower(), text)


# --------------------------------------------------------------------------
# Dates
# --------------------------------------------------------------------------
_DAYS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def _leap(year: int) -> bool:
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


def _decimal(year: int, month: int, day: int) -> float:
    length = 366.0 if _leap(year) else 365.0
    doy = sum(_DAYS[: month - 1]) + (1 if (_leap(year) and month > 2) else 0) + day - 1
    return year + doy / length


def parse_date(raw):
    """-> (precision, year, decimal_midpoint, lower_bound, upper_bound).

    A year-only record is returned as the interval [YYYY.0, YYYY+1.0) with its
    midpoint, not as 1 January. Collapsing an interval to its start biases every
    root-to-tip regression and every tip-dated clock towards older dates.
    """
    if is_blank(raw):
        return "unknown", None, None, None, None
    text = str(raw).strip()[:10]
    m = re.match(r"^(\d{4})-(\d{2})-(\d{2})$", text)
    if m:
        y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
        if 1 <= mo <= 12 and 1 <= d <= 31:
            dec = _decimal(y, mo, d)
            return "day", y, round(dec, 5), round(dec, 5), round(dec, 5)
    m = re.match(r"^(\d{4})-(\d{2})$", text)
    if m:
        y, mo = int(m.group(1)), int(m.group(2))
        if 1 <= mo <= 12:
            lo = _decimal(y, mo, 1)
            days = _DAYS[mo - 1] + (1 if (mo == 2 and _leap(y)) else 0)
            hi = _decimal(y, mo, days)
            return "month", y, round((lo + hi) / 2, 5), round(lo, 5), round(hi, 5)
    m = re.match(r"^(\d{4})$", text)
    if m:
        y = int(m.group(1))
        return "year", y, y + 0.5, float(y), y + 0.99999
    return "unparsed", None, None, None, None


# --------------------------------------------------------------------------
def build(args) -> int:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- tips -------------------------------------------------------------
    if args.tree:
        tips = read_newick_tips(Path(args.tree))
        tip_source = f"newick:{args.tree}"
    else:
        spec = args.tips_from_csv
        path, _, column = spec.partition(":")
        column = column or "tip"
        frame = pd.read_csv(path, dtype=str)
        if column not in frame.columns:
            sys.exit(f"ERROR: column '{column}' not in {path}")
        tips = frame[column].dropna().astype(str).tolist()
        tip_source = f"csv:{path}:{column}"

    if not tips:
        sys.exit("ERROR: no tip labels found")
    duplicate_tips = [t for t, c in Counter(tips).items() if c > 1]

    # ---- metadata ---------------------------------------------------------
    sep = "\t" if args.metadata.lower().endswith((".tsv", ".tab", ".txt")) else ","
    meta = pd.read_csv(args.metadata, sep=sep, dtype=str)
    meta.columns = [c.strip() for c in meta.columns]
    for column in (args.id_col, args.country_col):
        if column not in meta.columns:
            sys.exit(f"ERROR: --metadata has no column '{column}'")

    meta["_key"] = meta[args.id_col].map(join_key)
    key_dupes = meta.loc[meta["_key"].duplicated(keep=False), [args.id_col, "_key"]]
    if not key_dupes.empty:
        sys.exit(
            "ERROR: metadata sample_ids collide after normalisation -- the join "
            f"would be ambiguous:\n{key_dupes.to_string(index=False)}"
        )

    table = pd.DataFrame({"tip": tips})
    table["join_key"] = table["tip"].map(join_key)
    lookup = meta.set_index("_key")
    table["matched"] = table["join_key"].isin(lookup.index)

    unmatched = table.loc[~table["matched"], ["tip", "join_key"]]
    if not unmatched.empty:
        unmatched.to_csv(outdir / "unmatched_tips.tsv", sep="\t", index=False)
        if not args.allow_unmatched:
            sys.exit(
                f"ERROR: {len(unmatched)} of {len(table)} tips have no metadata row. "
                f"See {outdir / 'unmatched_tips.tsv'}. Fix the identifiers, or pass "
                "--allow-unmatched to proceed with those tips blank."
            )

    joined = table.merge(
        meta.rename(columns={"_key": "join_key"}), on="join_key", how="left"
    )

    # ---- geography --------------------------------------------------------
    joined["sample_id"] = joined[args.id_col]
    joined["country_acquired"] = joined[args.country_col].map(canonical_country)
    if args.diagnosis_col in joined.columns:
        joined["country_diagnosis"] = joined[args.diagnosis_col].map(canonical_country)
    else:
        joined["country_diagnosis"] = ""

    joined["country_ambiguous"] = joined["country_acquired"].str.lower().isin(
        AMBIGUOUS_COUNTRY
    )
    joined["travel_associated"] = (
        (joined["country_diagnosis"] != "")
        & (joined["country_acquired"] != "")
        & (joined["country_diagnosis"] != joined["country_acquired"])
    )
    joined["geo_region"] = joined["country_acquired"].map(
        lambda c: REGION.get(c, "" if not c else "UNMAPPED")
    )
    unmapped = sorted(set(joined.loc[joined["geo_region"] == "UNMAPPED", "country_acquired"]))
    if unmapped and not args.allow_unmapped_region:
        sys.exit(
            "ERROR: no region for: " + ", ".join(unmapped) + ". Add them to REGION in "
            "this script (or pass --allow-unmapped-region)."
        )

    # geo_state: the trait a discrete phylogeographic model actually consumes.
    joined["geo_state"] = joined["country_acquired"].where(~joined["country_ambiguous"], "")
    counts = joined.loc[joined["geo_state"] != "", "geo_state"].value_counts()
    rare = set(counts[counts < args.min_state_count].index)
    # A tip whose country is rare, or too vague to be a state at all ("Africa",
    # "Panama and Peru"), still has a defensible REGION. Pooling recovers it for
    # the coarse analysis instead of discarding the genome.
    joined["geo_state_pooled"] = [
        region if (state == "" or state in rare) and region not in ("", "UNMAPPED") else state
        for state, region in zip(joined["geo_state"], joined["geo_region"])
    ]

    # admin-1: prefer the ISO-3166-2 name, which aggregates districts to their
    # province; the free-text Subregion mixes administrative levels.
    def admin1(row):
        for column, label in ((args.admin1_col, "iso_3166_2"), (args.subregion_col, "subregion")):
            if column in row.index and not is_blank(row[column]):
                return pd.Series([str(row[column]).strip(), label])
        return pd.Series(["", ""])

    joined[["admin1", "admin1_source"]] = joined.apply(admin1, axis=1)
    joined["thai_province"] = joined["admin1"].where(
        joined["country_acquired"] == "Thailand", ""
    )

    # ---- dates ------------------------------------------------------------
    parsed = joined[args.date_col].map(parse_date) if args.date_col in joined.columns else None
    if parsed is None:
        sys.exit(f"ERROR: --metadata has no column '{args.date_col}'")
    joined["date_raw"] = joined[args.date_col]
    joined["date_precision"] = [p[0] for p in parsed]
    joined["year"] = [p[1] for p in parsed]
    joined["date_decimal"] = [p[2] for p in parsed]
    joined["date_lower"] = [p[3] for p in parsed]
    joined["date_upper"] = [p[4] for p in parsed]
    joined["decade"] = [f"{int(y) // 10 * 10}s" if pd.notna(y) else "" for y in joined["year"]]
    joined["clock_usable"] = joined["date_decimal"].notna()

    # ---- source -----------------------------------------------------------
    def source(row):
        if args.envi_col in row.index and not is_blank(row.get(args.envi_col)):
            return "Environmental"
        if args.patient_col in row.index and not is_blank(row.get(args.patient_col)):
            return "Clinical"
        sample = str(row["sample_id"])
        if sample.upper().startswith("IE-"):
            return "Environmental"
        if sample.upper().startswith("IP-"):
            return "Clinical"
        return "Unspecified"

    joined["source"] = joined.apply(source, axis=1)
    joined["collection"] = [
        "In-house" if s != "Unspecified" else "Public" for s in joined["source"]
    ]

    # ---- cluster ----------------------------------------------------------
    joined["cluster_id"] = ""
    if args.clusters:
        clusters = pd.read_csv(args.clusters, sep="\t", dtype=str)
        clusters["join_key"] = clusters[clusters.columns[1]].map(join_key)
        cluster_map = dict(zip(clusters["join_key"], clusters[clusters.columns[0]]))
        joined["cluster_id"] = joined["join_key"].map(cluster_map).fillna("")
    if args.cluster_qc:
        qc = pd.read_csv(args.cluster_qc, dtype=str)
        for column in ("confidence_tier", "n_recombination_blocks",
                       "n_filtered_polymorphic_sites", "n_isolates"):
            if column in qc.columns:
                joined["cluster_" + column] = joined["cluster_id"].map(
                    dict(zip(qc["cluster_id"], qc[column]))
                ).fillna("")

    # ---- carry-forward ----------------------------------------------------
    carried: list[str] = []
    if args.carry:
        path, _, columns = args.carry.partition(":")
        wanted = [c for c in columns.split(",") if c] or ["Lineage", "Lineage_label"]
        prior = pd.read_csv(path, dtype=str)
        tip_column = "tip" if "tip" in prior.columns else prior.columns[0]
        prior["join_key"] = prior[tip_column].map(join_key)
        # Namespaced, so a carried column can never silently overwrite one this
        # script derives (a prior `Source` vs the `source` recomputed here).
        for column in wanted:
            if column in prior.columns:
                joined["prior_" + column.lower()] = joined["join_key"].map(
                    dict(zip(prior["join_key"], prior[column]))
                ).fillna("")
        carried[:] = ["prior_" + c.lower() for c in wanted if c in prior.columns]

    # ---- write ------------------------------------------------------------
    columns = [
        "tip", "sample_id", "join_key", "matched",
        "country_acquired", "country_diagnosis", "travel_associated",
        "country_ambiguous", "geo_region", "geo_state", "geo_state_pooled",
        "admin1", "admin1_source", "thai_province",
        "date_raw", "date_precision", "year", "decade",
        "date_decimal", "date_lower", "date_upper", "clock_usable",
        "source", "collection", "cluster_id",
    ]
    columns += [c for c in joined.columns if c.startswith("cluster_") and c != "cluster_id"]
    columns += carried
    passthrough = [c for c in args.keep.split(",") if c and c in joined.columns]
    columns += passthrough
    out = joined[[c for c in columns if c in joined.columns]]
    out.to_csv(outdir / "tip_annotation.tsv", sep="\t", index=False)

    # TreeTime: midpoint plus bounds, so imprecise dates stay imprecise.
    clock = joined[joined["clock_usable"]]
    clock[["tip", "date_decimal", "date_lower", "date_upper"]].rename(
        columns={"tip": "name", "date_decimal": "date",
                 "date_lower": "date_lower_bound", "date_upper": "date_upper_bound"}
    ).to_csv(outdir / "tip_dates_treetime.csv", index=False)

    # LSD2: b(lo,hi) for anything coarser than a day.
    with open(outdir / "tip_dates_lsd2.txt", "w") as handle:
        handle.write(f"{len(clock)}\n")
        for _, row in clock.iterrows():
            if row["date_precision"] == "day":
                handle.write(f"{row['tip']}\t{row['date_decimal']:.5f}\n")
            else:
                handle.write(f"{row['tip']}\tb({row['date_lower']:.5f},{row['date_upper']:.5f})\n")

    # Discrete traits for mugration / DTA.
    joined[["tip", "geo_state", "geo_state_pooled", "geo_region", "thai_province"]].rename(
        columns={"tip": "name"}
    ).to_csv(outdir / "traits_geography.tsv", sep="\t", index=False)

    # Every isolate whose country the model will treat differently from the
    # country it was diagnosed in.
    changed = joined[joined["travel_associated"] | joined["country_ambiguous"]]
    changed[["tip", "sample_id", "country_diagnosis", "country_acquired",
             "geo_region", "country_ambiguous"]].to_csv(
        outdir / "country_change_log.tsv", sep="\t", index=False
    )

    # ---- QC ---------------------------------------------------------------
    lines = [
        "tip annotation QC",
        "=================",
        f"tip source              : {tip_source}",
        f"metadata                : {args.metadata}",
        f"country of acquisition  : {args.country_col}",
        f"country of diagnosis    : {args.diagnosis_col}",
        "",
        f"tips                    : {len(table)}",
        f"duplicate tip labels    : {len(duplicate_tips)}"
        + (f"  {duplicate_tips[:10]}" if duplicate_tips else ""),
        f"matched to metadata     : {int(table['matched'].sum())} "
        f"({100 * table['matched'].mean():.2f}%)",
        f"unmatched               : {len(unmatched)}",
        f"metadata rows           : {len(meta)}",
        f"metadata rows unused    : {len(meta) - int(table['matched'].sum())}",
        "",
        "geography",
        f"  countries (acquired)  : {joined['country_acquired'].replace('', pd.NA).nunique()}",
        f"  ambiguous country     : {int(joined['country_ambiguous'].sum())}",
        f"  travel-associated     : {int(joined['travel_associated'].sum())}",
        f"  usable geo_state      : {int((joined['geo_state'] != '').sum())}",
        f"  states after pooling  : {joined.loc[joined['geo_state_pooled'] != '', 'geo_state_pooled'].nunique()}"
        f" (min {args.min_state_count} tips)",
        f"  admin1 resolved       : {int((joined['admin1'] != '').sum())}",
        "",
        "dates",
    ]
    for precision, count in joined["date_precision"].value_counts().items():
        lines.append(f"  {precision:<21} : {count}")
    lines += [
        f"  clock-usable          : {int(joined['clock_usable'].sum())} "
        f"({100 * joined['clock_usable'].mean():.1f}%)",
        f"  year range            : {joined['year'].min()} - {joined['year'].max()}",
        "",
        "source",
    ]
    for label, count in joined["source"].value_counts().items():
        lines.append(f"  {label:<21} : {count}")
    if unmapped:
        lines += ["", "countries with no region mapping: " + ", ".join(unmapped)]
    lines += [
        "",
        "country distribution (acquisition):",
    ]
    for country, count in joined["country_acquired"].value_counts().items():
        lines.append(f"  {country:<24} {count}")

    (outdir / "tip_annotation_QC.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwrote -> {outdir}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Map metadata onto tree tips and emit phylogeography / clock inputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--tree", help="Newick tree whose tips are to be annotated")
    source.add_argument("--tips-from-csv", help="FILE[:COLUMN] holding the tip labels")

    parser.add_argument("--metadata", required=True, help="metadata TSV or CSV")
    parser.add_argument("--outdir", required=True)

    parser.add_argument("--id-col", default="sample_id")
    parser.add_argument("--country-col", default="Country_Final",
                        help="country of ACQUISITION -- the state a phylogeographic model should use")
    parser.add_argument("--diagnosis-col", default="Country_Diagnosis",
                        help="country of DIAGNOSIS -- retained for audit, not used as the state")
    parser.add_argument("--date-col", default="final_collection_dates")
    parser.add_argument("--admin1-col", default="iso_3166_2_name")
    parser.add_argument("--subregion-col", default="Subregion")
    parser.add_argument("--envi-col", default="EnviSampleID")
    parser.add_argument("--patient-col", default="Linked_PatientCase")

    parser.add_argument("--clusters", help="clusters.tsv (cluster_id, sample_id)")
    parser.add_argument("--cluster-qc", help="cluster_phylogeny_summary.csv")
    parser.add_argument("--carry", help="FILE:COL1,COL2 to carry forward by tip "
                                        "(e.g. annotated_tips.csv:Lineage,Lineage_label)")
    parser.add_argument("--keep", default="",
                        help="comma-separated extra metadata columns to pass through")

    parser.add_argument("--min-state-count", type=int, default=5,
                        help="countries with fewer tips are pooled to their region in geo_state_pooled")
    parser.add_argument("--allow-unmatched", action="store_true",
                        help="do not abort when a tip has no metadata row")
    parser.add_argument("--allow-unmapped-region", action="store_true")
    return build(parser.parse_args())


if __name__ == "__main__":
    sys.exit(main())
