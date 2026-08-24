#!/usr/bin/env python3
"""Build draft FASTQ/ENA metadata tables from legacy sequencing records.

The source FASTQ files are not present in this repository. This script gathers
the information that can be recovered from old metadata, submission forms, lane
statistics, and the previous supplementary QC table. Fields that require manual
curation before ENA submission are left blank and flagged.
"""

from __future__ import annotations

import re
from pathlib import Path
from zipfile import BadZipFile

import pandas as pd
from docx import Document


ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT.parent / "old"
OUTDIR = ROOT / "data" / "metadata" / "ena"


def project_source(path: Path) -> str:
    return str(path.relative_to(ROOT.parent))


def clean_string(value):
    if pd.isna(value):
        return pd.NA
    value = str(value).strip()
    if value == "" or value.lower() in {"nan", "none"}:
        return pd.NA
    return value


def read_sample_rename() -> pd.DataFrame:
    path = OLD / "sample_rename.txt"
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            rows.append({"legacy_analysis_id": parts[0], "sample_id": parts[1]})
    return pd.DataFrame(rows)


def read_old_snowy_metadata() -> pd.DataFrame:
    path = OLD / "Metadata_curated_owls.xlsx"
    df = pd.read_excel(path, sheet_name="snowy_owl", dtype=str)
    df = df.rename(
        columns={
            "Sample_ID": "sample_short_id",
            "Seq ID": "old_seq_id",
            "Sequenced (date)": "sequenced_date_old_metadata",
            "VCF_ID": "old_vcf_id",
            "Coverage 1st round": "coverage_first_round_old_metadata",
        }
    )
    keep = [
        "sample_short_id",
        "old_seq_id",
        "sequenced_date_old_metadata",
        "old_vcf_id",
        "coverage_first_round_old_metadata",
    ]
    return df[keep].map(clean_string)


def read_metadata2() -> pd.DataFrame:
    path = OLD / "script (1)" / "metadata2.xlsx"
    df = pd.read_excel(path, dtype=str)
    df = df.rename(
        columns={
            "Sample_ID": "sample_short_id",
            "ID": "metadata2_old_seq_id",
            "Seq_round": "sequencing_round_metadata2",
            "seq_depth": "seq_depth_metadata2",
        }
    )
    keep = [
        "sample_short_id",
        "metadata2_old_seq_id",
        "sequencing_round_metadata2",
        "seq_depth_metadata2",
    ]
    return df[keep].map(clean_string)


def read_lane_stats() -> pd.DataFrame:
    path = OLD / "lab_seq" / "snowy-owls_lanestats.xlsx"
    df = pd.read_excel(path, dtype=str)
    df.columns = [re.sub(r"\s+", "_", c.strip().lower()) for c in df.columns]
    df = df.rename(
        columns={
            "omrf_name": "lane_omrf_name",
            "q30": "lane_q30_fraction",
            "%_lane": "lane_fraction",
            "ttl_%_all": "lane_total_fraction_all",
            "num_reads)": "lane_num_reads",
            "gb": "lane_yield_gb",
            "submission_name": "lane_submission_name",
            "index1": "index1",
            "index2": "index2",
        }
    )
    for col in df.columns:
        df[col] = df[col].map(clean_string)
    df["lane_library_prefix"] = df["lane_submission_name"].map(
        lambda x: pd.NA if pd.isna(x) else str(x).strip()
    )
    df["source_file_lane_stats"] = project_source(path)
    return df


def read_submission_table(path: Path, batch_name: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="sample table template", dtype=str)
    df = df.rename(
        columns={
            "Sample Name (only letters, numbers and hyphen allowed, MAX 16 characters)": "submission_sample_name",
            "Type   (see guidelines for codes)": "submitted_material_or_library_type",
            "Conc. (ng/µl)": "submitted_concentration_ng_ul",
            "Volume provided (µl)": "submitted_volume_ul",
            "Total DNA / RNA (µg)": "submitted_total_dna_rna_ug",
            "Approx no. Reads, Gb or lanes requested": "requested_output",
            "Comments": "submission_comments",
        }
    )
    keep = [
        "submission_sample_name",
        "submitted_material_or_library_type",
        "submitted_concentration_ng_ul",
        "submitted_volume_ul",
        "submitted_total_dna_rna_ug",
        "requested_output",
        "submission_comments",
    ]
    df = df[keep].map(clean_string)
    df = df[df["submission_sample_name"].notna()].copy()
    df = df[
        ~df["submission_sample_name"].str.contains("example|buffer", case=False, na=False)
    ].copy()
    df["submission_batch"] = batch_name
    df["source_file_submission_table"] = project_source(path)
    df["submission_sample_code"] = (
        df["submission_sample_name"]
        .str.replace(r"^Sample\s+", "", regex=True)
        .str.replace(r"\s+", "", regex=True)
    )
    return df


def read_all_submission_tables() -> pd.DataFrame:
    tables = [
        (
            OLD / "lab_seq" / "Baalsrud_sample_table_SnowyOwl_22112023.xlsx",
            "NSC_2023_Baalsrud_snowy_owl",
        ),
        (
            OLD / "lab_seq" / "Boessenkool_sample_table_96_SnowyOwl_25052022.xlsx",
            "NSC_2022_Boessenkool_snowy_owl",
        ),
    ]
    return pd.concat([read_submission_table(*x) for x in tables], ignore_index=True)


def read_alignment_qc_from_docx() -> pd.DataFrame:
    path = OLD / "supporting_info" / "Supplementary data.docx"
    try:
        doc = Document(path)
    except (BadZipFile, FileNotFoundError):
        return pd.DataFrame()

    rows = []
    for table in doc.tables:
        table_rows = [
            [cell.text.strip().replace("\n", " ") for cell in row.cells]
            for row in table.rows
        ]
        if not table_rows:
            continue
        header = table_rows[0]
        if "Sample Name" not in header or "% Aligned" not in header:
            continue
        for row in table_rows[1:]:
            if len(row) != len(header) or not any(row):
                continue
            rows.append(dict(zip(header, row)))

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df = df.rename(
        columns={
            "Sample Name": "alignment_qc_id",
            "% GC": "gc_percent",
            "Ins. size": "insert_size_bp",
            "≥ 30X": "percent_ge_30x",
            "Median cov": "median_coverage",
            "Mean cov": "mean_coverage",
            "% Aligned": "percent_aligned",
        }
    )
    for col in df.columns:
        df[col] = df[col].map(clean_string)
    df["source_file_alignment_qc"] = project_source(path)
    return df


def infer_submission_match(row: pd.Series, lookup: pd.DataFrame) -> pd.NA | str:
    code = row["submission_sample_code"]
    if pd.isna(code):
        return pd.NA
    code = str(code)
    if re.fullmatch(r"X\d+", code):
        n = code[1:]
        patterns = [
            f"merged_X{n}",
            f"-X{n}_",
            f"-X{n}.",
            f"X{n}_",
        ]
    elif re.fullmatch(r"WRG\d+", code):
        patterns = [f"-{code}_", f"{code}_", f"{code}.", f"_{code}_"]
    else:
        return pd.NA

    hits = lookup[
        lookup[["library_id", "old_seq_id", "old_vcf_id"]]
        .fillna("")
        .apply(lambda cols: any(p in " ".join(cols) for p in patterns), axis=1)
    ]
    if len(hits) == 1:
        return hits.iloc[0]["sample_id"]
    if len(hits) > 1:
        return ";".join(hits["sample_id"].tolist())
    return pd.NA


def first_nonmissing(values):
    values = [v for v in values if pd.notna(v)]
    return values[0] if values else pd.NA


def collapse_values(values):
    values = [str(v) for v in values if pd.notna(v) and str(v).strip()]
    if not values:
        return pd.NA
    return "; ".join(dict.fromkeys(values))


def infer_folder_date(folder: str) -> str | pd.NA:
    match = re.search(r"(\d{4})-(\d{2})-(\d{2})$", folder)
    if not match:
        return pd.NA
    year, first, second = match.groups()
    if int(first) > 12:
        return f"{year}-{second}-{first}"
    return f"{year}-{first}-{second}"


def infer_trimmed_match_code(sample_file_id: str) -> str | pd.NA:
    if sample_file_id.startswith("Undetermined"):
        return "Undetermined"
    dag = re.search(r"(DAG-BS-\d+)", sample_file_id)
    if dag:
        return dag.group(1)
    wrg = re.search(r"(WRG\d+)", sample_file_id)
    if wrg:
        return wrg.group(1)
    x_code = re.search(r"-?(X\d+)_", sample_file_id)
    if x_code:
        return x_code.group(1)
    return pd.NA


def infer_sample_match_code(row: pd.Series) -> str | pd.NA:
    fields = [
        row.get("library_id"),
        row.get("old_seq_id"),
        row.get("metadata2_old_seq_id"),
        row.get("old_vcf_id"),
    ]
    text = " ".join(str(x) for x in fields if pd.notna(x))
    dag = re.search(r"(DAG-BS-\d+)", text)
    if dag:
        return dag.group(1)
    wrg = re.search(r"(WRG\d+)", text)
    if wrg:
        return wrg.group(1)
    x_code = re.search(r"(?:merged_)?X(\d+)", text)
    if x_code:
        return f"X{x_code.group(1)}"
    return pd.NA


def read_trimmed_fastq_listing() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = OUTDIR / "trimmed_fastq_filenames.txt"
    if not path.exists():
        return pd.DataFrame(), pd.DataFrame()

    rows = []
    pattern = re.compile(
        r"(?P<folder>0[1-5]_[^/\s]+)/(?P<file>[^\s]+?_R(?P<read>[12])_001_trim_(?P<pair_status>unpair|pair)\.fastq\.gz)"
    )
    for line in path.read_text().splitlines():
        match = pattern.search(line.strip())
        if not match:
            continue
        folder = match.group("folder")
        filename = match.group("file")
        sample_file_id = re.sub(
            r"_R[12]_001_trim_(?:unpair|pair)\.fastq\.gz$", "", filename
        )
        rows.append(
            {
                "trimmed_sequencing_round_folder": folder,
                "trimmed_sequencing_round_number": folder.split("_", 1)[0],
                "trimmed_sequencing_date_inferred": infer_folder_date(folder),
                "trimmed_sample_file_id": sample_file_id,
                "trimmed_match_code": infer_trimmed_match_code(sample_file_id),
                "trimmed_read": f"R{match.group('read')}",
                "trimmed_pair_status": match.group("pair_status"),
                "trimmed_fastq_file": f"{folder}/{filename}",
            }
        )

    files = pd.DataFrame(rows)
    if files.empty:
        return files, pd.DataFrame()

    def joined_filtered(group: pd.DataFrame, read: str, pair_status: str) -> str | pd.NA:
        values = group.loc[
            (group["trimmed_read"] == read)
            & (group["trimmed_pair_status"] == pair_status),
            "trimmed_fastq_file",
        ].tolist()
        return collapse_values(values)

    grouped_rows = []
    for code, group in files.groupby("trimmed_match_code", dropna=False):
        if pd.isna(code):
            continue
        pair_groups = group[group["trimmed_pair_status"] == "pair"].groupby(
            ["trimmed_sequencing_round_folder", "trimmed_sample_file_id"], dropna=False
        )
        grouped_rows.append(
            {
                "trimmed_match_code": code,
                "trimmed_sample_file_ids": collapse_values(
                    group["trimmed_sample_file_id"].tolist()
                ),
                "trimmed_sequencing_round_folders": collapse_values(
                    group["trimmed_sequencing_round_folder"].tolist()
                ),
                "trimmed_sequencing_round_numbers": collapse_values(
                    group["trimmed_sequencing_round_number"].tolist()
                ),
                "trimmed_sequencing_dates_inferred": collapse_values(
                    group["trimmed_sequencing_date_inferred"].tolist()
                ),
                "n_trimmed_rounds_found": group[
                    "trimmed_sequencing_round_folder"
                ].nunique(),
                "n_trimmed_fastq_pairs_found": sum(
                    int((sub["trimmed_read"] == "R1").any())
                    and int((sub["trimmed_read"] == "R2").any())
                    for _, sub in pair_groups
                ),
                "trimmed_pair_r1_files": joined_filtered(group, "R1", "pair"),
                "trimmed_pair_r2_files": joined_filtered(group, "R2", "pair"),
                "trimmed_unpaired_r1_files": joined_filtered(group, "R1", "unpair"),
                "trimmed_unpaired_r2_files": joined_filtered(group, "R2", "unpair"),
                "source_file_trimmed_listing": str(path.relative_to(ROOT)),
            }
        )

    return files, pd.DataFrame(grouped_rows)


def add_manual_flags(row: pd.Series) -> str:
    flags = ["add actual FASTQ R1/R2 filenames and ENA run/accession IDs"]
    if pd.notna(row.get("trimmed_pair_r1_files")):
        flags.append("trimmed FASTQ evidence found: replace/confirm with raw FASTQs for ENA")
    if pd.notna(row.get("n_trimmed_fastq_pairs_found")):
        try:
            if int(row.get("n_trimmed_fastq_pairs_found")) > 1:
                flags.append("multiple trimmed FASTQ pairs/rounds found for this sample")
        except ValueError:
            pass
    if str(row.get("library_id", "")).startswith("merged_") or str(
        row.get("old_vcf_id", "")
    ).startswith("merged_"):
        flags.append("merged analysis ID: confirm all constituent FASTQ/lane files")
    if pd.isna(row.get("source_file_submission_table")):
        flags.append("no individual submission-table row recovered")
    if pd.isna(row.get("lane_omrf_name")):
        flags.append("no lane-stat row recovered")
    if pd.isna(row.get("percent_aligned")):
        flags.append("no historical alignment-QC row or blank QC row")
    if str(row.get("submitted_material_or_library_type", "")).lower().find("gdna") >= 0:
        flags.append("gDNA submitted: confirm library prep kit/protocol")
    if pd.notna(row.get("submission_comments")) and "Thruplex" in str(
        row.get("submission_comments")
    ):
        flags.append("Thruplex noted in comments: confirm exact kit/version")
    if pd.notna(row.get("submitted_material_or_library_type")) and "Truseq" in str(
        row.get("submitted_material_or_library_type")
    ):
        flags.append("TruSeq library noted: confirm exact kit/version")
    return "; ".join(flags)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    samples = pd.read_csv(ROOT / "data" / "metadata" / "samples_curated.csv", dtype=str)
    samples = samples.map(clean_string)
    old_meta = read_old_snowy_metadata()
    metadata2 = read_metadata2()
    renames = read_sample_rename()
    lanes = read_lane_stats()
    submissions = read_all_submission_tables()
    alignment = read_alignment_qc_from_docx()
    trimmed_files, trimmed_by_code = read_trimmed_fastq_listing()

    lookup = samples.merge(old_meta, how="left", on="sample_short_id")
    submissions["sample_id"] = submissions.apply(
        lambda row: infer_submission_match(row, lookup), axis=1
    )
    submissions = submissions[
        submissions["sample_id"].notna() & ~submissions["sample_id"].str.contains(";")
    ].copy()

    submission_by_sample = (
        submissions.groupby("sample_id", dropna=False)
        .agg({col: collapse_values for col in submissions.columns if col != "sample_id"})
        .reset_index()
    )

    lanes["lane_library_id"] = lanes["lane_submission_name"].map(
        lambda x: pd.NA if pd.isna(x) else f"{x}_S{x.split('-')[-1]}" if re.fullmatch(r"DAG-BS-\d+", x) else x
    )
    lanes_by_library = lanes.copy()

    alignment = alignment.merge(
        renames.rename(columns={"legacy_analysis_id": "alignment_qc_id"}),
        how="left",
        on="alignment_qc_id",
    )

    final = (
        samples.merge(old_meta, how="left", on="sample_short_id")
        .merge(metadata2, how="left", on="sample_short_id")
        .merge(submission_by_sample, how="left", on="sample_id")
        .merge(
            lanes_by_library,
            how="left",
            left_on="library_id",
            right_on="lane_library_id",
        )
        .merge(
            alignment.drop(columns=["sample_id"]).rename(
                columns={"alignment_qc_id": "alignment_qc_legacy_id"}
            ),
            how="left",
            left_on="library_id",
            right_on="alignment_qc_legacy_id",
        )
    )
    final["trimmed_match_code"] = final.apply(infer_sample_match_code, axis=1)
    if not trimmed_by_code.empty:
        final = final.merge(trimmed_by_code, how="left", on="trimmed_match_code")

    final["fastq_r1_file"] = pd.NA
    final["fastq_r2_file"] = pd.NA
    final["ena_run_accession"] = pd.NA
    final["ena_experiment_accession"] = pd.NA
    final["paired_end_expected"] = "yes"
    final["instrument_platform_inferred"] = "Illumina"
    final["instrument_model_or_run"] = pd.NA
    final["sequencing_centre_or_provider_inferred"] = pd.NA

    nsc_mask = final["submission_batch"].fillna("").str.contains("NSC")
    final.loc[nsc_mask, "sequencing_centre_or_provider_inferred"] = (
        "Norwegian Sequencing Centre (from submission table)"
    )
    final.loc[final["lane_omrf_name"].notna(), "sequencing_centre_or_provider_inferred"] = (
        "OMRF/USA batch inferred from lane-stats file; verify provider"
    )

    final["manual_check_flags"] = final.apply(add_manual_flags, axis=1)

    ordered = [
        "sample_id",
        "sample_short_id",
        "individual_id",
        "country",
        "location",
        "population_final",
        "library_id",
        "old_seq_id",
        "metadata2_old_seq_id",
        "old_vcf_id",
        "sequenced_date_old_metadata",
        "sequencing_round_metadata2",
        "seq_depth_metadata2",
        "mean_depth",
        "submission_batch",
        "submission_sample_name",
        "submitted_material_or_library_type",
        "submission_comments",
        "requested_output",
        "submitted_concentration_ng_ul",
        "submitted_volume_ul",
        "submitted_total_dna_rna_ug",
        "sequencing_centre_or_provider_inferred",
        "instrument_platform_inferred",
        "instrument_model_or_run",
        "paired_end_expected",
        "fastq_r1_file",
        "fastq_r2_file",
        "trimmed_match_code",
        "trimmed_sample_file_ids",
        "trimmed_sequencing_round_folders",
        "trimmed_sequencing_round_numbers",
        "trimmed_sequencing_dates_inferred",
        "n_trimmed_rounds_found",
        "n_trimmed_fastq_pairs_found",
        "trimmed_pair_r1_files",
        "trimmed_pair_r2_files",
        "trimmed_unpaired_r1_files",
        "trimmed_unpaired_r2_files",
        "ena_run_accession",
        "ena_experiment_accession",
        "lane_omrf_name",
        "lane_submission_name",
        "lane_q30_fraction",
        "lane_fraction",
        "lane_num_reads",
        "lane_yield_gb",
        "index1",
        "index2",
        "alignment_qc_legacy_id",
        "gc_percent",
        "insert_size_bp",
        "percent_ge_30x",
        "median_coverage",
        "mean_coverage",
        "percent_aligned",
        "source_file_submission_table",
        "source_file_lane_stats",
        "source_file_alignment_qc",
        "source_file_trimmed_listing",
        "manual_check_flags",
    ]
    final = final[[c for c in ordered if c in final.columns]]

    final.to_csv(OUTDIR / "fastq_ena_overview_draft.csv", index=False)
    lanes.to_csv(OUTDIR / "sequencing_lane_stats_clean.csv", index=False)
    alignment.to_csv(OUTDIR / "historical_alignment_qc_from_old_supplement.csv", index=False)
    submissions.to_csv(OUTDIR / "sample_submission_rows_clean.csv", index=False)
    if not trimmed_files.empty:
        trimmed_files.to_csv(OUTDIR / "trimmed_fastq_files_clean.csv", index=False)
    if not trimmed_by_code.empty:
        trimmed_by_code.to_csv(
            OUTDIR / "trimmed_fastq_by_match_code.csv", index=False
        )

    print(f"Wrote {OUTDIR / 'fastq_ena_overview_draft.csv'} ({len(final)} samples)")
    print(f"Wrote {len(lanes)} lane-stat rows and {len(submissions)} matched submission rows")
    print(f"Wrote {len(alignment)} historical alignment-QC rows")
    if not trimmed_files.empty:
        print(
            f"Wrote {len(trimmed_files)} trimmed FASTQ file rows for "
            f"{len(trimmed_by_code)} inferred sample codes"
        )
    print("Samples without submission rows:", final["source_file_submission_table"].isna().sum())
    print("Samples without lane stats:", final["lane_omrf_name"].isna().sum())
    print("Samples without alignment-QC rows:", final["percent_aligned"].isna().sum())
    if "trimmed_pair_r1_files" in final.columns:
        print("Samples without trimmed FASTQ evidence:", final["trimmed_pair_r1_files"].isna().sum())


if __name__ == "__main__":
    main()
