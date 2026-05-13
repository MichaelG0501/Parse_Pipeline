#!/usr/bin/env python3
import csv
import gzip
import re
import shutil
import urllib.request
from pathlib import Path

import pandas as pd


WD = Path("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline")
OUT = WD / "parse_outs" / "Auto_velocity"
REF_IN = Path("/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz")
GENOME_PREFIX = "GRCh38-1-1-3c_"
TARGET_SAMPLES = ["T0", "T1", "T2", "T4", "R4", "eR4"]
EXCLUDED_SAMPLES = {"PDO", "SUR1090"}

SUBLIBRARIES = {
    "NACT1090_A": {
        "trailmaker_dir": WD / "data/trailmaker/output_NACT1090_A_MKDL260004725-1A_23JFJWLT3_L6_REP_CLEAN",
        "combined_suffix": "__s2",
        "bam": WD / "data/trailmaker/output_NACT1090_A_MKDL260004725-1A_23JFJWLT3_L6_REP_CLEAN/process/barcode_headAligned_anno.bam",
    },
    "NACT1090_B": {
        "trailmaker_dir": WD / "data/trailmaker/output_NACT1090_B_MKDL260004726-1A_23JFJWLT3_L6_REP_CLEAN",
        "combined_suffix": "__s1",
        "bam": WD / "data/trailmaker/output_NACT1090_B_MKDL260004726-1A_23JFJWLT3_L6_REP_CLEAN/process/barcode_headAligned_anno.bam",
    },
}


def map_chrom(chrom: str) -> str:
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
    if chrom.startswith("Un_"):
        chrom = chrom[3:]
    chrom = re.sub(r"v([0-9]+)$", r".\1", chrom)
    return f"{GENOME_PREFIX}{chrom}"


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def write_prefixed_gene_gtf() -> Path:
    out = OUT / "ref/genes.GRCh38-1-1-3c.gtf"
    if out.exists() and out.stat().st_size > 0:
        return out

    out.parent.mkdir(parents=True, exist_ok=True)
    with open_text(REF_IN) as inp, open(out, "w") as handle:
        for line in inp:
            if line.startswith("#"):
                handle.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            fields[0] = map_chrom(fields[0])
            handle.write("\t".join(fields) + "\n")
    return out


def download_repeatmasker() -> Path:
    raw = OUT / "ref/rmsk.hg38.txt.gz"
    if raw.exists() and raw.stat().st_size > 0:
        return raw

    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
    tmp = raw.with_suffix(".tmp")
    urllib.request.urlretrieve(url, tmp)
    shutil.move(tmp, raw)
    return raw


def write_repeatmasker_gtf() -> Path:
    out = OUT / "ref/repeatmasker.GRCh38-1-1-3c.gtf"
    if out.exists() and out.stat().st_size > 0:
        return out

    raw = download_repeatmasker()
    with gzip.open(raw, "rt") as inp, open(out, "w") as handle:
        handle.write("# UCSC hg38 rmsk converted to GTF and prefixed for Parse Split-pipe GRCh38-1-1-3c BAMs\n")
        for row in csv.reader(inp, delimiter="\t"):
            if len(row) < 17:
                continue
            chrom = map_chrom(row[5])
            start = int(row[6]) + 1
            end = int(row[7])
            strand = row[9] if row[9] in {"+", "-"} else "."
            rep_name = row[10].replace('"', "'")
            rep_class = row[11].replace('"', "'")
            rep_family = row[12].replace('"', "'")
            attrs = (
                f'gene_id "{rep_name}"; transcript_id "{rep_name}"; '
                f'rep_class "{rep_class}"; rep_family "{rep_family}";'
            )
            handle.write(f"{chrom}\tUCSC_rmsk\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")
    return out


def write_barcodes_and_metadata() -> None:
    rows = []
    for sublib, cfg in SUBLIBRARIES.items():
        meta_path = cfg["trailmaker_dir"] / "all-sample/DGE_filtered/cell_metadata.csv.gz"
        meta = pd.read_csv(meta_path)
        meta = meta[meta["sample"].isin(TARGET_SAMPLES)]
        meta = meta[~meta["sample"].isin(EXCLUDED_SAMPLES)].copy()
        meta["sublibrary"] = sublib
        meta["combined_cell_id"] = meta["bc_wells"].astype(str) + cfg["combined_suffix"]
        meta["bam"] = str(cfg["bam"])

        barcode_path = OUT / f"barcodes/{sublib}_nonPDO_barcodes.tsv"
        barcode_path.parent.mkdir(parents=True, exist_ok=True)
        meta["bc_wells"].drop_duplicates().sort_values().to_csv(barcode_path, index=False, header=False)

        meta_path_out = OUT / f"tables/{sublib}_nonPDO_cell_metadata.csv"
        meta_path_out.parent.mkdir(parents=True, exist_ok=True)
        meta.to_csv(meta_path_out, index=False)
        rows.append(meta)

    combined = pd.concat(rows, ignore_index=True)
    combined.to_csv(OUT / "tables/Auto_velocity_nonPDO_cell_metadata.csv", index=False)

    summary = (
        combined.groupby(["sublibrary", "sample"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values(["sublibrary", "sample"])
    )
    summary.to_csv(OUT / "tables/Auto_velocity_input_summary.csv", index=False)


def write_bam_symlinks() -> None:
    tmp = OUT / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    for sublib, cfg in SUBLIBRARIES.items():
        link = tmp / f"{sublib}.bam"
        if link.exists() or link.is_symlink():
            continue
        link.symlink_to(cfg["bam"])


def main() -> None:
    for subdir in ["barcodes", "logs", "looms", "ref", "figures", "tables", "tmp"]:
        (OUT / subdir).mkdir(parents=True, exist_ok=True)
    write_barcodes_and_metadata()
    write_bam_symlinks()
    gene_gtf = write_prefixed_gene_gtf()
    repeat_gtf = write_repeatmasker_gtf()
    print(f"Wrote barcode lists and metadata under {OUT}")
    print(f"Gene GTF: {gene_gtf}")
    print(f"RepeatMasker GTF: {repeat_gtf}")


if __name__ == "__main__":
    main()
