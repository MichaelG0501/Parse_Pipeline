"""Shared Parse_Pipeline configuration for Python workflows."""

import os
from pathlib import Path


def parse_project_root() -> Path:
    env_root = os.environ.get("AUTO_PARSE_ROOT_DIR")
    if env_root:
        return Path(env_root)
    candidates = [
        Path("/rds/general/project/spatialtranscriptomics/ephemeral/Parse_Pipeline"),
        Path("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not locate Parse_Pipeline root. Set AUTO_PARSE_ROOT_DIR.")


PROJECT_ROOT = parse_project_root()
PARSE_OUTS = PROJECT_ROOT / "parse_outs"
PARSE_SAMPLES = ["T0", "T1", "T2", "T4", "R4", "eR4"]
ALL_SAMPLES = PARSE_SAMPLES + ["PDO", "SUR1090_Untreated", "SUR1090_Treated"]
VELOCITY_OUT = PARSE_OUTS / "Auto_velocity"
VELOCITY_GTF = Path(
    "/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
)
OUTPUT_TIERS = {
    "intermediate": PARSE_OUTS / "intermediate",
    "tables": PARSE_OUTS / "tables",
    "figures": PARSE_OUTS / "figures",
    "logs": PARSE_OUTS / "logs",
    "reports": PARSE_OUTS / "reports",
}
