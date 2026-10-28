"""Loader for the Kaggle OpenVaccine (stanford-covid-vaccine) train.json file."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from aptalens.data.base import SequenceDataset
from aptalens.data.registry import register_dataset

TARGET_COLUMNS = ("reactivity", "deg_Mg_pH10", "deg_pH10", "deg_Mg_50C", "deg_50C")
SCORED_TARGET_COLUMNS = ("reactivity", "deg_Mg_pH10", "deg_Mg_50C")
ERROR_COLUMNS = {
    "reactivity": "reactivity_error",
    "deg_Mg_pH10": "deg_error_Mg_pH10",
    "deg_pH10": "deg_error_pH10",
    "deg_Mg_50C": "deg_error_Mg_50C",
    "deg_50C": "deg_error_50C",
}
REQUIRED_COLUMNS = (
    "id",
    "sequence",
    "structure",
    "predicted_loop_type",
    "signal_to_noise",
    "SN_filter",
    "seq_length",
    "seq_scored",
    *TARGET_COLUMNS,
)


def _read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _pad(values: list[float], record: dict) -> np.ndarray:
    """Pad a target list (length seq_scored) with NaN up to seq_length."""
    if len(values) != record["seq_scored"]:
        raise ValueError(f"{record['id']}: target length {len(values)} != seq_scored")
    out = np.full(record["seq_length"], np.nan)
    out[: record["seq_scored"]] = values
    return out


@register_dataset("openvaccine")
def load_openvaccine(path: str | Path, include_errors: bool = True) -> SequenceDataset:
    """Load train.json (newline-delimited JSON) into a SequenceDataset.

    The Kaggle `index` column is ignored. Targets exist only for positions
    0..seq_scored-1; later positions are NaN and flagged scored=False.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    records = _read_jsonl(path)
    if not records:
        raise ValueError(f"{path} contains no records")
    missing = [c for c in REQUIRED_COLUMNS if c not in records[0]]
    if missing:
        raise ValueError(f"train file is missing columns: {missing}")

    lengths = np.array([r["seq_length"] for r in records], dtype=int)

    # Per-molecule table (sequence, structure, loop type, metadata)
    scored = [np.concatenate([r[c] for c in SCORED_TARGET_COLUMNS]) for r in records]
    molecules = pd.DataFrame(
        {
            "id": [r["id"] for r in records],
            "sequence": [r["sequence"] for r in records],
            "structure": [r["structure"] for r in records],
            "loop_type": [r["predicted_loop_type"] for r in records],
            "seq_length": lengths,
            "seq_scored": [r["seq_scored"] for r in records],
            "signal_to_noise": [r["signal_to_noise"] for r in records],
            "sn_filter": [int(r["SN_filter"]) for r in records],
            "min_target": [float(np.min(v)) for v in scored],
            "max_abs_target": [float(np.max(np.abs(v))) for v in scored],
        }
    )

    # Per-position table, grouped by molecule in the same order
    positions = pd.DataFrame(
        {
            "id": np.repeat(molecules["id"].to_numpy(), lengths),
            "pos": np.concatenate([np.arange(n) for n in lengths]),
            "base": np.array(list("".join(molecules["sequence"]))),
            "struct": np.array(list("".join(molecules["structure"]))),
            "loop": np.array(list("".join(molecules["loop_type"]))),
        }
    )
    n_scored = np.repeat(molecules["seq_scored"].to_numpy(), lengths)
    positions["scored"] = positions["pos"].to_numpy() < n_scored
    for col in TARGET_COLUMNS:
        positions[col] = np.concatenate([_pad(r[col], r) for r in records])
    if include_errors:
        for col, err_col in ERROR_COLUMNS.items():
            if err_col in records[0]:
                positions[f"{col}_error"] = np.concatenate([_pad(r[err_col], r) for r in records])

    dataset = SequenceDataset(
        name="openvaccine",
        molecules=molecules,
        positions=positions,
        target_columns=TARGET_COLUMNS,
        scored_target_columns=SCORED_TARGET_COLUMNS,
    )
    dataset.validate()
    return dataset