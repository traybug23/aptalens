"""Dataset container and validation for per-nucleotide sequence data."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd

REQUIRED_MOLECULE_COLUMNS = ("id", "sequence", "structure", "loop_type", "seq_length", "seq_scored")
REQUIRED_POSITION_COLUMNS = ("id", "pos", "base", "struct", "paired", "loop_type", "scored")


@dataclass(frozen=True)
class SequenceDataset:
    """A set of molecules, stored two ways.

    molecules: one row per molecule (id, sequence, structure, loop_type, lengths, quality columns).
    positions: one row per nucleotide position, every row keeping its molecule `id`.
            Target columns are NaN outside the scored region.

    Splitting must always be done on molecule ids, never on rows of `positions`.
    """

    name: str
    molecules: pd.DataFrame
    positions: pd.DataFrame
    scored_targets: tuple[str, ...]
    unscored_targets: tuple[str, ...] = ()

    @property
    def targets(self) -> tuple[str, ...]:
        return self.scored_targets + self.unscored_targets

    @property
    def ids(self) -> np.ndarray:
        return self.molecules["id"].to_numpy()

    def __len__(self) -> int:
        return len(self.molecules)

    def scored_positions(self) -> pd.DataFrame:
        """Rows of `positions` that carry ground-truth targets."""
        return self.positions[self.positions["scored"]]

    def subset(self, ids: Iterable[str]) -> "SequenceDataset":
        """Return a dataset restricted to whole molecules (never partial molecules)."""
        wanted = pd.Index(list(ids))
        unknown = wanted.difference(pd.Index(self.molecules["id"]))
        if len(unknown) > 0:
            raise KeyError(f"unknown ids: {list(unknown[:5])}")
        mols = self.molecules[self.molecules["id"].isin(wanted)].reset_index(drop=True)
        pos = self.positions[self.positions["id"].isin(wanted)].reset_index(drop=True)
        return SequenceDataset(self.name, mols, pos, self.scored_targets, self.unscored_targets)


def validate(ds: SequenceDataset) -> None:
    """Raise ValueError if the dataset breaks the invariants the pipeline relies on."""
    mol, pos = ds.molecules, ds.positions

    for col in REQUIRED_MOLECULE_COLUMNS:
        if col not in mol.columns:
            raise ValueError(f"molecules is missing column '{col}'")
    for col in REQUIRED_POSITION_COLUMNS + ds.targets:
        if col not in pos.columns:
            raise ValueError(f"positions is missing column '{col}'")

    # One sequence per id, no missing ids.
    if mol["id"].isna().any():
        raise ValueError("missing ids in molecules")
    if mol["id"].duplicated().any():
        raise ValueError("duplicate ids in molecules (expected one row per molecule)")
    if pos["id"].isna().any():
        raise ValueError("missing ids in positions")
    if not set(pos["id"]).issubset(set(mol["id"])):
        raise ValueError("positions contain ids that are not in molecules")

    # Lengths agree everywhere.
    seq_len = mol["sequence"].str.len()
    if not (seq_len == mol["seq_length"]).all():
        raise ValueError("sequence length differs from seq_length for some molecules")

    per_id = pos.groupby("id")
    counts = per_id.size().reindex(mol["id"])
    if not (counts.to_numpy() == mol["seq_length"].to_numpy()).all():
        raise ValueError("number of position rows differs from seq_length for some molecules")
    if not (per_id["pos"].max().reindex(mol["id"]).to_numpy() == (mol["seq_length"] - 1).to_numpy()).all():
        raise ValueError("position index does not run 0..seq_length-1 for some molecules")
    if pos.duplicated(["id", "pos"]).any():
        raise ValueError("duplicate (id, pos) rows")

    # Scored region has the documented size, and targets are NaN outside it.
    scored_counts = per_id["scored"].sum().reindex(mol["id"])
    if not (scored_counts.to_numpy() == mol["seq_scored"].to_numpy()).all():
        raise ValueError("number of scored positions differs from seq_scored for some molecules")
    outside = pos.loc[~pos["scored"], list(ds.targets)]
    if outside.notna().any().any():
        raise ValueError("targets are present outside the scored region")