"""Core data container shared by every dataset loader."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
import pandas as pd

SEQUENCE_COLUMNS = ("sequence", "structure", "loop_type")


def assert_disjoint_ids(**id_groups: Iterable[str]) -> None:
    """Raise ValueError if any molecule id appears in more than one split.

    Usage: assert_disjoint_ids(train=train_ids, valid=valid_ids, test=test_ids)
    Splitting must always be by molecule id, never by row, so that positions
    of the same RNA never end up on both sides of a split.
    """
    seen: dict[str, str] = {}
    for group, ids in id_groups.items():
        for molecule_id in ids:
            if molecule_id in seen and seen[molecule_id] != group:
                raise ValueError(f"id {molecule_id} is in both '{seen[molecule_id]}' and '{group}'")
            seen[molecule_id] = group


@dataclass
class SequenceDataset:
    """A set of nucleic-acid molecules in two linked tables.

    molecules: one row per molecule (id, sequence, structure, loop_type + metadata).
    positions: one row per position, grouped by molecule in the same order as
        `molecules`. Target columns are NaN where a position is not scored.
    """

    name: str
    molecules: pd.DataFrame
    positions: pd.DataFrame
    target_columns: tuple[str, ...]
    scored_target_columns: tuple[str, ...]

    def __len__(self) -> int:
        return len(self.molecules)

    # ---- the six things a model or splitter usually asks for ----
    @property
    def ids(self) -> list[str]:
        return self.molecules["id"].tolist()

    @property
    def sequences(self) -> list[str]:
        return self.molecules["sequence"].tolist()

    @property
    def structures(self) -> list[str]:
        return self.molecules["structure"].tolist()

    @property
    def loop_types(self) -> list[str]:
        return self.molecules["loop_type"].tolist()

    @property
    def targets(self) -> pd.DataFrame:
        """All target columns for every position (NaN where unscored)."""
        return self.positions.loc[:, list(self.target_columns)]

    @property
    def metadata(self) -> pd.DataFrame:
        """Per-molecule table without the long sequence strings."""
        return self.molecules.drop(columns=list(SEQUENCE_COLUMNS))

    # ---- helpers ----
    def scored_positions(self) -> pd.DataFrame:
        return self.positions[self.positions["scored"]]

    def target_array(self, columns: Iterable[str] | None = None) -> np.ndarray:
        """Scored targets as an array of shape (n_molecules, n_scored, n_targets)."""
        cols = list(columns) if columns is not None else list(self.scored_target_columns)
        counts = self.molecules["seq_scored"].unique()
        if len(counts) != 1:
            raise ValueError("target_array needs the same seq_scored for every molecule")
        values = self.scored_positions()[cols].to_numpy(dtype=float)
        return values.reshape(len(self), int(counts[0]), len(cols))

    def subset(self, ids: Iterable[str]) -> SequenceDataset:
        """Return a dataset restricted to whole molecules (never partial molecules)."""
        wanted = pd.Index(list(ids))
        unknown = wanted.difference(self.molecules["id"])
        if len(unknown) > 0:
            raise KeyError(f"{len(unknown)} unknown ids, e.g. {list(unknown[:3])}")
        return SequenceDataset(
            name=self.name,
            molecules=self.molecules[self.molecules["id"].isin(wanted)].reset_index(drop=True),
            positions=self.positions[self.positions["id"].isin(wanted)].reset_index(drop=True),
            target_columns=self.target_columns,
            scored_target_columns=self.scored_target_columns,
        )

    def validate(self) -> None:
        """Raise ValueError if the tables are inconsistent."""
        mol, pos = self.molecules, self.positions
        if mol["id"].isna().any() or pos["id"].isna().any():
            raise ValueError("missing (null) ids")
        if not mol["id"].is_unique:
            raise ValueError("molecule ids are not unique")
        for col in ("structure", "loop_type"):
            bad = mol[mol[col].str.len() != mol["seq_length"]]
            if len(bad):
                raise ValueError(f"{col} length != seq_length for {bad['id'].iloc[0]}")
        bad = mol[mol["sequence"].str.len() != mol["seq_length"]]
        if len(bad):
            raise ValueError(f"sequence length != seq_length for {bad['id'].iloc[0]}")
        unbalanced = mol[mol["structure"].str.count(r"\(") != mol["structure"].str.count(r"\)")]
        if len(unbalanced):
            raise ValueError(f"unbalanced dot-bracket for {unbalanced['id'].iloc[0]}")
        if len(pos) != int(mol["seq_length"].sum()):
            raise ValueError("positions row count != sum of seq_length")
        n_pos = mol["seq_length"].to_numpy()
        if not np.array_equal(pos["id"].to_numpy(), np.repeat(mol["id"].to_numpy(), n_pos)):
            raise ValueError("positions are not grouped in molecule order")
        n_scored = np.repeat(mol["seq_scored"].to_numpy(), n_pos)
        expected_scored = pos["pos"].to_numpy() < n_scored
        if not np.array_equal(pos["scored"].to_numpy(), expected_scored):
            raise ValueError("'scored' flag does not match seq_scored")
        has_value = pos[list(self.target_columns)].notna().to_numpy()
        if has_value[~expected_scored].any():
            raise ValueError("targets present at unscored positions")
        if not has_value[expected_scored].all():
            raise ValueError("missing targets at scored positions")
        if set(pos["id"]) != set(mol["id"]):
            raise ValueError("ids in positions and molecules do not match")
        pairs = (("sequence", "base"), ("structure", "struct"), ("loop_type", "loop"))
        for mol_col, pos_col in pairs:
            rebuilt = pos.groupby("id", sort=False)[pos_col].agg("".join)
            expected = mol.set_index("id")[mol_col]
            if rebuilt.loc[expected.index].tolist() != expected.tolist():
                raise ValueError(f"positions do not rebuild exactly one {mol_col} per id")