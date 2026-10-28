import json
from pathlib import Path

import numpy as np
import pytest

from aptalens.data import assert_disjoint_ids, list_datasets, load_dataset
from aptalens.data.openvaccine import ERROR_COLUMNS, TARGET_COLUMNS

N_TOTAL, N_SCORED = 10, 6


def _record(i, structure="(((...)))."):
    rec = {
        "index": i,
        "id": f"id_{i}",
        "sequence": "GGGAAACCCA",
        "structure": structure,
        "predicted_loop_type": "SSSHHHSSSE",
        "signal_to_noise": 2.0 + i,
        "SN_filter": i % 2,
        "seq_length": N_TOTAL,
        "seq_scored": N_SCORED,
    }
    for k, col in enumerate(TARGET_COLUMNS):
        rec[col] = [0.1 * (i + k)] * N_SCORED
    for err in ERROR_COLUMNS.values():
        rec[err] = [0.05] * N_SCORED
    return rec


@pytest.fixture
def train_path(tmp_path):
    path = tmp_path / "train.json"
    records = [_record(i) for i in range(3)]
    records[2]["deg_Mg_50C"][0] = -20.0
    path.write_text("\n".join(json.dumps(r) for r in records))
    return path


def test_registered():
    assert "openvaccine" in list_datasets()


def test_accessors(train_path):
    ds = load_dataset("openvaccine", path=train_path)
    assert len(ds) == 3
    assert ds.ids == ["id_0", "id_1", "id_2"]
    assert ds.sequences[0] == "GGGAAACCCA"
    assert ds.structures[0] == "(((...)))."
    assert ds.loop_types[0] == "SSSHHHSSSE"
    assert "sequence" not in ds.metadata.columns
    assert {"sn_filter", "signal_to_noise", "min_target"} <= set(ds.metadata.columns)
    assert len(ds.positions) == 3 * N_TOTAL
    assert ds.targets.shape == (3 * N_TOTAL, 5)


def test_targets_only_on_scored_positions(train_path):
    ds = load_dataset("openvaccine", path=train_path)
    assert ds.positions["reactivity"].notna().sum() == 3 * N_SCORED
    assert ds.positions.loc[~ds.positions["scored"], "reactivity"].isna().all()
    arr = ds.target_array()
    assert arr.shape == (3, N_SCORED, 3)
    assert ds.molecules.loc[2, "min_target"] == -20.0


def test_subset_keeps_whole_molecules(train_path):
    ds = load_dataset("openvaccine", path=train_path).subset(["id_2", "id_0"])
    assert ds.ids == ["id_0", "id_2"]
    assert len(ds.positions) == 2 * N_TOTAL
    ds.validate()
    with pytest.raises(KeyError):
        ds.subset(["nope"])


def test_validate_catches_unbalanced_structure(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(_record(0, structure="(((...))..")))
    with pytest.raises(ValueError, match="unbalanced"):
        load_dataset("openvaccine", path=path)


def test_missing_file():
    with pytest.raises(FileNotFoundError):
        load_dataset("openvaccine", path="does_not_exist.json")


def test_expected_shapes(train_path):
    ds = load_dataset("openvaccine", path=train_path)
    assert len(ds.molecules) == 3
    assert len(ds.positions) == 3 * N_TOTAL
    assert ds.target_array().shape == (3, N_SCORED, 3)
    assert ds.positions.groupby("id").size().eq(N_TOTAL).all()


def test_no_missing_ids(train_path):
    ds = load_dataset("openvaccine", path=train_path)
    assert ds.molecules["id"].notna().all()
    assert ds.positions["id"].notna().all()
    assert set(ds.positions["id"]) == set(ds.molecules["id"])
    assert ds.molecules["id"].is_unique


def test_one_sequence_per_id(train_path):
    ds = load_dataset("openvaccine", path=train_path)
    rebuilt = ds.positions.groupby("id", sort=False)["base"].agg("".join)
    assert rebuilt.loc[ds.ids].tolist() == ds.sequences


def test_validate_catches_corruption(train_path):
    ds = load_dataset("openvaccine", path=train_path)
    ds.positions.loc[0, "base"] = "X"
    with pytest.raises(ValueError, match="exactly one sequence"):
        ds.validate()
    ds = load_dataset("openvaccine", path=train_path)
    ds.positions = ds.positions.drop(index=5)
    with pytest.raises(ValueError, match="row count"):
        ds.validate()
    ds = load_dataset("openvaccine", path=train_path)
    ds.molecules.loc[0, "id"] = None
    with pytest.raises(ValueError, match="null"):
        ds.validate()


def test_subset_never_splits_a_molecule(train_path):
    ds = load_dataset("openvaccine", path=train_path).subset(["id_1"])
    assert ds.positions["id"].eq("id_1").all()
    assert len(ds.positions) == N_TOTAL


def test_assert_disjoint_ids():
    assert_disjoint_ids(train=["a", "b"], valid=["c"])
    with pytest.raises(ValueError, match="both"):
        assert_disjoint_ids(train=["a", "b"], valid=["b"])


REAL = Path(__file__).resolve().parents[1] / "data" / "raw" / "train.json"


@pytest.mark.skipif(not REAL.exists(), reason="real train.json not present")
def test_real_data():
    ds = load_dataset("openvaccine", path=REAL)
    assert len(ds) == 2400
    assert len(ds.positions) == 2400 * 107
    assert ds.positions["reactivity"].notna().sum() == 2400 * 68
    assert np.isfinite(ds.target_array()).all()
    assert ds.molecules["id"].is_unique
    assert ds.positions.groupby("id").size().eq(107).all()