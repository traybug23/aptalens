"""Name -> loader registry so callers can write load_dataset("openvaccine", path=...)."""

from __future__ import annotations

from collections.abc import Callable

from aptalens.data.base import SequenceDataset

_REGISTRY: dict[str, Callable[..., SequenceDataset]] = {}


def register_dataset(name: str):
    def decorator(func: Callable[..., SequenceDataset]):
        if name in _REGISTRY:
            raise ValueError(f"dataset '{name}' is already registered")
        _REGISTRY[name] = func
        return func

    return decorator


def load_dataset(name: str, **kwargs) -> SequenceDataset:
    if name not in _REGISTRY:
        raise KeyError(f"unknown dataset '{name}'. Available: {list_datasets()}")
    return _REGISTRY[name](**kwargs)


def list_datasets() -> list[str]:
    return sorted(_REGISTRY)