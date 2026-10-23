import importlib


def test_import():
    for module_name in (
        "aptalens.data",
        "aptalens.splits",
        "aptalens.audit",
        "aptalens.features",
        "aptalens.models",
        "aptalens.eval",
    ):
        importlib.import_module(module_name)
