"""Reading exported matrices gives the same ModuleData the store route gives.

The export is what makes a deposited copy of the inputs enough to re-run an analysis, so
the two routes have to agree, including the DE mask at a threshold other than the default.
"""
import numpy as np

from mmc.data import module_data
from mmc.eval.holdout import ModuleData

GENES = ["A", "B", "C", "D"]
PERTS = ["A", "B"]
OBS = np.array([[0.0, 1.5, -0.4, 0.2],
                [2.0, 0.0, 0.1, -3.0]])
Q = np.array([[1.00, 0.01, 0.40, 0.08],
              [0.02, 1.00, 0.90, 0.00]])


def _write(tmp_path, module="mod", condition="Stim8hr"):
    np.savez_compressed(tmp_path / f"{module}__{condition}.npz",
                        genes=np.array(GENES, dtype=str),
                        perts=np.array(PERTS, dtype=str),
                        effect_size=OBS, fdr=Q)


def test_matrices_route_reproduces_the_store_route(tmp_path, monkeypatch):
    _write(tmp_path)
    monkeypatch.setenv(module_data.MATRICES_ENV, str(tmp_path))

    # the store route, standing in for a store that is not present in the test environment
    monkeypatch.setattr(module_data, "store_matrices",
                        lambda m, c: (GENES, PERTS, OBS, Q))
    from_store = ModuleData(GENES, PERTS, OBS, Q < module_data.DEFAULT_FDR, None)

    got = module_data.build_module_data("mod", "Stim8hr")
    assert got.genes == from_store.genes
    assert got.perts == from_store.perts
    assert np.array_equal(got.observed, from_store.observed)
    assert np.array_equal(got.de_mask, from_store.de_mask)


def test_the_threshold_is_applied_on_read_not_baked_in(tmp_path, monkeypatch):
    _write(tmp_path)
    monkeypatch.setenv(module_data.MATRICES_ENV, str(tmp_path))

    loose = module_data.build_module_data("mod", "Stim8hr", fdr=0.5)
    tight = module_data.build_module_data("mod", "Stim8hr", fdr=0.01)
    assert loose.de_mask.sum() > tight.de_mask.sum()
    assert np.array_equal(loose.de_mask, Q < 0.5)
    assert np.array_equal(tight.de_mask, Q < 0.01)


def test_a_missing_module_falls_back_to_the_store(tmp_path, monkeypatch):
    monkeypatch.setenv(module_data.MATRICES_ENV, str(tmp_path))
    called = []

    def fake(module, condition):
        called.append((module, condition))
        return GENES, PERTS, OBS, Q

    monkeypatch.setattr(module_data, "store_matrices", fake)
    module_data.build_module_data("absent", "Stim8hr")
    assert called == [("absent", "Stim8hr")]
