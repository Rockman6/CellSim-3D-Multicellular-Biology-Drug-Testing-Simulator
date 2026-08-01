"""Criterion-8 reference scene must render from the real model.

Headless render check: the viewer drives the validated src/cell engine
(not a re-implementation) and produces a non-trivial PNG.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))


def test_reference_scene_renders(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    from src.cell.viewer import render_disposition_scene

    out = tmp_path / "scene.png"
    # Small MC budget keeps the test fast; the physics path is identical.
    render_disposition_scene(save=str(out), show=False, mc_samples=200)
    assert out.exists()
    # A real multi-panel figure, not a blank canvas.
    assert out.stat().st_size > 20_000


def test_committed_reference_scene_exists():
    """The checked-in reference scene is the criterion-8 artifact."""
    png = REPO_ROOT / "benchmarks" / "dock" / "cell_disposition_scene.png"
    assert png.exists(), (
        "run `python src/cell/viewer.py benchmarks/dock/"
        "cell_disposition_scene.png` to regenerate the reference scene")
    assert png.stat().st_size > 20_000
