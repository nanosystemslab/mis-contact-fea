"""Guard against drift between the package source and the inlined copies
in `src/mis_contact_fea/gui/align_marimo.py`.

The marimo notebook is self-contained (it has to be, to run in Pyodide
WASM where the package isn't installable). That means the seven shape
generators and the polyline-layout helpers are duplicated. This test
extracts them from the notebook and compares polylines against the
package canonical version.

If you update either side, run the tests — if they fail, sync the other
side and re-run.
"""
from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pytest

from mis_contact_fea.mesh.layout import (
    bottom_body_vertices,
    top_body_vertices,
)
from mis_contact_fea.profiles import SHAPES


NOTEBOOK = Path(__file__).resolve().parent.parent / "src" / "mis_contact_fea" / "gui" / "align_marimo.py"


def _extract_cell_globals() -> dict:
    """Execute the marimo notebook's `app.cell` function bodies in
    isolation, harvesting the SHAPES dict and the layout helpers."""
    source = NOTEBOOK.read_text()
    tree = ast.parse(source)

    namespace: dict = {"np": np}
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "_":
            for stmt in node.body:
                if isinstance(stmt, ast.Return):
                    continue
                code = compile(ast.Module(body=[stmt], type_ignores=[]), str(NOTEBOOK), "exec")
                try:
                    exec(code, namespace)
                except Exception:  # noqa: BLE001
                    pass
    return namespace


@pytest.fixture(scope="module")
def notebook_namespace() -> dict:
    return _extract_cell_globals()


def test_notebook_defines_shapes(notebook_namespace: dict) -> None:
    """Sanity: the notebook should expose a SHAPES dict via cell return."""
    assert "SHAPES" in notebook_namespace, (
        "couldn't locate SHAPES dict in the marimo notebook — has the "
        "cell structure changed?"
    )
    nb_shapes = notebook_namespace["SHAPES"]
    assert set(nb_shapes) == set(SHAPES), (
        "notebook SHAPES keys differ from package SHAPES keys"
    )


def test_notebook_shapes_match_package(notebook_namespace: dict) -> None:
    nb_shapes = notebook_namespace["SHAPES"]
    for name in SHAPES:
        pkg_x, pkg_z = SHAPES[name](Hp=220.0)
        nb_x, nb_z = nb_shapes[name](Hp=220.0)
        np.testing.assert_array_equal(
            nb_x, pkg_x, err_msg=f"x mismatch for shape {name!r} — sync align_marimo.py"
        )
        np.testing.assert_array_equal(
            nb_z, pkg_z, err_msg=f"z mismatch for shape {name!r} — sync align_marimo.py"
        )


def test_notebook_layout_helpers_match_package(notebook_namespace: dict) -> None:
    assert "bottom_body_vertices" in notebook_namespace
    assert "top_body_vertices" in notebook_namespace
    nb_bb = notebook_namespace["bottom_body_vertices"]
    nb_tb = notebook_namespace["top_body_vertices"]

    half_x, half_z = SHAPES["sphere"](Hp=220.0)
    for fn_name, nb_fn, pkg_fn, args in [
        ("bottom_body_vertices", nb_bb, bottom_body_vertices, (half_x, half_z, 198.0, 30.0)),
        ("top_body_vertices",    nb_tb, top_body_vertices,    (half_x, half_z, 198.0, 30.0, -55.0)),
    ]:
        nb_x, nb_z, nb_tags = nb_fn(*args)
        pkg_x, pkg_z, pkg_tags = pkg_fn(*args)
        np.testing.assert_array_equal(nb_x, pkg_x, err_msg=f"x mismatch in {fn_name}")
        np.testing.assert_array_equal(nb_z, pkg_z, err_msg=f"z mismatch in {fn_name}")
        assert nb_tags == pkg_tags, f"tag mismatch in {fn_name}"
