"""mis-contact-fea — Mechanical Interlocking Structure Contact FEA.

A 2D plane-strain finite-element pipeline for characterizing the
insertion vs retention force asymmetry of axisymmetric nanostructure
shapes arranged as a mechanical interlock.

Subpackages:
    profiles    builtin shape generators + loaders (CSV / SVG / STEP …)
    mesh        Gmsh → dolfinx meshers
    solver      dolfinx + dolfinx_contact contact-FEA driver
    postproc    smoothing, plotting, animation, analysis
    gui         interactive alignment tools (Streamlit)
"""
from __future__ import annotations

__version__ = "0.1.0"
