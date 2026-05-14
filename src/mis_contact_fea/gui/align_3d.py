"""3D alignment GUI — PLANNED, not yet implemented.

When implemented, this will be a PyVista- or trame-based interactive
viewer that lets you position two 3D bodies relative to each other
(translation, rotation, mirror) and export the resulting `SimulationConfig`
for the 3D solver.

For 2D alignment today, use `mis_contact_fea.gui.align_2d` (Streamlit).
"""
from __future__ import annotations


def main():
    """Placeholder entry point."""
    raise NotImplementedError(
        "3D alignment GUI is not yet implemented. Use the 2D Streamlit "
        "GUI at mis_contact_fea.gui.align_2d (or `./scripts/launch_gui.sh`) "
        "for now."
    )


if __name__ == "__main__":
    raise SystemExit(
        "mis_contact_fea.gui.align_3d is a roadmap placeholder. "
        "Use mis_contact_fea.gui.align_2d for the current 2D alignment GUI."
    )
