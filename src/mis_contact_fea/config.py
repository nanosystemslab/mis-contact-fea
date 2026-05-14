"""Configuration schema for one mis-contact-fea simulation run.

A single YAML config fully specifies a simulation:
    geometry      pitch, pillar height, plate thickness, initial gap
    profile       which lobe shape to use (builtin / csv / svg / step …)
    mesh          characteristic element size
    solver        Newton tolerances, contact mode, direction, total disp
    output        where to write mesh + results + animation

Schema example (also see examples/01_seven_builtins/sphere_push.yaml):

    geometry:
      pitch_um: 198.0
      Hp_um: 220.0
      plate_um: 30.0
      initial_gap_um: -55.0
    profile:
      kind: builtin
      name: sphere
    mesh:
      characteristic_size_um: 4.0
    solver:
      steps: 425
      disp_um: 85.0
      direction: down
      contact_mode: raytracing
      fric: 0.0
      newton_rtol: 0.2
      newton_atol: 1000.0
      gamma_scale: 0.2
    output:
      out_dir: runs_all_shapes/sphere
"""
from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal, Optional

import yaml
from pydantic import BaseModel, ConfigDict, Field


class GeometryConfig(BaseModel):
    """Unit-cell dimensions and starting offset of the top body."""

    model_config = ConfigDict(extra="forbid")

    pitch_um: Annotated[float, Field(gt=0)] = 198.0
    Hp_um: Annotated[float, Field(gt=0)] = 220.0
    plate_um: Annotated[float, Field(gt=0)] = 30.0
    # Negative = bodies overlap at start (= engaged); positive = clearance.
    initial_gap_um: float


class ProfileConfig(BaseModel):
    """How to obtain the right-half axisymmetric polyline.

    `kind=builtin` uses one of the registered shapes in
    `mis_contact_fea.profiles.SHAPES`. The other kinds read a polyline
    from `path` — see the loaders in `mis_contact_fea.profiles.from_*`.

    `units` declares what length unit the values in `path` are in.
    Internally everything is stored in µm, so e.g. `units=mm` triggers
    a × 1000 conversion at load time."""

    model_config = ConfigDict(extra="forbid")

    kind: Literal["builtin", "csv", "svg", "step"] = "builtin"
    name: Optional[str] = None  # for kind=builtin
    path: Optional[Path] = None  # for kind in {csv, svg, step}
    units: Literal["um", "mm", "cm", "m"] = "um"
    # Optional: stitch the standard pillar to the loaded lobe outline.
    # Only meaningful for kind in {csv, svg}; STEP usually already
    # includes its own pillar.
    prepend_pillar: bool = False


class MeshConfig(BaseModel):
    """Gmsh meshing parameters."""

    model_config = ConfigDict(extra="forbid")

    characteristic_size_um: Annotated[float, Field(gt=0)] = 4.0


class SolverConfig(BaseModel):
    """dolfinx_contact step-control + tolerances."""

    model_config = ConfigDict(extra="forbid")

    steps: Annotated[int, Field(gt=0)]
    disp_um: Annotated[float, Field(gt=0)]
    direction: Literal["down", "up"] = "down"  # down = push, up = retention
    contact_mode: Literal["raytracing", "closest"] = "raytracing"
    fric: Annotated[float, Field(ge=0)] = 0.0
    newton_rtol: Annotated[float, Field(gt=0)] = 0.2
    newton_atol: Annotated[float, Field(gt=0)] = 1000.0
    gamma_scale: Annotated[float, Field(gt=0)] = 0.2
    ksp_type: str = "preonly"
    pc_type: str = "lu"


class OutputConfig(BaseModel):
    """Where the mesher and solver write their artifacts."""

    model_config = ConfigDict(extra="forbid")

    out_dir: Path


class SimulationConfig(BaseModel):
    """Top-level config — one of these fully specifies a single sim run.

    The `profile` field is the *bottom* body's lobe (and the top body's
    by default, for a symmetric pair). Set `top_profile` to use a
    different shape on the top body — asymmetric pair, e.g. sphere
    bottom vs. cap top. When `top_profile` is omitted (the common case)
    behavior is unchanged from earlier versions.
    """

    model_config = ConfigDict(extra="forbid")

    geometry: GeometryConfig
    profile: ProfileConfig
    top_profile: Optional[ProfileConfig] = None
    mesh: MeshConfig = MeshConfig()
    solver: SolverConfig
    output: OutputConfig

    @property
    def effective_top_profile(self) -> ProfileConfig:
        """Profile used for the top body (falls back to `profile`)."""
        return self.top_profile if self.top_profile is not None else self.profile

    @classmethod
    def from_yaml(cls, path: Path | str) -> "SimulationConfig":
        """Load and validate a config from a YAML file."""
        path = Path(path)
        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return cls.model_validate(data)

    def to_yaml(self, path: Path | str) -> None:
        """Round-trip the config back to YAML (useful for the GUI)."""
        path = Path(path)
        data = self.model_dump(mode="json", exclude_none=True)
        with path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(data, f, sort_keys=False, default_flow_style=False)

    # Convenience accessors so existing code doesn't have to traverse nested.
    @property
    def mesh_dir(self) -> Path:
        """Where mesher writes mesh.msh + mesh_tags.json."""
        return self.output.out_dir / "mesh"

    @property
    def results_dir(self) -> Path:
        """Where solver writes force_displacement.csv + contact_results.xdmf."""
        return self.output.out_dir / "results"

    @property
    def mesh_preview_path(self) -> Path:
        return self.output.out_dir / "mesh_preview.png"

    @property
    def animation_path(self) -> Path:
        return self.output.out_dir / "contact_animation_pv.gif"
