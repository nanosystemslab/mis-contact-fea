#!/usr/bin/env python3
"""Quasi-static 2D plane-strain contact for the offset-by-P/2 unit cell.

Modeled after src/run_contact_sim.py but reduced to gdim=2 with:
  - plane-strain Lame parameters
  - Dirichlet u_x=0 on sym_left/sym_right (mirror-symmetric layout =
    exactly equivalent to periodic BC for offset = P/2)
  - dolfinx_contact Nitsche between bottom_contact and top_contact

Run serially (single rank):
    python3 src/run_2d_slice_sim.py --mesh-dir mesh_out_2d ...

Run in parallel with MPI:
    mpirun -n 4 python3 src/run_2d_slice_sim.py --mesh-dir mesh_out_2d ...

Units are µm / MPa / µN (lengths in µm, E in MPa = µN/µm², reaction in
µN per unit out-of-plane depth; multiply by physical depth for total
force on a 3D realization).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI
from petsc4py.PETSc import InsertMode, ScatterMode  # type: ignore

from dolfinx import default_scalar_type, fem, io
from dolfinx.io import gmshio
from dolfinx.fem.petsc import (
    apply_lifting,
    assemble_matrix,
    assemble_vector,
    create_vector,
    set_bc,
)
from dolfinx.graph import adjacencylist
from dolfinx_contact.cpp import ContactMode
from dolfinx_contact.general_contact.contact_problem import ContactProblem, FrictionLaw
from dolfinx_contact.helpers import epsilon, lame_parameters, rigid_motions_nullspace_subdomains, sigma_func
from dolfinx_contact.newton_solver import NewtonSolver


def _load_meta(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="2D plane-strain displacement-controlled contact simulation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", type=Path, default=None,
                        help="YAML simulation config (provides defaults for the flags below).")
    parser.add_argument("--mesh-dir", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--steps", type=int, default=None)
    parser.add_argument("--disp-um", type=float, default=None,
                        help="Total top-plate downward displacement (µm). "
                             "If unset, defaults to initial_gap + 3*z_apex so "
                             "the top cap fully passes the bottom cap.")
    parser.add_argument("--E-mpa", type=float, default=200000.0,
                        help="Young's modulus in MPa.")
    parser.add_argument("--nu", type=float, default=0.3)
    parser.add_argument("--fric", type=float, default=None)
    parser.add_argument("--q-degree", type=int, default=5)
    parser.add_argument("--gamma-scale", type=float, default=None)
    parser.add_argument("--no-update-contact", action="store_true")
    parser.add_argument("--newton-max-it", type=int, default=40)
    parser.add_argument("--newton-rtol", type=float, default=None)
    parser.add_argument("--newton-atol", type=float, default=None)
    parser.add_argument("--relax", type=float, default=1.0)
    parser.add_argument("--ksp-rtol", type=float, default=1e-6)
    parser.add_argument("--ksp-atol", type=float, default=1e-10)
    parser.add_argument("--ksp-max-it", type=int, default=2000)
    parser.add_argument("--ksp-type", type=str, default=None)
    parser.add_argument("--pc-type", type=str, default=None)
    parser.add_argument("--contact-mode", choices=["closest", "raytracing"],
                        default=None,
                        help="dolfinx_contact pair-detection mode. raytracing "
                             "shoots normals to find genuine surface contacts; "
                             "closest can produce spurious cross-cell pairs.")
    parser.add_argument("--direction", choices=["down", "up"], default=None,
                        help="Direction the top plate moves. down = push-in "
                             "(default), up = retention/pull-out. Total disp "
                             "magnitude controlled by --disp-um.")
    args = parser.parse_args()

    # If --config is given, fill any unset flag from the config; explicit
    # flags always win. Otherwise apply the legacy defaults.
    if args.config is not None:
        from mis_contact_fea.config import SimulationConfig

        cfg = SimulationConfig.from_yaml(args.config)
        if args.mesh_dir is None:
            args.mesh_dir = cfg.mesh_dir
        if args.out_dir is None:
            args.out_dir = cfg.results_dir
        if args.steps is None:
            args.steps = cfg.solver.steps
        if args.disp_um is None:
            args.disp_um = cfg.solver.disp_um
        if args.fric is None:
            args.fric = cfg.solver.fric
        if args.gamma_scale is None:
            args.gamma_scale = cfg.solver.gamma_scale
        if args.newton_rtol is None:
            args.newton_rtol = cfg.solver.newton_rtol
        if args.newton_atol is None:
            args.newton_atol = cfg.solver.newton_atol
        if args.ksp_type is None:
            args.ksp_type = cfg.solver.ksp_type
        if args.pc_type is None:
            args.pc_type = cfg.solver.pc_type
        if args.contact_mode is None:
            args.contact_mode = cfg.solver.contact_mode
        if args.direction is None:
            args.direction = cfg.solver.direction

    # Legacy defaults for the no-config path (preserves old CLI behavior).
    if args.mesh_dir is None:
        args.mesh_dir = Path("mesh_out_2d")
    if args.out_dir is None:
        args.out_dir = Path("results_2d")
    if args.steps is None:
        args.steps = 20
    if args.fric is None:
        args.fric = 0.3
    if args.gamma_scale is None:
        args.gamma_scale = 10.0
    if args.newton_rtol is None:
        args.newton_rtol = 1e-7
    if args.newton_atol is None:
        args.newton_atol = 1e-7
    if args.ksp_type is None:
        args.ksp_type = "fgmres"
    if args.pc_type is None:
        args.pc_type = "gamg"
    if args.contact_mode is None:
        args.contact_mode = "raytracing"
    if args.direction is None:
        args.direction = "down"

    if args.E_mpa <= 0:
        raise ValueError(f"E must be > 0, got {args.E_mpa}")
    if not (0 < args.nu < 0.5):
        raise ValueError(f"nu must be in (0, 0.5), got {args.nu}")
    if args.fric < 0:
        raise ValueError(f"friction must be >= 0, got {args.fric}")

    mesh_dir = args.mesh_dir.resolve()
    msh_path = mesh_dir / "mesh.msh"
    tags_path = mesh_dir / "mesh_tags.json"
    if not msh_path.exists():
        raise FileNotFoundError(f"Mesh file not found: {msh_path}")
    if not tags_path.exists():
        raise FileNotFoundError(f"Tags file not found: {tags_path}")

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    meta = _load_meta(tags_path)
    tags = meta["tags"]

    comm = MPI.COMM_WORLD

    if args.disp_um is None:
        z_apex = float(meta.get("z_apex_um", 50.0))
        initial_gap = float(meta.get("initial_gap_um", 10.0))
        disp_um = initial_gap + 3.0 * z_apex
        if comm.rank == 0:
            print(f"Auto displacement: {disp_um:.2f} µm "
                  f"(initial_gap {initial_gap:g} + 3*z_apex {3*z_apex:g})")
    else:
        disp_um = args.disp_um

    if comm.rank == 0:
        print(f"Loading mesh {msh_path} on {comm.size} rank(s)")

    mesh, cell_tags, facet_tags = gmshio.read_from_msh(str(msh_path), comm, 0, gdim=2)

    tdim = mesh.topology.dim
    fdim = tdim - 1
    mesh.topology.create_connectivity(fdim, tdim)
    mesh.topology.create_connectivity(0, tdim)

    V = fem.functionspace(mesh, ("CG", 1, (mesh.geometry.dim,)))
    V0 = fem.functionspace(mesh, ("DG", 0))
    V_stress = fem.functionspace(mesh, ("DG", 0))

    mu_func, lambda_func = lame_parameters(True)  # plane strain
    mu = mu_func(args.E_mpa, args.nu)
    lmbda = lambda_func(args.E_mpa, args.nu)
    sigma = sigma_func(mu, lmbda)

    mu_dg = fem.Function(V0)
    lmbda_dg = fem.Function(V0)
    fric_dg = fem.Function(V0)
    mu_dg.interpolate(lambda x: np.full((1, x.shape[1]), mu))
    lmbda_dg.interpolate(lambda x: np.full((1, x.shape[1]), lmbda))
    fric_dg.interpolate(lambda x: np.full((1, x.shape[1]), args.fric))

    u = fem.Function(V, name="u")
    du = fem.Function(V, name="du")
    von_mises = fem.Function(V_stress, name="von_mises")

    w = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    dx = ufl.Measure("dx", domain=mesh, subdomain_data=cell_tags)
    F = ufl.inner(sigma(u), epsilon(v)) * dx
    F = ufl.replace(F, {u: u + du})
    J = ufl.derivative(F, du, w)

    F_compiled = fem.form(F)
    J_compiled = fem.form(J)

    contact_pairs = [(0, 1), (1, 0)]
    data = np.array([tags["bottom_contact"], tags["top_contact"]], dtype=np.int32)
    offsets = np.array([0, 2], dtype=np.int32)
    surfaces = adjacencylist(data, offsets)
    if args.contact_mode == "raytracing":
        _mode = ContactMode.Raytracing
    else:
        _mode = ContactMode.ClosestPoint
    search_mode = [_mode, _mode]

    contact_problem = ContactProblem(
        [facet_tags], surfaces, contact_pairs, mesh, args.q_degree, search_mode,
    )
    contact_problem.generate_contact_data(
        FrictionLaw.Coulomb,
        V,
        {"u": u, "du": du, "mu": mu_dg, "lambda": lmbda_dg, "fric": fric_dg},
        args.E_mpa * args.gamma_scale,
        -1.0,
    )

    bot_fixed_facets = facet_tags.find(tags["bottom_fixed"])
    top_disp_facets = facet_tags.find(tags["top_disp"])
    sym_left_facets = facet_tags.find(tags["sym_left"])
    sym_right_facets = facet_tags.find(tags["sym_right"])
    bot_contact_facets = facet_tags.find(tags["bottom_contact"])
    top_contact_facets = facet_tags.find(tags["top_contact"])

    missing = []
    for name, arr in [
        ("bottom_fixed", bot_fixed_facets),
        ("top_disp", top_disp_facets),
        ("sym_left", sym_left_facets),
        ("sym_right", sym_right_facets),
        ("bottom_contact", bot_contact_facets),
        ("top_contact", top_contact_facets),
    ]:
        if arr.size == 0:
            missing.append(name)
    if missing:
        if comm.rank == 0:
            print(f"Available facet tags: {np.unique(facet_tags.values)}")
        raise RuntimeError(f"Missing facet tags: {', '.join(missing)}")

    # BCs
    bot_dofs = fem.locate_dofs_topological(V, fdim, bot_fixed_facets)
    zero_vec = fem.Constant(mesh, default_scalar_type((0.0, 0.0)))
    bc_bot = fem.dirichletbc(zero_vec, bot_dofs, V)

    top_dofs_z = fem.locate_dofs_topological(V.sub(1), fdim, top_disp_facets)
    disp_value = fem.Constant(mesh, default_scalar_type(0.0))
    bc_top = fem.dirichletbc(disp_value, top_dofs_z, V.sub(1))

    sym_left_dofs = fem.locate_dofs_topological(V.sub(0), fdim, sym_left_facets)
    bc_sym_left = fem.dirichletbc(default_scalar_type(0.0), sym_left_dofs, V.sub(0))
    sym_right_dofs = fem.locate_dofs_topological(V.sub(0), fdim, sym_right_facets)
    bc_sym_right = fem.dirichletbc(default_scalar_type(0.0), sym_right_dofs, V.sub(0))

    bcs = [bc_bot, bc_top, bc_sym_left, bc_sym_right]

    def compute_coefficients(x, coeffs):
        du.x.scatter_forward()
        contact_problem.update_contact_data(du)

    def compute_residual(x, b, coeffs):
        b.zeroEntries()
        b.ghostUpdate(addv=InsertMode.INSERT, mode=ScatterMode.FORWARD)
        contact_problem.assemble_vector(b, V)
        assemble_vector(b, F_compiled)
        if bcs:
            apply_lifting(b, [J_compiled], bcs=[bcs], x0=[x], alpha=-1.0)
        b.ghostUpdate(addv=InsertMode.ADD, mode=ScatterMode.REVERSE)
        if bcs:
            set_bc(b, bcs, x, -1.0)

    def compute_jacobian_matrix(x, a_mat, coeffs):
        a_mat.zeroEntries()
        contact_problem.assemble_matrix(a_mat, V)
        assemble_matrix(a_mat, J_compiled, bcs=bcs)
        a_mat.assemble()

    a_mat = contact_problem.create_matrix(J_compiled)
    b = create_vector(F_compiled)

    newton_solver = NewtonSolver(mesh.comm, a_mat, b, contact_problem.coeffs)
    newton_solver.set_residual(compute_residual)
    newton_solver.set_jacobian(compute_jacobian_matrix)
    newton_solver.set_coefficients(compute_coefficients)

    null_space = rigid_motions_nullspace_subdomains(
        V, cell_tags, np.unique(cell_tags.values), num_domains=2
    )
    newton_solver.A.setNearNullSpace(null_space)

    newton_solver.set_newton_options({
        "relaxation_parameter": float(args.relax),
        "atol": float(args.newton_atol),
        "rtol": float(args.newton_rtol),
        "convergence_criterion": "residual",
        "max_it": int(args.newton_max_it),
        "error_on_nonconvergence": True,
    })
    newton_solver.set_krylov_options({
        "ksp_type": str(args.ksp_type),
        "ksp_rtol": float(args.ksp_rtol),
        "ksp_atol": float(args.ksp_atol),
        "ksp_max_it": int(args.ksp_max_it),
        "pc_type": str(args.pc_type),
        "pc_gamg_type": "agg",
        "pc_gamg_agg_nsmooths": 1,
        "pc_gamg_threshold": 0.02,
        "mg_levels_ksp_type": "chebyshev",
        "mg_levels_pc_type": "sor",
        "mg_levels_ksp_max_it": 2,
    })

    update_contact = not args.no_update_contact
    steps = max(1, int(args.steps))
    total_disp = float(disp_um)
    step_size = total_disp / steps

    # Von Mises (2D plane-strain-aware). σ_yy = ν(σ_xx + σ_zz) for plane strain.
    s2d = sigma(u)
    sxx = s2d[0, 0]
    szz = s2d[1, 1]
    sxz = s2d[0, 1]
    syy = args.nu * (sxx + szz)
    s_h = (sxx + syy + szz) / 3.0
    dxx = sxx - s_h
    dyy = syy - s_h
    dzz = szz - s_h
    von_mises_expr = ufl.sqrt(3.0 / 2.0 * (dxx * dxx + dyy * dyy + dzz * dzz + 2.0 * sxz * sxz))
    von_mises_compiled = fem.Expression(von_mises_expr, V_stress.element.interpolation_points())

    xdmf_path = out_dir / "contact_results.xdmf"
    force_data = []
    reaction_vec = create_vector(F_compiled)

    if comm.rank == 0:
        print(f"Stepping {steps} increments to total disp {total_disp:.3f} µm "
              f"(step {step_size:.4f} µm)")

    with io.XDMFFile(comm, xdmf_path, "w") as xdmf:
        xdmf.write_mesh(mesh)
        prev_target = 0.0
        for step in range(1, steps + 1):
            current_target = step * step_size
            delta = current_target - prev_target
            # Direction sign: down = -z (push in), up = +z (retention pull-out)
            sign = -1.0 if args.direction == "down" else +1.0
            disp_value.value = default_scalar_type(sign * delta)

            iters, converged = newton_solver.solve(du, write_solution=False)
            du.x.scatter_forward()
            u.x.array[:] += du.x.array[:]

            von_mises.interpolate(von_mises_compiled)

            reaction_vec.zeroEntries()
            assemble_vector(reaction_vec, F_compiled)
            reaction_vec.ghostUpdate(addv=InsertMode.ADD, mode=ScatterMode.REVERSE)
            # 2D blocked layout: every 2nd DOF starting at idx 1 is z-component.
            reaction_z = 0.0
            for dof in bot_dofs:
                if dof < len(reaction_vec.array) and dof % 2 == 1:
                    reaction_z += reaction_vec.array[dof]
            reaction_z_global = comm.allreduce(reaction_z, op=MPI.SUM)

            xdmf.write_function(u, float(step))
            xdmf.write_function(von_mises, float(step))

            force_data.append((current_target, reaction_z_global))

            if comm.rank == 0:
                print(
                    f"step {step}/{steps}: disp={current_target:.3f} µm, "
                    f"Fz={reaction_z_global:.3f} µN/µm-depth, "
                    f"iters={iters}, converged={converged}"
                )
            prev_target = current_target

            # End-of-step rebuild so the NEXT step's matrix reflects new
            # contact connectivity (pairings change as bodies approach).
            # This follows the box_key_3D / christmas_tree demo pattern.
            if step < steps:
                if update_contact:
                    contact_problem.update_contact_detection(u)
                a_mat_new = contact_problem.create_matrix(J_compiled)
                a_mat_new.setNearNullSpace(null_space)
                newton_solver.set_petsc_matrix(a_mat_new)
                # warm start: carry 5% of last increment as initial guess
                du.x.array[:] = 0.05 * du.x.array[:]
                contact_problem.update_contact_data(du)

    if comm.rank == 0:
        csv_path = out_dir / "force_displacement.csv"
        with csv_path.open("w", encoding="utf-8") as f:
            f.write("displacement_um,reaction_force_z_uN_per_um_depth\n")
            for d, fz in force_data:
                f.write(f"{d:.6f},{fz:.6f}\n")
        print(f"Wrote force data to {csv_path}")
        print(f"Wrote results to {xdmf_path}")


if __name__ == "__main__":
    main()
