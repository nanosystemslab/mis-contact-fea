#!/usr/bin/env python3
"""
Quasi-static two-body contact solve using dolfinx_contact.
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


def _load_tags(path: Path) -> dict[str, int]:
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    return data["tags"]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run a displacement-controlled contact simulation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--mesh-dir", type=Path, default=Path("mesh_out"))
    parser.add_argument("--mesh", type=Path, default=None, help="Use mesh.xdmf or mesh.msh")
    parser.add_argument("--facet-tags", type=Path, default=None)
    parser.add_argument("--tags", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=Path("results_contact"))
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--disp-mm", type=float, default=None)
    parser.add_argument("--E-mpa", type=float, default=200000.0)
    parser.add_argument("--nu", type=float, default=0.3)
    parser.add_argument("--fric", type=float, default=0.3)
    parser.add_argument("--q-degree", type=int, default=5)
    parser.add_argument("--gamma-scale", type=float, default=10.0)
    parser.add_argument("--no-update-contact", action="store_true")
    parser.add_argument("--newton-max-it", type=int, default=40)
    parser.add_argument("--newton-rtol", type=float, default=1e-7)
    parser.add_argument("--newton-atol", type=float, default=1e-7)
    parser.add_argument("--relax", type=float, default=1.0)
    parser.add_argument("--ksp-rtol", type=float, default=1e-6)
    parser.add_argument("--ksp-atol", type=float, default=1e-10)
    parser.add_argument("--ksp-max-it", type=int, default=2000)
    parser.add_argument("--ksp-type", type=str, default="fgmres")
    parser.add_argument("--pc-type", type=str, default="gamg")
    parser.add_argument(
        "--no-structure-pin",
        action="store_true",
        help="Disable extra pin DOFs on the structure (not recommended).",
    )
    args = parser.parse_args()

    # Validate material parameters
    if args.E_mpa <= 0:
        raise ValueError(f"Young's modulus must be positive, got {args.E_mpa}")
    if not (0 < args.nu < 0.5):
        raise ValueError(f"Poisson's ratio must be in (0, 0.5), got {args.nu}")
    if args.fric < 0:
        raise ValueError(f"Friction coefficient must be non-negative, got {args.fric}")

    mesh_dir = args.mesh_dir.resolve()
    default_mesh = mesh_dir / "mesh.msh"
    if not default_mesh.exists():
        default_mesh = mesh_dir / "mesh.xdmf"
    mesh_path = (args.mesh or default_mesh).resolve()
    facet_path = (args.facet_tags or mesh_dir / "facet_tags.xdmf").resolve()
    tags_path = (args.tags or mesh_dir / "mesh_tags.json").resolve()

    # Validate input files exist
    if not mesh_path.exists():
        raise FileNotFoundError(f"Mesh file not found: {mesh_path}")
    if not tags_path.exists():
        raise FileNotFoundError(f"Tags file not found: {tags_path}")

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    tags = _load_tags(tags_path)
    with tags_path.open("r", encoding="utf-8") as f:
        tag_meta = json.load(f)

    disp_mm = args.disp_mm
    if disp_mm is None:
        # Use pass_through_mm if available, otherwise fall back to clearance
        disp_mm = float(tag_meta.get("pass_through_mm", tag_meta.get("clearance_mm", 30.0)))
        if MPI.COMM_WORLD.rank == 0:
            print(f"Auto displacement: {disp_mm:.2f} mm (pass-through)")

    comm = MPI.COMM_WORLD
    if mesh_path.suffix == ".msh":
        mesh, cell_tags, facet_tags = gmshio.read_from_msh(
            str(mesh_path), comm, 0, gdim=3
        )
    else:
        with io.XDMFFile(comm, mesh_path, "r") as xdmf:
            mesh = xdmf.read_mesh(name="Grid")
            cell_tags = xdmf.read_meshtags(mesh, name="cell_marker")

        tdim = mesh.topology.dim
        fdim = tdim - 1
        mesh.topology.create_connectivity(fdim, tdim)
        with io.XDMFFile(comm, facet_path, "r") as xdmf:
            facet_tags = xdmf.read_meshtags(mesh, name="facet_marker")

    tdim = mesh.topology.dim
    fdim = tdim - 1
    mesh.topology.create_connectivity(fdim, tdim)
    mesh.topology.create_connectivity(0, tdim)

    V = fem.functionspace(mesh, ("CG", 1, (mesh.geometry.dim,)))
    V0 = fem.functionspace(mesh, ("DG", 0))
    V_stress = fem.functionspace(mesh, ("DG", 0))  # For von Mises stress output

    mu_func, lambda_func = lame_parameters(True)
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
    data = np.array([tags["base_contact"], tags["structure_contact"]], dtype=np.int32)
    offsets = np.array([0, 2], dtype=np.int32)
    surfaces = adjacencylist(data, offsets)
    search_mode = [ContactMode.ClosestPoint, ContactMode.ClosestPoint]

    contact_problem = ContactProblem(
        [facet_tags],
        surfaces,
        contact_pairs,
        mesh,
        args.q_degree,
        search_mode,
    )
    contact_problem.generate_contact_data(
        FrictionLaw.Coulomb,
        V,
        {"u": u, "du": du, "mu": mu_dg, "lambda": lmbda_dg, "fric": fric_dg},
        args.E_mpa * args.gamma_scale,
        -1.0,
    )

    base_facets = facet_tags.find(tags["base_fixed"])
    struct_contact_facets = facet_tags.find(tags["structure_contact"])
    base_contact_facets = facet_tags.find(tags["base_contact"])
    top_facets = facet_tags.find(tags["structure_top"])
    if base_facets.size == 0 or top_facets.size == 0 or base_contact_facets.size == 0 or struct_contact_facets.size == 0:
        if comm.rank == 0:
            available = np.unique(facet_tags.values)
            print(f"Available facet tags: {available}")
        missing = []
        if base_facets.size == 0:
            missing.append("base_fixed")
        if base_contact_facets.size == 0:
            missing.append("base_contact")
        if struct_contact_facets.size == 0:
            missing.append("structure_contact")
        if top_facets.size == 0:
            missing.append("structure_top")
        raise RuntimeError(
            f"Missing facet tags: {', '.join(missing)}; check mesh tags or increase --surface-tol-mm."
        )

    base_dofs = fem.locate_dofs_topological(V, fdim, base_facets)
    zero = fem.Constant(mesh, default_scalar_type((0.0, 0.0, 0.0)))
    bc_base = fem.dirichletbc(zero, base_dofs, V)

    if top_facets.size == 0:
        raise RuntimeError(
            "No facets found for structure_top tag; check mesh tags or increase --surface-tol-mm."
        )
    top_dofs = fem.locate_dofs_topological(V.sub(2), fdim, top_facets)
    disp = fem.Constant(mesh, default_scalar_type(0.0))
    bc_top = fem.dirichletbc(disp, top_dofs, V.sub(2))

    bcs = [bc_base, bc_top]
    if not args.no_structure_pin:
        mesh.topology.create_connectivity(fdim, 0)
        f_to_v = mesh.topology.connectivity(fdim, 0)
        top_vertices = np.unique(
            np.concatenate([f_to_v.links(int(f)) for f in np.array(top_facets, dtype=np.int32)])
        )
        if top_vertices.size >= 2:
            coords = mesh.geometry.x
            vx = coords[top_vertices, 0]
            vy = coords[top_vertices, 1]
            v_min = top_vertices[int(np.argmin(vx))]
            v_max = top_vertices[int(np.argmax(vx))]
            if v_min == v_max:
                v_max = top_vertices[int(np.argmax(vy))]
            dof_x_min = fem.locate_dofs_topological(V.sub(0), 0, np.array([v_min], dtype=np.int32))
            dof_y_min = fem.locate_dofs_topological(V.sub(1), 0, np.array([v_min], dtype=np.int32))
            dof_x_max = fem.locate_dofs_topological(V.sub(0), 0, np.array([v_max], dtype=np.int32))
            bc_pin_x = fem.dirichletbc(default_scalar_type(0.0), dof_x_min, V.sub(0))
            bc_pin_y = fem.dirichletbc(default_scalar_type(0.0), dof_y_min, V.sub(1))
            bc_pin_x2 = fem.dirichletbc(default_scalar_type(0.0), dof_x_max, V.sub(0))
            bcs.extend([bc_pin_x, bc_pin_y, bc_pin_x2])
        elif top_vertices.size == 1:
            v_min = int(top_vertices[0])
            dof_x_min = fem.locate_dofs_topological(V.sub(0), 0, np.array([v_min], dtype=np.int32))
            dof_y_min = fem.locate_dofs_topological(V.sub(1), 0, np.array([v_min], dtype=np.int32))
            bc_pin_x = fem.dirichletbc(default_scalar_type(0.0), dof_x_min, V.sub(0))
            bc_pin_y = fem.dirichletbc(default_scalar_type(0.0), dof_y_min, V.sub(1))
            bcs.extend([bc_pin_x, bc_pin_y])

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

    newton_solver.set_newton_options(
        {
            "relaxation_parameter": float(args.relax),
            "atol": float(args.newton_atol),
            "rtol": float(args.newton_rtol),
            "convergence_criterion": "residual",
            "max_it": int(args.newton_max_it),
            "error_on_nonconvergence": True,
        }
    )

    newton_solver.set_krylov_options(
        {
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
        }
    )

    update_contact = not args.no_update_contact
    steps = max(1, int(args.steps))
    total_disp = float(disp_mm)
    step_size = total_disp / steps

    # Von Mises stress expression
    s = sigma(u) - (1.0 / 3.0) * ufl.tr(sigma(u)) * ufl.Identity(mesh.geometry.dim)
    von_mises_expr = ufl.sqrt(3.0 / 2.0 * ufl.inner(s, s))
    von_mises_compiled = fem.Expression(von_mises_expr, V_stress.element.interpolation_points())

    xdmf_path = out_dir / "contact_results.xdmf"
    force_data = []  # Store (disp, Fz) for CSV output
    reaction_vec = create_vector(F_compiled)  # Reusable vector for reaction force
    with io.XDMFFile(comm, xdmf_path, "w") as xdmf:
        xdmf.write_mesh(mesh)
        prev_target = 0.0
        for step in range(1, steps + 1):
            current_target = step * step_size
            delta = current_target - prev_target
            disp.value = default_scalar_type(-delta)
            du.x.array[:] = 0.0

            if update_contact:
                contact_problem.update_contact_detection(u)

            iters, converged = newton_solver.solve(du, write_solution=False)
            du.x.scatter_forward()
            u.x.array[:] += du.x.array[:]

            # Compute von Mises stress
            von_mises.interpolate(von_mises_compiled)

            # Compute reaction force at base (z-component)
            reaction_vec.zeroEntries()
            assemble_vector(reaction_vec, F_compiled)
            reaction_vec.ghostUpdate(addv=InsertMode.ADD, mode=ScatterMode.REVERSE)
            # Sum z-component of reaction at base_fixed DOFs
            reaction_z = 0.0
            for dof in base_dofs:
                if dof < len(reaction_vec.array):
                    # Every 3rd DOF starting from index 2 is z-component
                    if dof % 3 == 2:
                        reaction_z += reaction_vec.array[dof]
            reaction_z_global = comm.allreduce(reaction_z, op=MPI.SUM)

            xdmf.write_function(u, float(step))
            xdmf.write_function(von_mises, float(step))

            force_data.append((current_target, reaction_z_global))

            if comm.rank == 0:
                print(
                    f"step {step}/{steps}: disp={current_target:.3f} mm, "
                    f"Fz={reaction_z_global:.2f} N, iters={iters}, converged={converged}"
                )
            prev_target = current_target

    # Write force-displacement data to CSV
    if comm.rank == 0:
        csv_path = out_dir / "force_displacement.csv"
        with csv_path.open("w", encoding="utf-8") as f:
            f.write("displacement_mm,reaction_force_z_N\n")
            for d, fz in force_data:
                f.write(f"{d:.6f},{fz:.6f}\n")
        print(f"Wrote force data to {csv_path}")
        print(f"Wrote results to {xdmf_path}")


if __name__ == "__main__":
    main()
