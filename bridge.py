import firedrake as fd
import firedrake_adjoint as fda
from firedrake import (
    inner,
    sqrt,
    jump,
    dx,
    ds,
    dS,
    sym,
    nabla_grad,
    tr,
    Identity,
)
from pyMMAopt import MMASolver, ReducedInequality
import itertools
import argparse
from penalization import ramp
from solver_parameters import gamg_parameters, hypre_params


def compliance_bridge():
    DECK = 0
    DOMAIN = 1
    SUPPORT = 2
    LOAD = 3
    SYMMETRY_X = 4
    SYMMETRY_Y = 5

    parser = argparse.ArgumentParser(description="Bridge design")
    parser.add_argument(
        "--output_dir",
        dest="output_dir",
        type=str,
        action="store",
        default="./",
        help="Output directory",
    )
    opts = parser.parse_args()
    output_dir = opts.output_dir
    # Elasticity parameters
    E, nu = 1.0, 0.3
    mu, lmbda = (
        fd.Constant(E / (2 * (1 + nu))),
        fd.Constant(E * nu / ((1 + nu) * (1 - 2 * nu))),
    )

    mesh = fd.Mesh("./bridge.msh")

    RHO = fd.FunctionSpace(mesh, "DG", 0)
    rho = fd.Function(RHO).assign(fd.Constant(0.5))

    af, b = fd.TrialFunction(RHO), fd.TestFunction(RHO)

    filter_radius = fd.Constant(2e-4)
    x, y, z = fd.SpatialCoordinate(mesh)
    with fda.stop_annotating():
        x_ = fd.interpolate(x, RHO)
        y_ = fd.interpolate(y, RHO)
        z_ = fd.interpolate(z, RHO)
    Delta_h = sqrt(jump(x_) ** 2 + jump(y_) ** 2 + jump(z_) ** 2)
    aH = filter_radius * jump(af) / Delta_h * jump(b) * dS + af * b * dx
    LH = rho * b * dx

    rhof = fd.Function(RHO)
    problem_filter = fd.LinearVariationalProblem(aH, LH, rhof)
    solver_filter = fd.LinearVariationalSolver(
        problem_filter, solver_parameters=hypre_params
    )
    solver_filter.solve()
    rhofControl = fda.Control(rhof)

    H1 = fd.VectorElement("CG", mesh.ufl_cell(), 1)
    W = fd.FunctionSpace(mesh, H1)

    x, y, z = fd.SpatialCoordinate(mesh)
    modes = [fd.Function(W) for _ in range(6)]
    modes[0].interpolate(fd.Constant([1, 0, 0]))
    modes[1].interpolate(fd.Constant([0, 1, 0]))
    modes[2].interpolate(fd.Constant([0, 0, 1]))
    modes[3].interpolate(fd.as_vector([0, z, -y]))
    modes[4].interpolate(fd.as_vector([-z, 0, x]))
    modes[5].interpolate(fd.as_vector([y, -x, 0]))
    nullmodes = fd.VectorSpaceBasis(modes)
    # Make sure they're orthonormal.
    nullmodes.orthonormalize()

    u = fd.TrialFunction(W)
    v = fd.TestFunction(W)

    def epsilon(u):
        return sym(nabla_grad(u))

    def sigma(v):
        return 2.0 * mu * epsilon(v) + lmbda * tr(epsilon(v)) * Identity(3)

    # Variational forms
    a = inner(ramp(rhof, ramp_p=10.0, val_0=1e-5) * sigma(u), nabla_grad(v)) * dx(
        DOMAIN
    ) + inner(sigma(u), nabla_grad(v)) * dx(DECK)
    t = fd.Constant((0.0, 0.0, -1.0))
    L = inner(t, v) * ds(LOAD)

    # Dirichlet BCs
    bc1 = fd.DirichletBC(W, fd.Constant((0.0, 0.0, 0.0)), SUPPORT)
    bc2 = fd.DirichletBC(W.sub(0), fd.Constant(0.0), SYMMETRY_X)
    bc3 = fd.DirichletBC(W.sub(1), fd.Constant(0.0), SYMMETRY_Y)

    u_sol = fd.Function(W)
    problem = fd.LinearVariationalProblem(a, L, u_sol, bcs=[bc1, bc2, bc3])
    solver = fd.LinearVariationalSolver(
        problem, near_nullspace=nullmodes, solver_parameters=gamg_parameters
    )
    solver.solve()

    # Cost function
    J = fd.assemble(fd.Constant(1e1) * inner(t, u_sol) * ds(LOAD))
    # Constraint
    VolPen = fd.assemble(rhof * dx(DOMAIN))
    with fda.stop_annotating():
        total_vol = fd.assemble(fd.Constant(1.0) * dx(DOMAIN, domain=mesh))
        Vlimit = 0.15 * total_vol
    VolControl = fda.Control(VolPen)

    # Plotting
    global_counter1 = itertools.count()
    phi_pvd = fd.File(f"{output_dir}/bridge_evolution.pvd", target_continuity=fd.H1)
    rho_viz_f = fd.Function(RHO, name="rho")

    def deriv_cb(design):
        iter = next(global_counter1)
        if iter % 10 == 0:
            rho_viz_f.assign(rhofControl.tape_value())
            fd.par_loop(
                ("{[i] : 0 <= i < f.dofs}", "f[i, 0] = 1.0"),
                dx(DECK),
                {"f": (rho_viz_f, fd.WRITE)},
                is_loopy_kernel=True,
            )
            phi_pvd.write(rho_viz_f)

    c = fda.Control(rho)
    Jhat = fda.ReducedFunctional(J, c, derivative_cb_pre=deriv_cb)
    Volhat = fda.ReducedFunctional(VolPen, c)

    lb = 0.0
    ub = 1.0
    problem = fda.MinimizationProblem(
        Jhat,
        bounds=(lb, ub),
        constraints=[ReducedInequality(Volhat, Vlimit, VolControl)],
    )

    parameters_mma = {
        "move": 0.2,
        "maximum_iterations": 200,
        "m": 1,
        "IP": 0,
        "tol": 1e-6,
        "accepted_tol": 1e-4,
        "norm": "L2",
        "gcmma": False,
    }
    solver = MMASolver(problem, parameters=parameters_mma)

    rho_opt = solver.solve()


if __name__ == "__main__":
    compliance_bridge()
