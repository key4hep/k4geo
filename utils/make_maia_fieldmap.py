#!/usr/bin/env python3
"""Generate an approximate realistic (Br, Bz) field map for the MAIA solenoid.

The output is a ROOT n-tuple in exactly the format k4geo's ``FieldBrBz`` plugin
expects (``detector/other/FieldMapBrBz.cpp``): a TTree with four ``float``
branches on a complete, uniformly spaced (rho, z) grid with rho varying fastest.

Motivation
----------
MAIA currently uses DD4hep's built-in ``type="solenoid"``, which is a two-step
function in Bz only: exactly +5 T inside the coil radius, exactly -0.937 T in
the return annulus, and exactly zero beyond |z| > zmax.  There is no Br
anywhere and no fringe field, so track states at the calorimeter endcaps sit on
a hard 5 T -> 0 T cliff.

Physics model
-------------
Axisymmetric magnetostatics solved for the stream function psi = r * A_phi::

    B_r = -(1/r) dpsi/dz            B_z = (1/r) dpsi/dr

so that div(B) = 0 holds identically.  The phi component of curl(H) = J gives
the self-adjoint elliptic equation

    d/dr [ (nu_z / r) dpsi/dr ] + d/dz [ (nu_r / r) dpsi/dz ] = -J_phi

with nu = 1/mu.  Note that nu_z, the reluctivity seen by B_z, multiplies the
*r* derivative.  It is discretised with a 5-point finite volume scheme using
harmonic-mean face reluctivities, which handles the material jumps.

The source is a uniform azimuthal current density over the coil cross section,
normalised so that Bz(0, 0) is exactly 5 T.

Iron is taken *literally from the MAIA geometry*: the only magnetic material is
the ``Steel235`` absorber of the HCal barrel and endcap.  The MAIA yoke
absorber slices are ``Air`` (see YokeBarrel_o1_v01_01.xml), so by default the
yoke is magnetically transparent.

Laminated calorimeter stacks are homogenised anisotropically: flux running
along the laminations sees mu_par = f*mu_Fe + (1-f)*mu0, flux crossing them
sees mu_perp = [f/mu_Fe + (1-f)/mu0]^-1 ~ mu0/(1-f).  That distinction matters
here: the HCal barrel is a good axial flux conductor while the HCal endcap,
whose plates are perpendicular to z, is a poor one.
"""

import argparse
import os

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

MU0 = 4.0e-7 * np.pi  # H/m

# ----------------------------------------------------------------------------
# MAIA_v0 geometry, in mm, taken from MuColl/MAIA/compact/MAIA_v0/
# ----------------------------------------------------------------------------
SOLENOID_INNER_RADIUS = 1500.0  # MAIA_v0.xml
VACUUM_TANK_THICKNESS = 30.0  # MAIA_v0.xml
COIL_THICKNESS = 296.7  # MAIA_v0.xml
SOLENOID_HALF_LENGTH = 2307.0  # MAIA_v0.xml

COIL_R_MIN = SOLENOID_INNER_RADIUS + VACUUM_TANK_THICKNESS  # 1530.0
COIL_R_MAX = COIL_R_MIN + COIL_THICKNESS  # 1826.7
COIL_Z_MAX = SOLENOID_HALF_LENGTH

HCAL_BARREL_R_MIN = 2126.0  # MAIA_v0.xml
HCAL_BARREL_R_MAX = 4113.5  # MAIA_v0.xml
HCAL_BARREL_Z_MAX = 2574.5  # MAIA_v0.xml
HCAL_ENDCAP_Z_MIN = 2575.5  # MAIA_v0.xml
HCAL_ENDCAP_Z_MAX = 4600.0  # MAIA_v0.xml
HCAL_ENDCAP_R_MIN = 310.0  # MAIA_v0.xml

YOKE_BARREL_R_MIN = 4150.0  # MAIA_v0.xml
YOKE_BARREL_R_MAX = 5895.0  # MAIA_v0.xml
YOKE_BARREL_Z_MAX = 4600.3  # MAIA_v0.xml (HCalEndcap_max_z + 3*env_safety)
YOKE_ENDCAP_Z_MIN = YOKE_BARREL_Z_MAX
YOKE_ENDCAP_Z_MAX = 6200.0  # MAIA_v0.xml
YOKE_ENDCAP_R_MIN = 450.0  # MAIA_v0.xml

# Steel fraction of an HCal layer: 0.5 + 19 + 0.5 mm of Steel235 out of a
# 0.5+19+3+0.1+0.7+0.5+2.7 = 26.5 mm layer.  Identical in barrel and endcap.
# 75 layers x 26.5 mm = 1987.5 mm = 4113.5 - 2126.0 exactly.
HCAL_STEEL_FRACTION = 20.0 / 26.5

# If the yoke absorber gaps are taken to be steel: 40.6 cm out of a 43.6 cm
# barrel layer, 19.7 cm out of a 22.7 cm endcap layer.
YOKE_BARREL_STEEL_FRACTION = 40.6 / 43.6
YOKE_ENDCAP_STEEL_FRACTION = 19.7 / 22.7

# Material tags
VACUUM = 0
LAMINATED_R = 1  # laminations stacked along r: normal is r_hat (barrel-like)
LAMINATED_Z = 2  # laminations stacked along z: normal is z_hat (endcap-like)

# ----------------------------------------------------------------------------
# Steel B(H)
# ----------------------------------------------------------------------------
# Froelich-Kennelly relation, B = mu0*H + Js*H/(H + H0), a standard closed-form
# model for soft magnetic materials.  Parameters chosen to be representative of
# mild structural steel (S235JR / AISI 1010): saturation polarisation
# Js = 2.05 T and initial relative permeability MUR_INIT = 4000, which fixes
# H0 = Js / (mu0 * (MUR_INIT - 1)) ~ 408 A/m.  This is a generic mild-steel
# curve, not measured S235 data; it reproduces the usual behaviour (mu_r halved
# at ~1.0 T, knee at ~1.5 T, B ~ 2.0 T at 2e4 A/m).
STEEL_JS = 2.05  # T
STEEL_MUR_INIT = 4000.0
STEEL_H0 = STEEL_JS / (MU0 * (STEEL_MUR_INIT - 1.0))  # A/m


def steel_nu(b):
    """Reluctivity nu = H/B [m/H] of the steel at induction ``b`` [T].

    Inverts B = mu0*H + Js*H/(H + H0) for H, which is a quadratic in H:
        mu0*H^2 + (mu0*H0 + Js - B)*H - B*H0 = 0
    """
    b = np.abs(b)
    beta = MU0 * STEEL_H0 + STEEL_JS - b
    h = (-beta + np.sqrt(beta * beta + 4.0 * MU0 * b * STEEL_H0)) / (2.0 * MU0)
    # H/B -> 1/(mu0*mur_init) as B -> 0
    return np.where(b > 1e-9, h / np.maximum(b, 1e-12), 1.0 / (MU0 * STEEL_MUR_INIT))


# ----------------------------------------------------------------------------
# Grid
# ----------------------------------------------------------------------------
def build_axis(fine_max, step, far_max, growth):
    """Uniform ``step`` spacing up to ``fine_max``, then geometric stretching."""
    fine = np.arange(0.0, fine_max + 0.5 * step, step)
    nodes = list(fine)
    delta = step
    while nodes[-1] < far_max:
        delta *= growth
        nodes.append(nodes[-1] + delta)
    return np.array(nodes)


def cell_widths(nodes):
    """Finite-volume cell widths; the first cell is a half cell (symmetry/axis)."""
    widths = np.empty_like(nodes)
    widths[1:-1] = 0.5 * (nodes[2:] - nodes[:-2])
    widths[0] = 0.5 * (nodes[1] - nodes[0])
    widths[-1] = 0.5 * (nodes[-1] - nodes[-2])
    return widths


# ----------------------------------------------------------------------------
# Materials and source
# ----------------------------------------------------------------------------
def build_regions(r_mm, z_mm, yoke_is_steel):
    """Return (material tag, steel fill fraction) arrays of shape (nr, nz)."""
    rr, zz = np.meshgrid(r_mm, z_mm, indexing="ij")
    tag = np.full(rr.shape, VACUUM, dtype=np.int8)
    frac = np.zeros(rr.shape)

    def fill(mask, material, fraction):
        tag[mask] = material
        frac[mask] = fraction

    fill(
        (rr >= HCAL_BARREL_R_MIN) & (rr <= HCAL_BARREL_R_MAX) & (zz <= HCAL_BARREL_Z_MAX),
        LAMINATED_R,
        HCAL_STEEL_FRACTION,
    )
    fill(
        (rr >= HCAL_ENDCAP_R_MIN)
        & (rr <= HCAL_BARREL_R_MAX)
        & (zz >= HCAL_ENDCAP_Z_MIN)
        & (zz <= HCAL_ENDCAP_Z_MAX),
        LAMINATED_Z,
        HCAL_STEEL_FRACTION,
    )

    if yoke_is_steel:
        fill(
            (rr >= YOKE_BARREL_R_MIN) & (rr <= YOKE_BARREL_R_MAX) & (zz <= YOKE_BARREL_Z_MAX),
            LAMINATED_R,
            YOKE_BARREL_STEEL_FRACTION,
        )
        fill(
            (rr >= YOKE_ENDCAP_R_MIN)
            & (rr <= YOKE_BARREL_R_MAX)
            & (zz >= YOKE_ENDCAP_Z_MIN)
            & (zz <= YOKE_ENDCAP_Z_MAX),
            LAMINATED_Z,
            YOKE_ENDCAP_STEEL_FRACTION,
        )

    return tag, frac


def build_source(r_mm, z_mm, current_density):
    """Uniform azimuthal current density [A/m^2] over the coil cross section."""
    rr, zz = np.meshgrid(r_mm, z_mm, indexing="ij")
    coil = (rr >= COIL_R_MIN) & (rr <= COIL_R_MAX) & (zz <= COIL_Z_MAX)
    return np.where(coil, current_density, 0.0)


def steel_induction(b_par, b_perp, f, iterations=10):
    """Induction inside the steel laminae given the homogenised field.

    Across the laminae the induction is continuous, B_Fe_perp = B_perp.  Along
    them the flux concentrates into the steel by

        B_Fe_par = B_par / (f + (1 - f) * mu0 / mu_Fe)

    which tends to B_par / f while the steel is permeable and to B_par once it
    saturates, so it cannot manufacture the unphysical over-saturation that a
    fixed 1/f factor would.  mu_Fe itself depends on the result, so the scalar
    relation is closed by a few fixed-point sweeps.
    """
    b_fe = np.hypot(b_par / f, b_perp)
    for _ in range(iterations):
        mu_fe = 1.0 / steel_nu(b_fe)
        b_fe_par = b_par / (f + (1.0 - f) * MU0 / mu_fe)
        b_fe = np.hypot(b_fe_par, b_perp)
    return b_fe


def update_reluctivity(b_r, b_z, tag, frac):
    """Anisotropic homogenised (nu_r, nu_z) for the laminated steel stacks."""
    nu_r = np.full(tag.shape, 1.0 / MU0)
    nu_z = np.full(tag.shape, 1.0 / MU0)

    for material in (LAMINATED_R, LAMINATED_Z):
        sel = tag == material
        if not np.any(sel):
            continue
        f = frac[sel]
        if material == LAMINATED_R:  # normal along r: B_z is parallel
            b_par, b_perp = b_z[sel], b_r[sel]
        else:  # normal along z: B_r is parallel
            b_par, b_perp = b_r[sel], b_z[sel]

        mu_fe = 1.0 / steel_nu(steel_induction(b_par, b_perp, f))
        mu_par = f * mu_fe + (1.0 - f) * MU0
        mu_perp = 1.0 / (f / mu_fe + (1.0 - f) / MU0)

        if material == LAMINATED_R:
            nu_z[sel] = 1.0 / mu_par
            nu_r[sel] = 1.0 / mu_perp
        else:
            nu_r[sel] = 1.0 / mu_par
            nu_z[sel] = 1.0 / mu_perp

    return nu_r, nu_z


# ----------------------------------------------------------------------------
# Solver
# ----------------------------------------------------------------------------
class Solver:
    """Finite-volume discretisation of the psi = r*A_phi equation.

    Boundary conditions: psi = 0 on the axis r = 0 and on the far r and z
    boundaries; dpsi/dz = 0 on the z = 0 symmetry plane (implicit, the first
    cell is a half cell with no south face).
    """

    def __init__(self, r, z):
        self.r, self.z = r, z
        self.nr, self.nz = len(r), len(z)
        self.drc = cell_widths(r)
        self.dzc = cell_widths(z)
        self.dr_f = np.diff(r)
        self.dz_f = np.diff(z)
        self.rf = 0.5 * (r[:-1] + r[1:])
        # r[0] = 0 is a Dirichlet node; its row is never used, but the 1/r
        # factor of the z-face coefficients would still be evaluated there.
        self.r_safe = r.copy()
        self.r_safe[0] = 0.5 * r[1]

        self.index = np.arange(self.nr * self.nz).reshape(self.nr, self.nz)
        self.dirichlet = np.zeros((self.nr, self.nz), dtype=bool)
        self.dirichlet[0, :] = True  # axis
        self.dirichlet[-1, :] = True  # far r
        self.dirichlet[:, -1] = True  # far z
        self.free = ~self.dirichlet

    def assemble(self, nu_r, nu_z):
        # Harmonic-mean face reluctivities (exact for the uniform fine grid;
        # the stretched outer region is all vacuum so the weighting is moot).
        nu_z_f = 2.0 * nu_z[:-1, :] * nu_z[1:, :] / (nu_z[:-1, :] + nu_z[1:, :])
        nu_r_f = 2.0 * nu_r[:, :-1] * nu_r[:, 1:] / (nu_r[:, :-1] + nu_r[:, 1:])

        a_e = self.dzc[None, :] * (nu_z_f / self.rf[:, None]) / self.dr_f[:, None]
        a_n = self.drc[:, None] * (nu_r_f / self.r_safe[:, None]) / self.dz_f[None, :]

        diag = np.zeros((self.nr, self.nz))
        diag[:-1, :] += a_e
        diag[1:, :] += a_e
        diag[:, :-1] += a_n
        diag[:, 1:] += a_n

        rows = [self.index[self.free], self.index[self.dirichlet]]
        cols = [self.index[self.free], self.index[self.dirichlet]]
        vals = [diag[self.free], np.ones(int(self.dirichlet.sum()))]

        # Couplings are dropped where either end is a Dirichlet node: the row
        # is replaced by the identity, and psi = 0 there contributes nothing to
        # the right-hand side of its neighbours.
        pair_r = self.free[:-1, :] & self.free[1:, :]
        rows += [self.index[:-1, :][pair_r], self.index[1:, :][pair_r]]
        cols += [self.index[1:, :][pair_r], self.index[:-1, :][pair_r]]
        vals += [-a_e[pair_r], -a_e[pair_r]]

        pair_z = self.free[:, :-1] & self.free[:, 1:]
        rows += [self.index[:, :-1][pair_z], self.index[:, 1:][pair_z]]
        cols += [self.index[:, 1:][pair_z], self.index[:, :-1][pair_z]]
        vals += [-a_n[pair_z], -a_n[pair_z]]

        n = self.nr * self.nz
        return sp.csc_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)
        )

    def solve(self, nu_r, nu_z, j_phi):
        rhs = j_phi * self.drc[:, None] * self.dzc[None, :]
        rhs[self.dirichlet] = 0.0
        psi = spla.spsolve(self.assemble(nu_r, nu_z), rhs.ravel())
        return psi.reshape(self.nr, self.nz)

    def field(self, psi):
        """B from psi; on the axis Bz = lim 2*psi/r^2 and Br = 0."""
        b_z = np.gradient(psi, self.r, axis=0)
        b_r = -np.gradient(psi, self.z, axis=1)
        b_z[1:, :] /= self.r[1:, None]
        b_r[1:, :] /= self.r[1:, None]
        b_z[0, :] = 2.0 * psi[1, :] / self.r[1] ** 2
        b_r[0, :] = 0.0
        # np.gradient falls back to a one-sided difference on the z = 0 edge,
        # but psi is even in z there, so the central difference and hence Br
        # vanish identically.
        b_r[:, 0] = 0.0
        return b_r, b_z


def solve_field(solver, tag, frac, current_density, target_bz, args):
    """Picard loop on nu(B), wrapped in a loop that normalises Bz(0, 0).

    The relaxation is applied to the induction that the constitutive law is
    evaluated at rather than to nu itself.  nu is an extremely stiff function of
    B around the knee, so relaxing nu directly oscillates instead of converging.
    """
    j_phi = build_source(solver.r * 1e3, solver.z * 1e3, current_density)

    # Start from the vacuum solution.  Initialising the working field to zero
    # would put the steel at its initial permeability everywhere, which drags
    # the whole return flux into the calorimeter on the first solve and takes
    # many iterations to unwind.
    nu_vacuum = np.full(tag.shape, 1.0 / MU0)
    psi = solver.solve(nu_vacuum, nu_vacuum, j_phi)
    b_r, b_z = solver.field(psi)
    b_eff_r, b_eff_z = b_r.copy(), b_z.copy()
    converged = False

    for scaling in range(args.max_scalings):
        for it in range(args.max_picard):
            nu_r, nu_z = update_reluctivity(b_eff_r, b_eff_z, tag, frac)
            psi = solver.solve(nu_r, nu_z, j_phi)
            b_r, b_z = solver.field(psi)
            # The convergence measure is the fixed-point residual, the gap
            # between the field the constitutive law was evaluated at and the
            # field that came back.  Measuring the step between successive
            # solves instead would report false convergence, since heavy
            # under-relaxation makes consecutive solves nearly identical long
            # before the iteration has actually settled.
            residual = np.maximum(np.abs(b_r - b_eff_r), np.abs(b_z - b_eff_z))
            change = float(residual.max())
            hot = np.unravel_index(int(np.argmax(residual)), residual.shape)
            b_eff_r += args.relax * (b_r - b_eff_r)
            b_eff_z += args.relax * (b_z - b_eff_z)
            converged = change < args.tolerance
            if converged or it % 10 == 0:
                print(
                    f"    picard {it:3d}  max |dB| = {change * 1e3:9.4f} mT "
                    f"at r = {solver.r[hot[0]] * 1e3:.0f} mm, z = {solver.z[hot[1]] * 1e3:.0f} mm"
                )
            if converged:
                break
        if not converged:
            print(
                f"    WARNING: Picard did not reach {args.tolerance * 1e3:g} mT "
                f"in {args.max_picard} iterations"
            )

        bz0 = b_z[0, 0]
        print(f"  scaling pass {scaling}: Bz(0,0) = {bz0:.6f} T")
        if abs(bz0 - target_bz) < 1e-4:
            break
        scale = target_bz / bz0
        j_phi *= scale
        current_density *= scale

    if not converged:
        raise SystemExit("solver did not converge; lower --relax or raise --max-picard")
    return psi, b_r, b_z, current_density


# ----------------------------------------------------------------------------
# Diagnostics and output
# ----------------------------------------------------------------------------
def report(solver, psi, b_r, b_z, tag, frac, args):
    r_mm, z_mm = solver.r * 1e3, solver.z * 1e3
    b_mag = np.hypot(b_r, b_z)

    def at(r_target, z_target):
        i = int(np.argmin(np.abs(r_mm - r_target)))
        j = int(np.argmin(np.abs(z_mm - z_target)))
        return b_r[i, j], b_z[i, j]

    print("\n--- field summary -----------------------------------------------")
    print(f"  Bz(r=0, z=0)            = {b_z[0, 0]:8.4f} T")
    print(f"  Bz(r=0, z=1000 mm)      = {at(0, 1000)[1]:8.4f} T")
    print(f"  Bz(r=0, z=2307 mm)      = {at(0, 2307)[1]:8.4f} T   (coil end)")
    print(f"  Bz(r=0, z=3000 mm)      = {at(0, 3000)[1]:8.4f} T")
    print(f"  Bz(r=0, z=6000 mm)      = {at(0, 6000)[1]:8.4f} T")
    br, bz = at(1400, 2307)
    print(f"  B(r=1400, z=2307 mm)    = Br {br:7.4f} T, Bz {bz:7.4f} T   (coil end)")
    br, bz = at(1400, 2000)
    print(f"  B(r=1400, z=2000 mm)    = Br {br:7.4f} T, Bz {bz:7.4f} T")

    # Mean return field over the annulus the analytic model covers, at z = 0.
    ann = (r_mm >= 1678.35) & (r_mm <= 4113.5)
    area = np.pi * (4113.5e-3**2 - 1678.35e-3**2)
    flux = 2.0 * np.pi * (psi[ann, 0][-1] - psi[ann, 0][0])
    print(f"  <Bz> over return annulus = {flux / area:8.4f} T   (analytic: -0.9369 T)")

    # Flux balance: 2*pi*psi is the flux through the z = 0 disc of radius r.
    print(f"  peak flux through z=0    = {2.0 * np.pi * psi[:, 0].max():8.4f} Wb")
    print(f"  residual flux at r_max   = {2.0 * np.pi * psi[-2, 0]:8.4f} Wb   (should be ~0)")

    steel = tag != VACUUM
    if np.any(steel):
        b_par = np.where(tag == LAMINATED_R, b_z, b_r)
        b_perp = np.where(tag == LAMINATED_R, b_r, b_z)
        b_fe = np.zeros(tag.shape)
        b_fe[steel] = steel_induction(b_par[steel], b_perp[steel], frac[steel])
        i, j = np.unravel_index(int(np.argmax(b_fe)), b_fe.shape)
        # Js is the saturation polarisation, not a ceiling on B: beyond it the
        # steel simply stops contributing and the excess is carried by mu0*H.
        # What matters is how much of the steel is pushed past the knee.
        above = float((b_fe[steel] > STEEL_JS).mean())
        print(f"  mean |B| in the steel    = {b_fe[steel].mean():8.4f} T")
        print(
            f"  peak |B| in the steel    = {b_fe.max():8.4f} T   at "
            f"r = {r_mm[i]:.0f} mm, z = {z_mm[j]:.0f} mm, "
            f"mu_r = {1.0 / (steel_nu(b_fe[i, j]) * MU0):.1f}"
        )
        print(f"  steel past saturation    = {above * 100:8.1f} %   (Js = {STEEL_JS} T)")

    # The FieldBrBz plugin returns zero outside the map, so the field at the
    # map edge is the size of the discontinuity that introduces.
    i_max = int(np.argmin(np.abs(r_mm - args.r_max)))
    j_max = int(np.argmin(np.abs(z_mm - args.z_max)))
    edge = max(b_mag[i_max, : j_max + 1].max(), b_mag[: i_max + 1, j_max].max())
    print(f"  peak |B| on the map edge = {edge * 1e3:8.2f} mT  (step at the map boundary)")
    print("-----------------------------------------------------------------\n")


def write_map(path, r_mm, z_mm, b_r, b_z):
    """Write the TTree that FieldBrBz reads: float branches, rho varying fastest."""
    import uproot

    rr, zz = np.meshgrid(r_mm, z_mm, indexing="ij")
    with uproot.recreate(path) as handle:
        handle["ntuple"] = {
            "rho_mm": rr.ravel(order="F").astype(np.float32),
            "z_mm": zz.ravel(order="F").astype(np.float32),
            "Brho": b_r.ravel(order="F").astype(np.float32),
            "Bz": b_z.ravel(order="F").astype(np.float32),
        }
    print(f"wrote {path}  ({len(r_mm)} x {len(z_mm)} = {r_mm.size * z_mm.size} entries)")


def write_plots(directory, r_mm, z_mm, b_r, b_z):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(directory, exist_ok=True)
    envelopes = [
        ((COIL_R_MIN, COIL_R_MAX), (0.0, COIL_Z_MAX), "coil"),
        ((HCAL_BARREL_R_MIN, HCAL_BARREL_R_MAX), (0.0, HCAL_BARREL_Z_MAX), "HCal barrel"),
        (
            (HCAL_ENDCAP_R_MIN, HCAL_BARREL_R_MAX),
            (HCAL_ENDCAP_Z_MIN, HCAL_ENDCAP_Z_MAX),
            "HCal endcap",
        ),
        ((YOKE_BARREL_R_MIN, YOKE_BARREL_R_MAX), (0.0, YOKE_BARREL_Z_MAX), "yoke barrel"),
        (
            (YOKE_ENDCAP_R_MIN, YOKE_BARREL_R_MAX),
            (YOKE_ENDCAP_Z_MIN, YOKE_ENDCAP_Z_MAX),
            "yoke endcap",
        ),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].plot(z_mm, b_z[0, :], label="Bz(r=0)")
    axes[0].axvline(COIL_Z_MAX, ls=":", c="grey")
    axes[0].axhline(0.0, lw=0.5, c="k")
    axes[0].set_xlabel("z [mm]")
    axes[0].set_ylabel("Bz [T]")
    axes[0].set_title("on-axis field")
    axes[0].legend()
    axes[1].plot(r_mm, b_z[:, 0], label="Bz(z=0)")
    axes[1].axhline(-0.936941, ls="--", c="r", label="analytic return field")
    axes[1].axhline(0.0, lw=0.5, c="k")
    axes[1].set_xlabel("r [mm]")
    axes[1].set_ylabel("Bz [T]")
    axes[1].set_title("midplane field")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(directory, "maia_field_profiles.png"), dpi=130)
    plt.close(fig)

    for name, values, limit in (
        ("Bz", b_z, 5.0),
        ("Br", b_r, 1.5),
        ("Bmag", np.hypot(b_r, b_z), 5.0),
    ):
        fig, ax = plt.subplots(figsize=(8, 5))
        mesh = ax.pcolormesh(
            z_mm, r_mm, values, cmap="RdBu_r", vmin=-limit, vmax=limit, shading="nearest"
        )
        for (r0, r1), (z0, z1), _ in envelopes:
            ax.add_patch(plt.Rectangle((z0, r0), z1 - z0, r1 - r0, fill=False, ec="k", lw=0.8))
        fig.colorbar(mesh, ax=ax, label=f"{name} [T]")
        ax.set_xlabel("z [mm]")
        ax.set_ylabel("r [mm]")
        ax.set_title(f"MAIA {name}")
        fig.tight_layout()
        fig.savefig(os.path.join(directory, f"maia_field_{name}.png"), dpi=130)
        plt.close(fig)

    print(f"wrote diagnostic plots to {directory}")


# ----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("-o", "--output", default="fieldmaps/maia_fieldMap_Solenoid5T_v1.root")
    parser.add_argument("--plots", default=None, help="directory for diagnostic PNGs")
    parser.add_argument("--central-field", type=float, default=5.0, help="Bz(0,0) [T]")
    parser.add_argument("--step", type=float, default=25.0, help="map grid step [mm]")
    parser.add_argument("--r-max", type=float, default=7000.0, help="map extent in r [mm]")
    parser.add_argument("--z-max", type=float, default=8000.0, help="map extent in z [mm]")
    parser.add_argument(
        "--solve-r-max", type=float, default=8000.0, help="fine grid extent in r [mm]"
    )
    parser.add_argument(
        "--solve-z-max", type=float, default=8500.0, help="fine grid extent in z [mm]"
    )
    parser.add_argument("--far-max", type=float, default=40000.0, help="outer boundary [mm]")
    parser.add_argument("--growth", type=float, default=1.10, help="stretched-grid growth ratio")
    # The steel boundary at the HCal inner radius drives a Picard limit cycle
    # for relaxation above ~0.1, so the default is deliberately heavy.
    parser.add_argument("--relax", type=float, default=0.05, help="Picard under-relaxation on B")
    parser.add_argument(
        "--tolerance", type=float, default=1e-3, help="Picard convergence on |dB| [T]"
    )
    parser.add_argument("--max-picard", type=int, default=400)
    parser.add_argument("--max-scalings", type=int, default=6)
    parser.add_argument(
        "--yoke-is-steel",
        action="store_true",
        help="treat the yoke absorber gaps as steel instead of the Air the geometry declares",
    )
    args = parser.parse_args()

    r_mm = build_axis(args.solve_r_max, args.step, args.far_max, args.growth)
    z_mm = build_axis(args.solve_z_max, args.step, args.far_max, args.growth)
    solver = Solver(r_mm * 1e-3, z_mm * 1e-3)
    print(
        f"grid: {solver.nr} x {solver.nz} = {solver.nr * solver.nz} nodes, "
        f"fine step {args.step} mm, outer boundary {args.far_max / 1e3:.0f} m"
    )

    tag, frac = build_regions(r_mm, z_mm, args.yoke_is_steel)
    print(
        f"magnetic material: {'HCal steel + yoke steel' if args.yoke_is_steel else 'HCal steel only'}"
        f" ({int((tag != VACUUM).sum())} nodes)"
    )

    # Infinite-solenoid estimate as the starting current density.
    current_density = args.central_field / (MU0 * COIL_THICKNESS * 1e-3)
    psi, b_r, b_z, current_density = solve_field(
        solver, tag, frac, current_density, args.central_field, args
    )
    print(
        f"coil current density = {current_density * 1e-6:.3f} A/mm^2, "
        f"total current = {current_density * COIL_THICKNESS * 2 * COIL_Z_MAX * 1e-6 * 1e-6:.3f} MA"
    )

    report(solver, psi, b_r, b_z, tag, frac, args)

    n_r = int(round(args.r_max / args.step)) + 1
    n_z = int(round(args.z_max / args.step)) + 1
    out_r, out_z = r_mm[:n_r], z_mm[:n_z]
    assert np.isclose(out_r[-1], args.r_max) and np.isclose(out_z[-1], args.z_max)
    write_map(args.output, out_r, out_z, b_r[:n_r, :n_z], b_z[:n_r, :n_z])

    if args.plots:
        write_plots(args.plots, out_r, out_z, b_r[:n_r, :n_z], b_z[:n_r, :n_z])


if __name__ == "__main__":
    main()
