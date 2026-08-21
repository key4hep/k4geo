# MAIA solenoid field map

Generator: [`utils/make_maia_fieldmap.py`](../../../../utils/make_maia_fieldmap.py).
Used by the compact [`MAIA_v1.xml`](./MAIA_v1.xml) via
[`Field_Solenoid_Map_5T.xml`](./Field_Solenoid_Map_5T.xml).

## 1. Motivation

The MAIA_v0 field is DD4hep's built-in `type="solenoid"`
([`Field_Solenoid_Analytic_5T.xml`](../MAIA_v0/Field_Solenoid_Analytic_5T.xml)), whose
implementation (`dd4hep::SolenoidField`) is a two-step function in `Bz` alone:

```
Bz = +5.000 T          for r < 1678.35 mm    and |z| < zmax
Bz = -0.936941 T       for 1678.35 < r < 4113.5 mm and |z| < zmax
Bz = 0                 everywhere else
Br = 0                 everywhere
```

The map replaces it with a physically self-consistent field that has a smooth
fringe and a genuine radial component.

## 2. Methodology

### 2.1 Formulation

The problem is axisymmetric, so the field is derived from a single scalar. The
generator solves for the **stream function** `psi = r * A_phi`, where `A_phi` is
the azimuthal vector potential. The flux through a disc of radius `r` is
`Phi = 2*pi*psi`, and

```
Br = -(1/r) dpsi/dz            Bz = (1/r) dpsi/dr
```

This choice matters: `div(B) = 0` is then satisfied *identically*, by
construction, rather than approximately. Taking the phi component of
`curl(H) = J` gives a self-adjoint elliptic equation,

```
d/dr [ (nu_z / r) dpsi/dr ] + d/dz [ (nu_r / r) dpsi/dz ] = -J_phi
```

with reluctivity `nu = 1/mu`. Note the pairing: `nu_z`, the reluctivity seen by
`Bz`, multiplies the *r* derivative.

### 2.2 Discretisation and domain

* **Scheme.** Five-point finite volume. Face coefficients use harmonic-mean
  reluctivities, which is what makes the material discontinuities come out
  right rather than smeared.
* **Grid.** Uniform 25 mm out to r = 8 m and z = 8.5 m, then geometrically
  stretched (ratio 1.10) to a 40 m outer boundary. 371 x 391 = 145 061 nodes.
  The far boundary is pushed well beyond the detector because the stray field
  is not negligible at the detector edge.
* **Symmetry.** Only z >= 0 is solved. The z = 0 plane is a Neumann boundary
  (`dpsi/dz = 0`, hence `Br = 0` there exactly).
* **Boundary conditions.** `psi = 0` on the axis and on the far r and z
  boundaries.
* **Linear solve.** Sparse direct factorisation (`scipy.sparse.linalg.spsolve`).

### 2.3 Source and normalisation

The coil is a uniform azimuthal current density over its true cross section,
taken from the compact: 1530 < r < 1826.7 mm, |z| < 2307 mm (that is
`Solenoid_inner_radius + SolenoidVacuumTank_thickness`, thickness
`SolenoidCoil_thickness`, half-length `Solenoid_half_length`).

An outer loop rescales the current density until `Bz(0,0)` is exactly 5.000 T,
so the map is anchored to the same central field as the analytic description.
The converged value is **14.427 A/mm^2**, i.e. 19.75 MA total ampere-turns —
consistent with the infinite-solenoid estimate `B*L/mu0` = 18.4 MA for a 4.6 m
long 5 T solenoid.

### 2.4 Magnetic material

Iron is taken **literally from the MAIA geometry**. The only magnetic material
in MAIA_v0 is the `Steel235` absorber of the HCal:

| Region | Extent [mm] | Laminations stacked along |
| ------ | ----------- | ------------------------- |
| HCal barrel | 2126.0 < r < 4113.5, \|z\| < 2574.5 | r (cylindrical shells) |
| HCal endcap | 2575.5 < \|z\| < 4600, 310 < r < 4113.5 | z (disks) |

The steel fill factor follows from the slice list: `0.5 + 19 + 0.5` mm of
`Steel235` in a 26.5 mm layer, so **f = 0.7547** (and 75 x 26.5 mm = 1987.5 mm
= 4113.5 - 2126.0, exactly).

> **The muon system carries no iron in this model.** The region is magnetically 
> transparent and some flux returns through air. 
> The generator has a `--yoke-is-steel` flag that treats those slices as steel 
> (fill factors 0.931 barrel, 0.868 endcap) for anyone wanting the comparison.

**Constitutive law.** Froelich-Kennelly, `B = mu0*H + Js*H/(H + H0)`, a standard
closed-form model for soft magnetic materials, with `Js = 2.05 T` and initial
relative permeability 4000 (hence `H0 = 408 A/m`). This is a *generic mild-steel
curve representative of S235JR / AISI 1010, not measured S235 data*. It is
inverted analytically for `H(B)` (a quadratic in `H`).

**Anisotropic homogenisation.** A laminated stack is not an isotropic medium.
Flux running along the laminations sees the steel and the gaps in parallel;
flux crossing them sees them in series:

```
mu_par  = f*mu_Fe + (1-f)*mu0                      (~3000*mu0)
mu_perp = [ f/mu_Fe + (1-f)/mu0 ]^-1  ~ mu0/(1-f)  (~4*mu0)
```

So the HCal **barrel** is an excellent axial flux conductor while the HCal
**endcap**, whose plates are perpendicular to z, is a poor one — a factor ~750
difference that an isotropic average would erase entirely.

The steel is evaluated at its own induction, not the homogenised average.
Across the laminae `B` is continuous; along them the flux concentrates by

```
B_Fe_par = B_par / ( f + (1-f)*mu0/mu_Fe )
```

which tends to `B_par/f` while the steel is permeable and to `B_par` once it
saturates. A fixed `1/f` factor would manufacture unphysical over-saturation.
Since `mu_Fe` itself depends on the answer, the scalar relation is closed by a
few local fixed-point sweeps.

### 2.5 Nonlinear solution

Picard iteration on `nu(B)`, started from the vacuum solution and wrapped in the
current-normalisation loop of section 2.3. Two details are load-bearing:

* **Relaxation is applied to the induction the constitutive law is evaluated
  at, not to `nu`.** `nu` is an extremely stiff function of `B` around the knee,
  and relaxing it directly oscillates rather than converges.
* **Convergence is measured as the fixed-point residual** — the gap between the
  field `nu` was evaluated at and the field that came back — not as the step
  between successive solves. Under heavy under-relaxation consecutive solves are
  nearly identical long before the iteration has settled, so a step-size test
  reports false convergence.

The relaxation factor must stay at or below ~0.1; above that the steel boundary
at the HCal inner radius drives a stable limit cycle. The default is 0.05, which
converges to 1 mT. Full run time is about 4 minutes.

## 3. The ROOT file

### 3.1 Format

The file is written to the contract that k4geo's `FieldBrBz` plugin
(`detector/other/FieldMapBrBz.cpp`) expects. That plugin *infers* the whole grid
geometry from the tree — nothing is declared in XML — so the layout is not
optional.

| Property | Value |
| -------- | ----- |
| File | `fieldmaps/MAIA_fieldMap_Solenoid5T_v1.root` (757 KB) |
| MD5 | `d856800a560b56e4a7bfea0c23c8db61` |
| TTree | `ntuple` |
| Branches | `rho_mm`, `z_mm`, `Brho`, `Bz`, all `Float_t` |
| Grid | 281 (r) x 321 (z) = 90 201 entries, complete and uniform |
| Range | 0 <= r <= 7000 mm, 0 <= z <= 8000 mm, step 25 mm |
| Ordering | `RZ` — r varies fastest, both axes low-to-high |
| Units | mm and tesla (`coorUnits="mm"`, `BfieldUnits="tesla"`) |

Constraints worth restating, because violating them fails at geometry-build
time or, worse, silently:

* Branches **must** be `Float_t`; `SetBranchAddress` against a `float*` fails
  otherwise.
* The grid must be **complete** — the plugin asserts `nRho * nZ == GetEntries()`.
* Only **z >= 0** is stored. The plugin mirrors z < 0 itself (`z *= -1;
  phi += M_PI`), which is the correct symmetry for this field.
* Outside the map the plugin returns **zero**, so the map has to extend to where
  the field is negligible (see section 4).

## 4. Validation

### 4.1 Internal consistency

Printed by the generator on every run:

| Check | Result | Expected |
| ----- | ------ | -------- |
| Central field `Bz(0,0)` | 5.0000 T | 5 T by construction |
| `div(B) = 0` | exact | guaranteed by the `psi` formulation |
| Flux through z = 0, peak | 45.99 Wb | `5 T * pi * (1.71 m)^2` = 45.9 Wb |
| Residual flux at the outer boundary | 0.133 Wb | ~0; 0.3% of the total |
| Mean \|B\| in the HCal steel | 1.175 T | ~1.2 T from flux balance |
| Fraction of steel past `Js` | 16% | localised, forward region only |
| Peak \|B\| at the map edge | 9.7 mT | small, see below |

The peak steel field is 3.64 T at r = 350 mm, z = 2600 mm, where `mu_r` has
fallen to 2.3. That is not a failure: `Js` is the saturation *polarisation*, not
a ceiling on `B` — beyond it the steel stops contributing and the excess is
carried by `mu0*H`. The affected region is the small forward volume directly
behind the coil end.

The 9.7 mT at the map boundary is the size of the discontinuity introduced by
the plugin returning zero outside the map. It is a consequence of the
no-yoke-iron assumption; with an iron yoke it would be smaller.

**Grid convergence.** A 100 mm run gives `Bz(0, 2307) = 3.441 T` and a return
field of -0.909 T against 3.407 T and -0.915 T at 25 mm, i.e. ~1%.

### 4.2 Against the analytic field, through DD4hep

Sampled from the built geometry via `dd4hep.Detector.field().magneticField()`,
so this exercises the real plugin and the real interpolation, not the generator:

| r, z [mm] | Br analytic | Bz analytic | Br map | Bz map |
| --------- | ----------- | ----------- | ------ | ------ |
| 0, 0 | 0 | 5.0000 | 0 | 5.0000 |
| 0, 1000 | 0 | 5.0000 | 0 | 4.7470 |
| 0, 2000 | 0 | 5.0000 | 0 | 3.8646 |
| 0, 2307 (coil end) | 0 | 5.0000 | 0 | 3.3930 |
| 0, 2500 | 0 | 0 | 0 | 2.8831 |
| 0, 3000 | 0 | 0 | 0 | 1.2303 |
| 0, 4000 | 0 | 0 | 0 | 0.0477 |
| 500, 2307 | 0 | 5.0000 | 0.3412 | 3.5495 |
| 1400, 2000 | 0 | 5.0000 | 1.0716 | 4.2514 |
| 1400, 2307 | 0 | 5.0000 | 1.3229 | 3.3343 |
| 1400, 2600 | 0 | 0 | 1.5803 | 2.5003 |
| 2500, 0 | 0 | -0.9369 | 0 | -1.5341 |
| 3000, 0 | 0 | -0.9369 | 0 | -1.2918 |
| 4000, 0 | 0 | -0.9369 | 0 | -0.9197 |
| 5000, 0 | 0 | 0 | 0 | -0.0053 |

Reading of the table:

* The 5 T to 0 T cliff at `zmax` is replaced by a smooth fringe decaying over
  roughly 2 m.
* `Br` reaches 1.3-1.6 T near the coil ends, where the analytic field has none.
* The return field is no longer flat: it is strongest at the inner radius of the
  HCal steel (-1.53 T) and falls to -0.92 T at the outer edge, since the flux
  concentrates in the iron.
* `Br` at z = 0 is exactly zero, as the midplane symmetry requires.

**Independent cross-check on the return field.** Averaged over the annulus the
analytic model covers (1678.35 to 4113.5 mm), the map gives **-0.9148 T**
against the hand-set **-0.936941 T** — 2.4% apart. These are independent
derivations: one is a flux-balance number baked into the compact years ago, the
other falls out of the magnetostatic solve.

### 4.3 Plugin and simulation

* The plugin's inferred grid parameters were checked against what was written:
  `nRho = 281`, `rhoMax = 700 cm`, `rhoStep = 2.5 cm`, ordering `RZ`,
  both axes low-to-high.
* `ddsim` runs on both `MAIA_v0.xml` and `MAIA_v1.xml` complete with
  zero lines matching ctest's failure regex (` Exception|EXCEPTION|ERROR|Error`),
  reporting `GlobalSolenoid [solenoid]` and `GlobalSolenoidMap [FieldBrBz]`
  respectively.

### 4.4 Not yet done

A track-level comparison — single muons at several eta, analytic versus map,
comparing pT resolution and the forward track states that motivated the work —
has not been run. That is the test that would tell you whether the map changes
physics performance, and it is the natural next step.

## 5. Limitations

* **No iron in muon system**. This is the largest modelling choice in the
  document. It leaves a stray field outside the detector and essentially no
  field in the muon system.
* **Generic BH curve**, not measured S235 data.
* **Axisymmetric.** The 12-fold segmentation of the calorimeters and yoke, and
  the gaps between staves, are averaged away.
* **No coil-end detail.** The solenoid endcap side plates are not in the
  geometry (`Solenoid_o1_v01_01.xml` says so explicitly) and so are not modelled.
* **Nothing beyond the solenoid.** No anti-DID, no compensating or shielding
  solenoids, no final-focus quadrupoles.

## 6. Regenerating

```bash
scripts/run_in_container.sh python3 utils/make_maia_fieldmap.py \
    -o fieldmaps/MAIA_fieldMap_Solenoid5T_v1.root --plots /tmp/plots
```

Useful options: `--yoke-is-steel` for the iron-yoke variant, `--step` for the
grid pitch, `--central-field` for a different `Bz(0,0)`, and `--relax` /
`--max-picard` for the nonlinear iteration. Run `--help` for the full list.
