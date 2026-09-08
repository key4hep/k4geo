# Geometry Review: MuSIC_v3

Review date: 2026-05-14
Reviewer: Claude (Opus 4.7)
Scope: all XML files under `MuColl/MuSIC/compact/MuSIC_v3/`

This review extends `../MuSIC_v2/geometry_review_first.md` (2026-05-13) and
checks what is new in v3 (Forward Tracker, Nozzle_v2 with gaps, regenerated
`text_description.md`).


## Overall structure

`MuSIC_v3.xml` orchestrates 13 sub-detector XMLs plus elements/materials. XML
is well-formed. All sub-includes resolve. Geometry compiles and exports with
the local k4geo plugin libraries pre-loaded (see `export_geometry.py`).

Component summary (envelope dimensions from `MuSIC_v3.xml`):

| Subdetector    | r [mm]              | half-z / z-range [mm] | Notes                                                                |
|----------------|---------------------|-----------------------|----------------------------------------------------------------------|
| Beampipe inner | 0–17.8              | ±6000                 | Vacuum, tapered bore                                                  |
| Beampipe outer | 1–22.81 mm thick    | ±132.3                | Be wall, tapered                                                      |
| Nozzle         | 1–550 (z-dep.)      | ±60 to ±5950          | W core + BCH₂ + W cladding, 3 disk gaps per side                      |
| Forward Tracker| 25–542              | 5450–5852 (and −side) | 3 disks × 2 layers, Si in nozzle gaps                                 |
| Vertex         | 28–115              | 380                   | 5 barrel + 4 endcap layers                                            |
| Inner Tracker  | 61–580              | 2306                  | 3 barrel layers (164/354/554), 7 endcap disks                         |
| Outer Tracker  | 580–1499.9          | 2306                  | 3 barrel layers (819/1153/1486), 4 endcap disks                       |
| ECalBarrel     | 1690–1960           | 2210                  | CRILIN 6 × (40 mm PbF₂ + 1 mm Si + 3 mm siPCB + 1 mm Air)             |
| ECalEndcap     | 310–1960            | 2307–2577             | same CRILIN stack                                                     |
| Solenoid tank  | 2055–2862           | 2509                  | 80 mm steel walls, 423 mm Al coil at r=2459                           |
| HCalBarrel     | 2902–4756           | 2509                  | 70 × (20 mm Fe + 3 mm polystyrene + readout)                          |
| HCalEndcap     | 320/499–4756        | 2579–4434             | 70 layers, stepped inner cutout (5 steps)                             |
| YokeBarrel     | 4806–6800           | 4444                  | 7 × (RPC + RPC + 244 mm Air)                                          |
| YokeEndcap     | 493/622–6800        | 4444–5903             | 6 layers (197 mm Air + 2 × RPC)                                       |

**Field:** +5 T inside r < HCalBarrel_inner_radius (2902 mm) and |z| < 2509 mm;
−2 T return in 2902 < r < 4756 mm, |z| < 2509 mm; B=0 elsewhere.


## Changes since MuSIC_v2

1. **Forward Tracker** added (`ForwardTracker.xml`): 3 Si disks per side in
   nozzle gaps, each implemented as two `VertexEndcap_o1_v07` layers (Ring A /
   Ring B) with 8 modules per ring, in azimuthal brick-wall pattern.
2. **Nozzle v2** (`Nozzle_10deg_v2.xml`): segmented into 4 z-segments per side
   with 5 cm gaps at |z| = 542.5–547.5, 562.5–567.5, 582.5–587.5 cm to host
   the Forward Tracker disks. Radii at gap boundaries computed by linear
   interpolation along the original cone.
3. **`text_description.md`** has been regenerated and is now broadly
   consistent with the XML (resolves issue 1 of the v2 review).
4. **Beampipe**: `BeampipeOuter` now uses `material="Beryllium"` (resolves
   issue 2 of the v2 review).
5. **Yoke iron removed**: the 24.4 cm "iron" slab in `YokeBarrel` is now Air;
   `YokeEndcap` uses 19.7 cm Air slabs. This is intentional per
   `text_description.md` — the field return is approximated by the HCAL barrel
   iron and the yoke retains only its RPC chambers in geometric positions.
6. **Solenoid field zmax** now uses `HCalBarrel_half_length` (2509 mm) instead
   of `Solenoid_Coil_half_length` (2280 mm) — OuterTracker is now fully inside
   the field volume (resolves issue 6 of the v2 review).


## Issues found in v3

### 1. Stale title/comment/author in `MuSIC_v3.xml`

```xml
<info name="MuSIC_v3"
      title="MuSIC geometry v2 with 5 layers in VXD barrel"
      author="Davide Zuliani"
      ...>
    <comment>MuSIC v2 geometry with 5 layers in VXD barrel</comment>
</info>
```

The `name` field has been bumped to `MuSIC_v3` but `title` and the `<comment>`
still say "v2". With the Forward Tracker added and several material/field
changes, the v3 name now means something distinct — the title should reflect
that. The README.md is also stale (refers to `MuSIC_v2.xml` as the compact
file in the v3 directory).

### 2. `BoronOxide` material is missing boron (`materials.xml:203–206`)

```xml
<material name="BoronOxide">
    <D type="density" value="2.46" unit="g/cm3"/>
    <composite n="3" ref="O"/>
</material>
```

B₂O₃ has 2 B + 3 O; this definition has only O. PyrexGlass uses 13%
BoronOxide by mass (`materials.xml:222–226`), so the simulated PyrexGlass is
effectively SiO₂ + 13% extra "oxygen", missing the boron entirely. RPC
chambers (used in both Yoke barrel/endcap, every single layer) are built
from PyrexGlass, so this affects the yoke material budget. Fix:

```xml
<material name="BoronOxide">
    <D type="density" value="2.46" unit="g/cm3"/>
    <composite n="2" ref="B"/>
    <composite n="3" ref="O"/>
</material>
```

(This bug is inherited from MuSIC_v2 and was not flagged in the v2 review.)

### 3. Three mixtures have fraction sums ≠ 1.0 (runtime warnings)

dd4hep emits at load time:
- `GroundOrHVMix`: sum = 0.998
- `siPCBMix`: sum = 0.997
- `PCB`: sum = 0.998

ROOT renormalises silently, so the simulation runs, but the composition is
ambiguous — better to fix the typos. Concretely:

- `GroundOrHVMix` (`materials.xml:125–132`): 0.003+0.094+0.010+0.028+0.863 = 0.998. Bump one to 0.864 or 0.030.
- `siPCBMix` (`materials.xml:134–142`): 0.014+0.083+0.065+0.003+0.014+0.818 = 0.997. Bump Cu to 0.821 or one of the small fractions.
- `PCB` (`materials.xml:188–195`): 0.180+0.405+0.278+0.068+0.067 = 0.998. Bump O to 0.407 or Si to 0.182.

`GroundOrHVMix` is defined but does **not** appear in any v3 compact file —
it can simply be removed.

### 4. Unused / dead constants and materials

These declarations are present but never referenced:

| Symbol                 | Defined in                                |
|------------------------|-------------------------------------------|
| `NozzleCut_halfgap`    | `Nozzle_10deg_v2.xml:5`                   |
| `tracker_region_zmax`  | `MuSIC_v3.xml:138`                        |
| `tracker_region_rmax`  | `MuSIC_v3.xml:139`                        |
| `DetID_NOTUSED`        | `MuSIC_v3.xml:25`                         |
| `BCH2` material        | `materials.xml:144–152` (Nozzle uses `BCH2_composite`)|
| `TungstenDens24`       | `materials.xml:154–159`                   |
| `Tungsten` (pure W)    | `materials.xml:161–166` (Nozzle uses `Tungsten_light`)|
| `Carbon` (ρ = 2)       | `materials.xml:52–55`                     |
| `Iron`                 | `materials.xml:36–39` (HCal uses `Steel235`)|
| `beam` material        | `materials.xml:28–34`                     |

All harmless; cleanup is cosmetic.

### 5. BCH₂ has two definitions with inconsistent formulae

`materials.xml:144–152` defines `BCH2` by **mass fractions**:
H 11.6%, B 5%, C 61.2%, O 22.2% — this is borated polyethylene with the boron
content listed by mass, plus a 22% oxygen content that BCH₂ should not have.

`materials.xml:168–173` defines `BCH2_composite` by **atomic composites**:
H₂ B₁ C₁ — i.e. (BCH₂)ₙ — pure borated polyethylene by stoichiometry.

The Nozzle uses `BCH2_composite`. If `BCH2` is meant to model a different
borated-poly variant (e.g. with oxide content for fire-retardant grades), the
naming is misleading. If not, delete `BCH2` and keep `BCH2_composite`.

### 6. Vertex barrel uses non-contiguous layer IDs … now made contiguous

Re-checked from v2 review: in v3 `Vertex_o2_v06_01.xml`, the 5 barrel layers
now use IDs **0, 1, 2, 3, 4** (lines 89, 94, 99, 104, 109). The earlier
non-contiguous numbering (0, 2, 4, 6, 8) has been fixed. ✓

### 7. Yoke second RPC chamber still not flagged sensitive

`YokeBarrel_o1_v01_01.xml:51` and the equivalent line in every layer of
`YokeEndcap_o1_v01_01.xml` (e.g. lines 56, 74, 92, …) declare the second RPC
gas slice without `sensitive="yes"`. `text_description.md` (lines 213–219)
now documents that this is **intentional**, matching the MuColl/MAIA
convention, and warns against changing it without consulting the muon-system
contact. ✓ (now documented)

### 8. `Solenoid_o1_v01_01.xml` endcap is only half the wall thickness, leaves an axial gap

```xml
<detector name="SolenoidEndcaps" type="DD4hep_DiskTracker" ... reflect="true">
    <layer id="1"
           inner_z="Solenoid_half_length-SolenoidVacuumTank_thickness"
           inner_r="..." outer_r="...">
        <slice material="Steel235" thickness="SolenoidVacuumTank_thickness/2.0" />
    </layer>
</detector>
```

The endcap disk extends from `z = Solenoid_half_length - 80 = 2429 mm` to
`z = 2429 + 40 = 2469 mm` (thickness = 80/2 = 40 mm). The cylindrical wall
goes out to `z = Solenoid_half_length = 2509 mm`. There is a 40 mm-long
ring at 2469 < z < 2509 mm where the tank has cylindrical wall but no endcap
closure. Probably inherited intent, but worth confirming — either the
endcap should be the full 80 mm thick (closing flush with the cylindrical
wall) or the cylindrical wall should be shortened.

### 9. `dd4hep::CheckOverlaps` reports 72 overlaps

Running `gGeoManager->CheckOverlaps(0.001)` on the full MuSIC_v3 geometry
yields 72 overlaps. They group into four classes:

**9a. Vertex envelope vs first Nozzle W segment (2.29 mm, worst overlap)**

```
ov00006: world_volume/NozzleW_right_2 overlapping world_volume/Vertex_27  ovlp=2.29116
ov00007: world_volume/NozzleW_left_6  overlapping world_volume/Vertex_27  ovlp=2.29116
```

(`NozzleW_right_2` is ROOT's name for the first NozzleW_right detector
placement — index 2 in world_volume.) The Vertex SubdetectorAssembly's
envelope is a plain `Tube` with `rmin = Vertex_inner_radius + env_safety =
28.1 mm` and `dz = 380 mm`. The first NozzleW segment is a cone from
(z=60 mm, rmax=10 mm) to (z=5425 mm, rmax=390 mm) — at z=163 mm its outer
radius equals the Vertex envelope's rmin, and from there on it pokes into
the envelope, reaching up to r≈66 mm at z=380 mm. The intrusion is into
the **Air** of the Vertex envelope (no overlap with actual sensitive
silicon), but `CheckOverlaps` still flags it. Fix: replace the Tube
envelope with a Polycone whose inner profile follows the nozzle outer
profile, or simply remove the explicit envelope and let the assembly
compute its own.

**9b. HCal endcap layers extruding their envelope (0.16–0.69 mm)**

```
ov00000: endcap extruded by: endcap/layerType4_4  ovlp=0.694319
ov00001: endcap extruded by: endcap/layerType3_3  ovlp=0.561877
ov00002: endcap extruded by: endcap/layerType2_2  ovlp=0.429435
ov00003: endcap extruded by: endcap/layerType5_5  ovlp=0.342211
ov00004: endcap extruded by: endcap/layerType1_1  ovlp=0.296992
ov00005: endcap extruded by: endcap/layerType0_0  ovlp=0.16455
```

The five stepped HCal endcap layer blocks each extrude their envelope by up
to 0.7 mm. The HCalEndcap envelope has `env_safety = 0.1 mm` and the layer
total thickness is exactly `(max_z − min_z) / 70 = 26.5 mm` × 70 layers =
1855 mm, fitted to the available z space with no safety margin between the
last layer's outer face and `HCalEndcap_max_z`. Add ~1 mm safety to
`HCalEndcap_max_z` or reduce the layer count by one.

**9c. Beampipe ↔ Nozzle W cone interfaces (1 µm to 19 µm)**

```
ov00008: BeampipeOuter_1 overlapping NozzleW_left_6   ovlp=0.0189592   (19 µm)
ov00009: BeampipeOuter_1 overlapping NozzleW_right_2  ovlp=0.0189592   (19 µm)
ov00070: BeampipeInner_0 overlapping NozzleW_right_2  ovlp=0.0010101   (1 µm)
ov00071: BeampipeInner_0 overlapping NozzleW_left_6   ovlp=0.0010101   (1 µm)
```

The 19 µm Be↔W overlap matches the analytic prediction (slope mismatch
between the W outer cone, parameterised from z=6 cm to z=15 cm, and the Be
inner cone, parameterised from z=6.25 to z=13.23 cm). Fix: align the Be
inner rmin at each zplane to the W outer rmax computed at the same z.

**9d. Inner-Tracker endcap modules vs IT outer barrel support (~16 µm × 52)**

52 entries of the form

```
ovNNN: InnerTrackers/.../layer_pos4_8/InnerTrackerEndcapModule_9x9_Out_NN
       overlapping InnerTrackers/.../InnerTrackerBarrelSupport_layer3_3
       ovlp=0.0163496
```

All ~16 µm. The `InnerTrackerBarrelSupport` layer "3" here corresponds to
the XML's `<layer id="4">` (the outermost barrel-support cylinder at
`inner_r = InnerTracker_outer_radius − 1 cm = 570 mm`, `outer_z =
InnerTracker_half_length = 2306 mm`). IT endcap layer 4 modules at
z = ±1741 mm — specifically the outermost 9×9_Out modules — clip this
cylinder at ~16 µm depth. Likely the 9×9 module bounding-box outer edge
just barely touches the 570 mm cylinder. Small enough to ignore for
physics, but easy to remove by either nudging the OT barrel support
inner_r outwards by 0.1 mm or reducing the 9×9 module width slightly.

**9e. Vertex endcap module-to-module overlaps (~7 µm × 8)**

```
ov00062-69: SiVertexEndcapModule{1..4}_0_0 overlapping SiVertexEndcapModule{1..4}_0_1
            ovlp=0.00746578  (7.5 µm)
```

Adjacent azimuthal modules in each vertex endcap ring overlap by 7.5 µm. This
is built into the design (`VertexEndcapOverlap = 1 mm` is added to the
trapezoid x1/x2 expressions in `Vertex_o2_v06_01.xml:129–149` to guarantee
azimuthal coverage). The reported 7.5 µm is just the depth at which the
trapezoid corners actually clip — the design overlap of 1 mm is *radial*,
not bulk. Probably intentional; verify no double-counting in reconstruction.

### 10. `tracker_region_*` constants suggest a tracker region was planned but never installed

`MuSIC_v3.xml:138–139` defines `tracker_region_zmax` and `tracker_region_rmax`
but no `<region name="TrackerRegion">` exists in the `<regions>` block. The
trackers each have their own per-detector regions
(`VertexBarrelRegion`, `InnerTrackerEndcapRegion`, …), so the global tracker
region constants are vestigial. Either install the region or remove the
constants.

### 11. Module mass: text description undercount

`text_description.md` describes Inner-/Outer-tracker modules as "200 µm Si +
250 µm Al", but the actual module stacks
(`InnerTrackerBarrelModuleDown.xml`, `OuterTrackerBarrelModuleUp/Down.xml`,
`TrackerDiskModuleIn/Out.xml`) also include substantial CarbonFiber, Epoxy,
Kapton, Rohacell, Allcomp_K9 and Water layers (≈1–4 mm total material
per module). FLUKA users who follow only the Si+Al line will underestimate
the tracker material budget by an order of magnitude. Worth noting
explicitly in `text_description.md` that the figures listed there are the
sensor-stack only and that the full module thicknesses are 1.66 mm (IT
barrel), 4.62 mm (OT barrel, with rohacell core), and similar for endcap
modules.

### 12. ForwardTracker shares `NozzleRegion`

`ForwardTracker.xml:49` puts the disks in the same region as the nozzle.
This is fine but means that any production cuts or step limits associated
with NozzleRegion (currently none — region body is empty) will apply equally
to the silicon disks. Worth verifying whether a dedicated
`ForwardTrackerRegion` would be more appropriate, especially once cuts get
attached.

### 13. ForwardTracker disk thickness is 0.33 mm, fits comfortably in 5 cm gap

Verified: at each disk, ring A is at z=zstart and ring B at z=zstart+2 mm,
each with 0.33 mm total Si stack thickness; both lie well inside the
respective 5 cm gap. No clearance issue.


## Sanity checks that passed

- Envelope radial ordering is monotonic: 28 < 61 < 115 < 580 < 1500 < 1690 <
  1960 < 2055 < 2862 < 2902 < 4756 < 4806 < 6800. No gaps inverted, no
  envelopes nested out of order.
- All `_xx_v##` plugin types referenced by the XML are present in the local
  `libk4geo.so` build (`VertexEndcap_o1_v07`, `TrackerBarrel_o1_v05`,
  `TrackerEndcap_o2_v07`, `GenericCalBarrel_o1_v01`, etc.).
- Forward Tracker disk midpoint radii (512.197, 527.318, 542.439 mm) match
  the linearly-interpolated nozzle-cladding outer radius at z = 545, 565,
  585 cm to better than 1 mm — disks reach exactly to the cladding outer
  edge.
- Solenoid field zmax = HCalBarrel_half_length = 2509 mm covers both
  InnerTracker_half_length (2306 mm) and OuterTracker_half_length (2306 mm).
  Tracker fully magnetised.
- All `_test_*` and backup files have been removed from the directory.
- All elements referenced in `materials.xml` exist in `elements.xml`.


## Suggested actions (priority order)

1. **High** — Fix `BoronOxide` to actually contain boron (`materials.xml:204`).
   Affects every RPC chamber (PyrexGlass plates) via Pyrex composition.
2. **High** — Fix the Vertex envelope shape (issue 9a, 2.29 mm overlap with
   the nozzle). Replace the `Tube` envelope with a `Polycone` whose inner
   radius widens to clear the nozzle, or drop the explicit envelope so the
   `SubdetectorAssembly` computes its own bounds.
3. **Medium** — Fix the three mixture fraction sums (`GroundOrHVMix`, `siPCBMix`,
   `PCB`) to total 1.0. Remove `GroundOrHVMix` since it's unused.
4. **Medium** — Add ~1 mm safety to `HCalEndcap_max_z` (or trim one HCal endcap
   layer) to resolve the 0.16–0.69 mm layer extrusions (issue 9b).
5. **Medium** — Match the Be inner cone to the W outer cone in
   `Beampipe_o1_v01_02.xml` to clear the 19 µm overlap (issue 9c).
6. **Low** — Update `MuSIC_v3.xml` `title` and `<comment>` and `README.md` to
   reference v3.
7. **Low** — Delete dead constants and unused materials (issue 4), and
   clarify the `BCH2` vs `BCH2_composite` duplication (issue 5).
8. **Low** — Document the sensor-only convention in `text_description.md`'s
   per-module mass figures, or list the full stack (issue 11).
9. **Low** — Nudge the IT outer barrel support inner_r outwards by 0.1 mm to
   clear the 52 × 16 µm overlaps with IT endcap 9×9 modules (issue 9d).
10. **Optional** — Decide whether to introduce a `ForwardTrackerRegion` rather
    than reusing `NozzleRegion`.
11. **Optional** — Decide whether the Solenoid endcap should close flush with
    the cylindrical wall (currently 40 mm short).
