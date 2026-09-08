## Textual description of the geometry
A simplified overview to easily recreate the main elements of the geometry in FLUKA for radiation studies.
Regenerated from the XML compact files on 2026-05-14.

All silicon "support" thicknesses below are the Silicon material assigned to passive
module/ladder elements in the XML (used as a stand-in for the carbon-fibre support
mass in some MuColl-derived models).


### Beampipe
Two coaxial volumes around the interaction point:

|Volume         | z range [cm]      | Rmin [cm]  | Rmax [cm]  | Material  |
|---------------|-------------------|------------|------------|-----------|
|Inner bore     | −600 to +600      | 0          | 0.3 → 1.78 (tapered) | Vacuum    |
|Be wall        | −13.23 to +13.23  | 1.0 → 2.281| 2.2 → 2.281| Beryllium |

The inner-bore taper follows the nozzle profile: rmax = 1 cm at the IP, narrows
to 0.3 cm at |z| = 100 cm, then widens back to 1.78 cm at |z| = 600 cm. Beyond
|z| > 13.23 cm there is no separate Be tube — the tungsten nozzle bore acts as
the beampipe wall.


### Nozzle (10° MAP-style)
Symmetric on both sides of the IP. Inner tungsten cone + borated-polyethylene
middle + tungsten cladding. Each of the three radial components is interrupted
by three 5 cm wide air gaps at |z| = 542.5–547.5, 562.5–567.5, and 582.5–587.5
cm to accommodate the Forward Tracker silicon disks; the radii at each gap
boundary are computed by exact linear interpolation along the nozzle cone.

| Component       | z range [cm]                          | Material       |
|-----------------|---------------------------------------|----------------|
| Tungsten core   | ±6 to ±595 (4 segments per side)      | Tungsten_light |
| BCH2 layer      | ±100 to ±595 (4 segments per side)    | BCH2_composite |
| Tungsten clad   | ±100 to ±595 (4 segments per side)    | Tungsten_light |

Tungsten core inner radius starts at 1 cm (z = ±6 cm) and tapers down to 0.3 cm
at |z| = 100 cm, then opens to 1.78 cm at |z| = 595 cm (matching the beampipe
inner-bore profile).

The outer radius of the tungsten cladding at the three disk positions (midpoint
of each gap) is approximately 512, 527, and 542 mm — these set the outer radius
of the corresponding Forward Tracker disks.


### Forward Tracker (silicon disks in nozzle gaps)
Three double-layer silicon disk stations placed in the 1 cm nozzle gaps on each
side. Each station consists of two staggered rings (A and B) registered as
separate detector layers so that every layer contains exactly one ring
(VertexEndcap_o1_v07 constraint). Each ring carries 8 trapezoidal modules; the
two rings together tile the full azimuth with 16 modules covering 22.5° each, in
a brick-wall arrangement (ring A covers even slots, ring B covers odd slots, offset
by 2 mm in z).

Module stack per module: 280 µm Si (passive support) + 50 µm Si (sensitive).
Total silicon thickness per layer: 330 µm.

| Disk | |z| centre [mm] | Rmin [mm] | Rmax [mm] | Layer IDs (A/B) |
|------|-----------------|-----------|-----------|-----------------|
| 0    | 5450            | 25        | 512.197   | 0 / 1           |
| 1    | 5650            | 25        | 527.318   | 2 / 3           |
| 2    | 5850            | 25        | 542.439   | 4 / 5           |

Each disk has 8 modules in ring A (phi0 = π/16) and 8 modules in ring B
(phi0 = 3π/16, zstart + 2 mm). The detector is reflected to cover both ±z sides
(`reflect="true"`), giving 6 layers × 2 sides = 12 active layers in total.
Readout: `ForwardTrackerCollection` (same bit-field as GlobalTrackerReadoutID).


### Vertex Barrel layers (5 single layers)
Silicon thickness: 190e-3 mm (50e-3 sensor + 140e-3 support)

|ID | R [mm]    | Zmax [mm] |
|---|-----------|-----------|
|0  | 29        | 130       |
|1  | 40        | 130       |
|2  | 50        | 130       |
|3  | 73        | 130       |
|4  | 101       | 130       |


### Vertex Endcap disks (4 single disks per side)
Silicon thickness: 330e-3 mm (50e-3 sensor + 280e-3 support)

|ID | Z [mm]    | Rmin [mm] | Rmax [mm] |
|---|-----------|-----------|-----------|
|0  | 180       | 34.7      | 112       |
|1  | 230       | 43.5      | 112       |
|2  | 298       | 55.5      | 112       |
|3  | 366       | 67.5      | 112       |


### Inner Tracker Barrel layers
Silicon thickness: 200e-3 mm (100e-3 sensor + 100e-3 r/o ASIC)
Aluminium thickness: 250e-3 mm (200e-3 power-bus conductor + 50e-3 FPC metal)

|ID | R [mm]    | Zmax [mm] |
|---|-----------|-----------|
|0  | 164       | 481.6     |
|1  | 354       | 481.6     |
|2  | 554       | 692.3     |


### Inner Tracker Endcap disks (7 per side)
Silicon thickness: 200e-3 mm
Aluminium thickness: 250e-3 mm

|ID | Z [mm]    | Rmin [mm] | Rmax [mm] |
|---|-----------|-----------|-----------|
|0  | 604       | 116       | 405       |
|1  | 888       | 166       | 555       |
|2  | 1173      | 201       | 555       |
|3  | 1457      | 225       | 555       |
|4  | 1741      | 249       | 555       |
|5  | 1946      | 266       | 555       |
|6  | 2190      | 287       | 555       |


### Outer Tracker Barrel layers
Silicon thickness: 200e-3 mm
Aluminium thickness: 250e-3 mm

|ID | R [mm]    | Zmax [mm] |
|---|-----------|-----------|
|0  | 819       | 1264.2    |
|1  | 1153      | 1264.2    |
|2  | 1486      | 1264.2    |


### Outer Tracker Endcap disks (4 per side)
Silicon thickness: 200e-3 mm
Aluminium thickness: 250e-3 mm

|ID | Z [mm]    | Rmin [mm] | Rmax [mm] |
|---|-----------|-----------|-----------|
|0  | 1410      | 617.5     | 1430.2    |
|1  | 1667      | 617.5     | 1430.2    |
|2  | 1933      | 617.5     | 1430.2    |
|3  | 2190      | 617.5     | 1430.2    |


### Calorimeter Barrel
Polyhedron: 12 sides (Rmin is the inscribed-circle radius).

|ID    | Rmin [mm] | Rmax [mm] | Zmax [mm] | Layers | Sensor (per layer)             | Other material (per layer)               |
|------|-----------|-----------|-----------|--------|--------------------------------|------------------------------------------|
|ECAL  | 1690      | 1960      | 2210      | 6      | 40 mm PbF2 (LeadDifluoride)      | 1 mm Si + 3 mm siPCBMix + 1 mm Air       |
|HCAL  | 2902      | 4756      | 2509      | 70     | 3 mm Polystyrene               | 20 mm Steel235 + 0.1 mm Cu + 0.7 mm PCB + 0.5 mm Steel235 + 2.7 mm Air |

ECAL is a homogeneous PbF2 (CRILIN) calorimeter: the 40 mm crystal is the
sensitive medium and there is no tungsten radiator. Total ECAL radial depth
= 6 × 45 mm = 270 mm.
Total HCAL radial depth = 70 × 26.5 mm = 1855 mm.

The HCAL iron absorbers (20 mm Steel235 per layer + 0.5 mm trim slices) also
provide the magnetic-field return: ~1400 mm of iron per barrel stave.


### Calorimeter Endcap
Polyhedron: 12 sides.

|ID    | Rmin [mm]           | Rmax [mm] | Zmin [mm] | Zmax [mm] | Layers | Sensor (per layer)          | Other material (per layer)               |
|------|---------------------|-----------|-----------|-----------|--------|-----------------------------|------------------------------------------|
|ECAL  | 310                 | 1960      | 2307      | 2577      | 6      | 40 mm PbF2 (LeadDifluoride)   | 1 mm Si + 3 mm siPCBMix + 1 mm Air       |
|HCAL  | 320 → 499 (5 steps) | 4756      | 2579      | 4434      | 70     | 3 mm Polystyrene            | 20 mm Steel235 + 0.1 mm Cu + 0.7 mm PCB + 0.5 mm Steel235 + 2.7 mm Air |

HCAL endcap has a conical inner cutout to clear the nozzle: the inner radius
steps from 320 mm at the front face (z = 2579 mm) to 499 mm at the back face
(z = 4434 mm) in 5 equal radial steps (≈35.8 mm each). 14 layers occupy each
step.


### Solenoid
Vacuum tank with steel walls and a single aluminium coil:

|ID                  | Rmin [mm] | Rmax [mm] | Zmax [mm] | Material               |
|--------------------|-----------|-----------|-----------|------------------------|
|Tank inner wall     | 2055      | 2135      | 2509      | 80 mm Steel235         |
|Coil                | 2247.5    | 2670.5    | 2280      | 423 mm Aluminium       |
|Tank outer wall     | 2782      | 2862      | 2509      | 80 mm Steel235         |

The volume between coil and tank walls is filled with Vacuum.

**Magnetic field** (modelled as an axial DD4hep `solenoid` field with two
cylindrical regions, no field outside):

| Region                                                          | B    |
|-----------------------------------------------------------------|------|
| r < 2902 mm and \|z\| < 2509 mm (tracker, ECAL, solenoid bore)  | +5 T |
| 2902 < r < 4756 mm and \|z\| < 2509 mm (HCAL barrel iron)       | −2 T |
| Elsewhere (HCAL endcap, yoke volume, beyond tracker z)          | 0    |

Because the iron has been removed from the yoke, the magnetic-flux return
is approximated by the HCAL barrel iron, which is treated as the return-yoke
material. The HCAL endcap and the muon-yoke volumes carry no field in this
simplified model.


### Muon System (formerly Return Yoke)
Polyhedron: 12 sides. **No iron absorber** in this version — the iron mass has
been removed from the simulation; the RPC chambers are kept in their original
geometric positions by replacing each former iron slab with an equal-thickness
air slice.

|ID     | Rmin [mm]    | Rmax [mm] | Zmin [mm] | Zmax [mm] | Layers | Per-layer structure                        |
|-------|--------------|-----------|-----------|-----------|--------|--------------------------------------------|
|Barrel | 4806         | 6800      | —         | 4444      | 7      | 2 × RPC chamber + air gap + 244 mm of air  |
|Endcap | 493 → 622    | 6800      | 4444      | 5903      | 6      | 197 mm of air + 2 × RPC chamber            |

The endcap layers have a stepped inner radius following the nozzle
(rmin increases by ≈20 mm per layer from 493 to 613 mm).

Each RPC chamber is built as:
PyrexGlass 2 mm + RPCGas 2 mm + PyrexGlass 2 mm, sandwiched between 1 mm Al
covers with 3.5 mm air gaps (≈15 mm per chamber, two stacked back-to-back).
Only the first RPC chamber of each pair carries the `sensitive="yes"` flag;
the second chamber is geometrically present but does not write hits to the
readout collection. This matches the convention used in all MuColl/MAIA-derived
geometries (consult the muon-system contact before changing it).
