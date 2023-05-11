## Geometries of the Muon Collider detector

Current baseline version: **`MuColl_v1.1`**

| Geometry name             | Description |
|---------------------------|-------------|
|                           | **Baseline geometries** |
| `MuColl_v0`               | 1st version with `MuColl` name, based on `CLIC` detector, adapted to MDI designed by MAP |
| `MuColl_v1`               | Copy of `MuColl_v0` with fixed asymmetry in thickness of Tracker Endcap Support structures (frozen for SnowMass 2021 studies) |
| `MuColl_v1.1` +           | Cleaned-up version of `MuColl_v1` with resolved overlaps and better code layout |
|                           | **Alternative geometries** |
| `MuColl_v1.0.1`           | Copy of `MuColl_v1` with Crilin design of ECAL |
| `MuColl_v1.0.2`           | Copy of `MuColl_v1` with MPGD material in Yoke |
| `MuColl_v1.1.1`           | Copy of `MuColl_v1.1` with 500um of Si added to passive material thickness in VTX detector |
| `MuColl_v1.1.2`           | Copy of `MuColl_v1.1` with a stronger magnetic field: 5.0 T |
| `MuColl_v1.1.3`           | Copy of `MuColl_v1.1` with increased double-layer gap in VTX: 4.0 mm |


### Geometry versioning

Full name of a geometry version consists of two identifiers:
* detector name: `MuColl`;
* semantic version with up to 3 dot-separated digits, e.g. `v1.0.1`, `v1`.

> NOTE: Only the major version number should be used in the root `XML` file for easier use in configuration files, e.g. `MuColl_v1.1.1/MuColl_v1.xml`.

**Major version** (1-st digit) of a new geometry is increased only if it is meant to be a new baseline geometry and is significantly different from the previous baseline geometry.
Geometries with different major versions are not guaranteed to be compatible between each other.
Normally switching to a geometry with a higher major version requires re-simulation of all the Monte Carlo samples.

**Minor version** (2-nd digit) of a new geometry is increased only if it is meant to be a new baseline geometry and only contains minor changes or bug fixes that do not change the layout of sensitive elements.
Geometries differing only by their minor version number should be technically interchangeable, allowing to use a newer geometry for reconstruction of events simulated with an older geometry, without rerunning the GEANT4 step.

**Patch version** (3-rd digit) of a new geometry is increased if it is an experimental variation of the baseline geometry, which doesn't have to be backward compatible with any previous version.

> NOTE: Every new geometry version should be placed in a separate folder with its name, adding a corresponding entry in the table above, summarising briefly its purpose.
> A more detailed explanation of the geometry should be given in the dedicated `README.md` file inside the geometry folder.


### Deprecated versions

These are old geometries used in the early days of full-simulation studies, listed for historic reference.
`CLIC_o3_v14` geometry was used as a starting point with several adjustments gradually implemented to cope with the beam-induced background.

| Geometry name             | Description |
|---------------------------|-------------|
| `CLIC_o3_v14_mod1`        | Modified CLIC geometry with inner radii of all forward detectors increased to accomodate the MAP nozzles for sqrt(s) = 1.5TeV Muon Collider. Vertex Endcaps have no propeller structure. |
| `CLIC_o3_v14_mod2`        | Fixed ECAL Endcaps to fit tight around the nozzles |
| `CLIC_o3_v14_mod3`        | Vertex detector with shorter barrel and +1 layer/disk to be similar to the MAP design |
| `CLIC_o3_v14_mod4`        | Vertex Barrel segmented in 5 modules (L=26mm) along Z |
| `CLIC_o3_v14_mod5`        | Better Theta coverage in Vertex Endcap. MAP magnetic field |
