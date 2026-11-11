# MAIA detector models - Overview

The following MAIA detector models are available in k4geo ( current production models **MAIA_v0** )

| Model         |  Description               | Hcal   |  Ecal   | Muon     | Status            |
| ------------- | ---------------------------|--------|---------|----------|-------------------|
| MAIA_v0	    | Initial MAIA design        | sci    | Si      | RPC      |  validated        |
| MAIA_v1	    | MAIA_v0 geometry, computed field map | sci    | Si      | RPC      |  development      |
| ------------- | ---------------------------|--------|---------|----------|-------------------|

## Details

### MAIA_v0
- Initial detector model
- Used for ESPPU 2026 update and 2502.00181
- Compact file: [./MAIA_v0/MAIA_v0.xml](./MAIA_v0/MAIA_v0.xml)
- Solenoid field: [Field_Solenoid_Analytic_5T.xml](./MAIA_v0/Field_Solenoid_Analytic_5T.xml) —
  DD4hep's idealised `solenoid`, a uniform 5 T inside the coil and a uniform
  -0.937 T in the return annulus, with no Br and a hard cut at `zmax`

### MAIA_v1
- Same geometry as MAIA_v0, included directly from `../MAIA_v0/`, so the two
  cannot diverge in geometry; only the field description differs
- Compact file: [./MAIA_v1/MAIA_v1.xml](./MAIA_v1/MAIA_v1.xml)
- Solenoid field: [Field_Solenoid_Map_5T.xml](./MAIA_v1/Field_Solenoid_Map_5T.xml) —
  a computed (Br, Bz) map with a realistic fringe field and a non-zero radial
  component, so track states in the calorimeter endcaps no longer sit on a
  5 T to 0 T cliff
- The map is `fieldmaps/MAIA_fieldMap_Solenoid5T_v1.root`, regenerated with
  [utils/make_maia_fieldmap.py](../../../utils/make_maia_fieldmap.py). Its
  methodology, ROOT file format and validation are documented in
  [./MAIA_v1/fieldmap_description.md](./MAIA_v1/fieldmap_description.md)
