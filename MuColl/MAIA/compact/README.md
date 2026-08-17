# MAIA detector models - Overview

The following MAIA detector models are available in k4geo ( current production models **MAIA_v0** )

| Model         |  Description               | Hcal   |  Ecal   | Muon     | Status            |
| ------------- | ---------------------------|--------|---------|----------|-------------------|
| MAIA_v0	    | Initial MAIA design        | sci    | Si      | RPC      |  validated        |
| ------------- | ---------------------------|--------|---------|----------|-------------------|

## Details

### MAIA_v0
- Initial detector model
- Used for ESPPU 2026 update and 2502.00181
- Compact file: [./MAIA_v0/MAIA_v0.xml](./MAIA_v0/MAIA_v0.xml)

## Magnetic field

MAIA_v0 comes with two interchangeable descriptions of the solenoid field. The
geometry is identical in both; only the `<fields>` include differs.

| Compact file | Field | When to use |
| ------------ | ----- | ----------- |
| [./MAIA_v0/MAIA_v0.xml](./MAIA_v0/MAIA_v0.xml) | [Field_Solenoid_Analytic_5T.xml](./MAIA_v0/Field_Solenoid_Analytic_5T.xml) | Default. DD4hep's idealised `solenoid`: uniform 5 T inside the coil, uniform -0.937 T in the return annulus, no Br, hard cut at `zmax`. |
| [./MAIA_v0/MAIA_v0_fieldmap.xml](./MAIA_v0/MAIA_v0_fieldmap.xml) | [Field_Solenoid_Map_5T.xml](./MAIA_v0/Field_Solenoid_Map_5T.xml) | Computed (Br, Bz) map with a realistic fringe field and a non-zero radial component. |

`MAIA_v0_fieldmap.xml` is a copy of `MAIA_v0.xml` with only that one include
changed. Edits to either must be mirrored in the other, which the
`t_MAIA_v0_fieldmap_in_sync` test enforces.
