Cleaned-up version of `MuColl_v1`:
- removed unused subdetectors and constants
- added reduced geometries for faster subdetector studies

Resolved 82 overlaps in the geometry:
- removed overlap margin in VTX Endcap Si sensors (was introduced by CLIC for tilted sensors)
- reduced radius of the 1st section of `BeampipeShell` support, away from IT Si sensors

Layout of sensitive surfaces should remain exactly as in `MuColl_v1`
