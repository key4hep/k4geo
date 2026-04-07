This is the FCChh configuration, copied from FCCDetectors as of March, 2026
and adjusted to work with k4geo.

Further notes:
  - Endcaps_coneCryo.xml and Forward_coneCryo.xml both defined the variable
    CryoThicknessInner with differing values: 100mm for the former and 20mm
    for the latter.  It appears that in previous configurations, the first
    one was what was used.  These variables have been adjusted to make them
    unique.
  - The simulation and reconstruction segmentations have been made consistent,
    adding additional CellID fields where needed.
