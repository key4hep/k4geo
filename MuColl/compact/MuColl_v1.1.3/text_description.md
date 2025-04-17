## Textual description of the geometry
A simplified overview to easily recreate the main elements of the geometry in FLUKA for radiation studies

### Vertex Barrel layers
Silicon thickness: 190e-3 mm

|ID | R [mm]    | Zmax [mm] |
|---|-----------|-----------|
|0  | 30        | 65        |
|1  | 32        | 65        |
|2  | 51        | 65        |
|3  | 53        | 65        |
|4  | 74        | 65        |
|5  | 76        | 65        |
|6  | 102       | 65        |
|7  | 104       | 65        |


### Vertex Endcap disks
Silicon thickness: 330e-3 mm

|ID | Z [mm]    | Rmin [mm] | Rmax [mm] |
|---|-----------|-----------|-----------|
|0  | 80        | 25        | 112       |
|1  | 82        | 25        | 112       |
|2  | 120       | 31        | 112       |
|3  | 122       | 31        | 112       |
|4  | 200       | 38        | 112       |
|5  | 202       | 38        | 112       |
|6  | 280       | 53        | 112       |
|7  | 282       | 53        | 112       |


### Inner Barrel layers
Silicon thickness: 200e-3 mm
Aluminium thickness: 250e-3 mm

|ID | R [mm]    | Zmax [mm] |
|---|-----------|-----------|
|0  | 127       | 481.6     |
|1  | 340       | 481.6     |
|2  | 554       | 692.3     |


### Inner Endcap disks
Silicon thickness: 200e-3 mm
Aluminium thickness: 250e-3 mm

|ID | Z [mm]    | Rmin [mm] | Rmax [mm] |
|---|-----------|-----------|-----------|
|0  | 524       | 95        | 405       |
|1  | 808       | 147       | 555       |
|2  | 1093      | 190       | 555       |
|3  | 1377      | 212       | 555       |
|4  | 1661      | 237       | 555       |
|5  | 1946      | 264       | 555       |
|6  | 2190      | 284       | 555       |


### Outer Barrel layers
Silicon thickness: 200e-3 mm
Aluminium thickness: 250e-3 mm

|ID | R [mm]    | Zmax [mm] |
|---|-----------|-----------|
|0  | 819       | 1249.15   |
|1  | 1153      | 1249.15   |
|2  | 1486      | 1249.15   |


### Outer Endcap disks
Silicon thickness: 200e-3 mm
Aluminium thickness: 250e-3 mm

|ID | Z [mm]    | Rmin [mm] | Rmax [mm] |
|---|-----------|-----------|-----------|
|0  | 1310      | 617.5     | 1430.2    |
|1  | 1617      | 617.5     | 1430.2    |
|2  | 1883      | 617.5     | 1430.2    |
|3  | 2190      | 617.5     | 1430.2    |


### Calorimeter Barrel
Polyhedron: 12 sides

|ID     | Rmin  [mm]| Rmax [mm] | Zmax [mm] | Radiator [mm] | Sensor [mm]           |
|-------|-----------|-----------|-----------|---------------|-----------------------|
|ECAL   | 1500      | 1702      | 2210      | 40 x 1.9 (W)  | 40 x 0.5 (Si)         |
|HCAL   | 1740      | 3330      | 2210      | 60 x 20 (Fe)  | 60 x 3 (Polystyrene)  |

### Calorimeter Endcap
Polyhedron: 12 sides

|ID     | Rmin  [mm]| Rmax [mm] | Zmin [mm] | Zmax [mm] | Radiator [mm] | Sensor [mm]           |
|-------|-----------|-----------|-----------|-----------|---------------|-----------------------|
|ECAL   | 310       | 1700      | 2307      | 2509      | 40 x 1.9 (W)  | 40 x 0.5 (Si)         |
|HCAL   | 307       | 3246      | 2539      | 4129      | 60 x 20 (Fe)  | 60 x 3 (Polystyrene)  |
|HCAL_R | 1738      | 3246      | 2353.5    | 2539      | 7 x 20 (Fe)   | 7 x 3 (Polystyrene)   |


### Solenoid
Magnetic field: 3.57 T (inside);  -1.34 T (outside);

|ID     | Rmin  [mm]| Rmax [mm] | Zmax [mm] | Thickness [mm] |
|-------|-----------|-----------|-----------|----------------|
|Tank   | 3483      | 4290      | 4129      | 40 (Fe)        |
|Coil   | 3649      | 3993      | 3900      | 344 (Al)       |


### Return Yoke
Polyhedron: 12 sides

|ID     | Rmin  [mm]| Rmax [mm] | Zmin [mm] | Zmax [mm] | Thickness [mm]    |
|-------|-----------|-----------|-----------|-----------|-------------------|
|Barrel | 4461      | 6450      | -         | 4179      | 7 x 244 (Fe)      |
|Endcap | 575       | 6450      | 4179      | 5638      | 6 x 197 + 98 (Fe) |

