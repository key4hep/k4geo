# k4geo (Key4hep Detector Geometries)

[![DOI](https://zenodo.org/badge/60772160.svg)](https://doi.org/10.5281/zenodo.596333)
[![Key4hep build](https://github.com/key4hep/k4geo/actions/workflows/key4hep-build.yaml/badge.svg)](https://github.com/key4hep/k4geo/actions/workflows/key4hep-build.yaml)

Implementation of detector models in DD4hep, for use within the [Key4hep](https://github.com/key4hep) software stack.

k4geo is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Requirements

To build `k4geo` the following packages are required:
- DD4hep
- Geant4
- ROOT

They are available for instance via the [Key4hep software stack](https://key4hep.github.io/key4hep-doc/), available on machines which have CVMFS mounted (e.g. lxplus).

## Download and Installation

```bash
git clone https://github.com/key4hep/k4geo.git
cd k4geo
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j4
cd ..
```
<details>

<summary>For people using the key4hep stack</summary>

The required packages for installation are available via the [Key4hep software stack](https://key4hep.github.io/key4hep-doc/), available on machines which have CVMFS mounted (e.g. lxplus). <br>

### Source the Key4hep stack

For development and to pick up the latest changes from the rest of the stack, source the nightly build:

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```
Alternatively, source the stable `key4hep` release:
```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh
```
Be aware that the version of `k4geo` you cloned may not be compatible with the stable release.

Sourcing the key4hep stack gives you access to the central version of `k4geo` via the `$K4GEO` environment variable.
When using the nightly stack, this will point to a `k4geo` version with all the latest updates (up until the previous night), whereas using the stable stack will point to the latest release version. <br>
The local installation of `k4geo` is thus only required if you plan to make local changes for development. <br>
See section [After Installation](#after-installation) for instructions how to use your local `k4geo` with the rest of the key4hep stack. 

### Download and Installation

Same steps as above

### After Installation

To make the `$K4GEO` variable point to your local version of `k4geo` use the command

```bash
k4_local_repo
```

in the top-level `k4geo` directory (where you should end up being after following the installation instructions above). 
This command adds your local `k4geo` path to the environment variable `$K4GEO` such that other processes will use your version rather than the central one.

</details>

## Simulating some events

To simulate events for a given detector use the `ddsim` command:

```bash
ddsim -G -N 10 path/to/DETECTOR.xml
``` 
`-G` enables the particle gun (10 GeV electrons per default) and `-N 10` sets the number of events to 10.

Note that this uses the `ddsim` default settings which are likely not appropriate for your detector.
For detector specific configuration for running ddsim, see the following repositories:

### FCCee detectors

- `ALLEGRO`: [`${FCCCONFIG}/FullSim/ALLEGRO/`](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/ALLEGRO)
- `CLD`: [`${CLDCONFIG}/share/CLDConfig/`](https://github.com/key4hep/CLDConfig)
- `IDEA`: [`${FCCCONFIG}/FullSim/IDEA/`](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/IDEA)
- `ILD@FCCee`: [`${FCCCONFIG}/FullSim/ILD_FCCee/`](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/ILD_FCCee)

### Linear Collider detectors

- `CLIC`: [CLICPerformance/clicConfig](https://github.com/iLCSoft/CLICPerformance/tree/master/clicConfig)
- `ILD`: [ILDConfig](https://github.com/iLCSoft/ILDConfig)
- `SiD`: [SiDPerformance](https://github.com/iLCSoft/SiDPerformance)

### Muon Collider detectors

#### Simulation

The simulation workflow is the same for all muon collider detector concepts, kept centrally at
[MuonColliderSoft/mucoll-benchmarks](https://github.com/MuonColliderSoft/mucoll-benchmarks).

#### Reconstruction

- `MAIA`: [MAIAConfig](https://github.com/MuonColliderSoft/MAIAConfig)
- `MuSIC`: [MuSICConfig](https://github.com/MuonColliderSoft/MuSICConfig)
- `MuColl (3 TeV concept)`: [MuCollConfig](https://github.com/MuonColliderSoft/MuCollConfig)

## Visualisation

- `k4CEDViewer/CED`:
    - "CLIENT" - k4CEDViewer Codebase: [https://github.com/key4hep/k4CEDViewer](https://github.com/key4hep/k4CEDViewer)
    - "SERVER" - CED Codebase: [https://github.com/iLCSoft/CED](https://github.com/iLCSoft/CED)
- for FCC detectors: `Phoenix@FCC`
    - Codebase: [https://github.com/HEP-FCC/phoenix-at-fcc](https://github.com/HEP-FCC/phoenix-at-fcc)
    - Deployment: [https://hep-fcc.github.io/phoenix-at-fcc/#/](https://hep-fcc.github.io/phoenix-at-fcc/#/)

## License and Copyright

Copyright (C), k4geo Authors

k4geo is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
