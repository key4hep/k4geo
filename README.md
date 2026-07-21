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
- LCIO
- ROOT

They are available for instance via the [Key4hep software stack](https://key4hep.github.io/key4hep-doc/), available on machines which have CVMFS mounted (e.g. lxplus). <br>
In the following instructions we will assume you have access to key4hep.

## Download and Installation

### Source the Key4hep stack

For development and to pick up the latest changes from the rest of the stack, source the nightly build:

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```
Alternatively, source the stable `key4hep` release.
Be aware that the version of `k4geo` you cloned may not be compatible with the stable release.

### Download

```bash
git clone https://github.com/key4hep/k4geo.git
cd k4geo
```

### Build

```bash
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j4
cd ..
```

### Set up the local repository in the environment

```bash
k4_local_repo
```

(in the top-level `k4geo` directory). 
This ensures that key4hep uses your local version of `k4geo` instead of the central one in CVMFS.

## Simulating some events

To simulate events for a given detector use the `ddsim` command:

```bash
ddsim -G -N 10 path/to/DETECTOR.xml
``` 
`-G` enables the particle gun (10 GeV electrons per default) and `-N 10` sets the number of events to 10.

Note that this uses the `ddsim` default settings which are likely not appropriate for your detector.
For detector specific configuration for running ddsim, see the following repositories:

### FCCee detectors

- `ALLEGRO`: [`${FCCCONFIG}/FCCee/FullSim/ALLEGRO/`](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/ALLEGRO)
- `IDEA`: [`${FCCCONFIG}/FCCee/FullSim/IDEA/`](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/IDEA)
- `ILD@FCCee`: [`${FCCCONFIG}/FCCee/FullSim/ILD/`](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/ILD)

### Linear Collider detectors

- `CLD`: [`${CLDCONFIG}/share/CLDConfig/`](https://github.com/key4hep/CLDConfig)
- `ILD`: [ILDConfig](https://github.com/iLCSoft/ILDConfig)

## Visualisation

- `k4CEDViewer/CED`:
    - "CLIENT" - k4CEDViewer Codebase: [https://github.com/key4hep/k4CEDViewer](https://github.com/key4hep/k4CEDViewer)
    - "SERVER" - CED Codebase: [https://github.com/iLCSoft/CED](https://github.com/iLCSoft/CED)
- for FCC detectors: `Phoenix@FCC`
    - Codebase: [https://github.com/HEP-FCC/phoenix-at-fcc](https://github.com/HEP-FCC/phoenix-at-fcc)
    - Deployment: [https://hep-fcc.github.io/phoenix-at-fcc/#/](https://hep-fcc.github.io/phoenix-at-fcc/#/) 

Detailed instructions for running visualtion will follow.

## License and Copyright

Copyright (C), k4geo Authors

k4geo is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
