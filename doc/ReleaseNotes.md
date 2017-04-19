# v00-11

* 2017-04-19 Frank Gaede ([PR#94](https://github.com/iLCSoft/lcgeo/pull/94))
  - added versions of the ILD Hcal drivers (with Tesla geometry ) that allow to use multiple readout collections:
        - SHcalSc04_Barrel_v04
        - SHcalSc04_Endcaps_v01
        - SHcalSc04_EndcapRing_v01
  - use these drivers in ILD_l4_v01

* 2017-04-13 Frank Gaede ([PR#93](https://github.com/iLCSoft/lcgeo/pull/93))
  - update documentation for ILD detector models ( see ./ILD/compact/README.md )

* 2017-04-13 Shaojun Lu ([PR#92](https://github.com/iLCSoft/lcgeo/pull/92))
  - Remove this demonstration model, we have one complete ILD_l4_v01 now.

* 2017-04-13 Shaojun Lu ([PR#91](https://github.com/iLCSoft/lcgeo/pull/91))
  - Added ILD_l4_v01 as a multi-technology simulation model.
     - It use the common part of the ILD_l1_v01, and replaced the HCAL by a generic HCAL.
     - used multi-segmentation for RPCHits and SciHits.
     - and saved into two collections: HCalBarrelRPCHits/HCalBarrelSciHits,  HCalEndcapRPCHits/HCalEndcapSciHits, HCalECRingRPCHits/HCalECRingSciHits.
     - reconstruction should use them for their HCAL in this model.

* 2017-04-13 Shaojun Lu ([PR#88](https://github.com/iLCSoft/lcgeo/pull/88))
  - Replaced the ILD_o1_v05(obsolete model) with ILD_l1_v01(optimisation model) in lcgeoTests.
     - ILD will focus on  the optimisation models ILD_l1/s1, and maintain/support for optimisation studies.

* 2017-04-13 Frank Gaede ([PR#93](https://github.com/iLCSoft/lcgeo/pull/93))
  - update documentation for ILD detector models ( see ./ILD/compact/README.md )

* 2017-04-13 Shaojun Lu ([PR#92](https://github.com/iLCSoft/lcgeo/pull/92))
  - Remove this demonstration model, we have one complete ILD_l4_v01 now.

* 2017-04-13 Shaojun Lu ([PR#91](https://github.com/iLCSoft/lcgeo/pull/91))
  - Added ILD_l4_v01 as a multi-technology simulation model.
     - It use the common part of the ILD_l1_v01, and replaced the HCAL by a generic HCAL.
     - used multi-segmentation for RPCHits and SciHits.
     - and saved into two collections: HCalBarrelRPCHits/HCalBarrelSciHits,  HCalEndcapRPCHits/HCalEndcapSciHits, HCalECRingRPCHits/HCalECRingSciHits.
     - reconstruction should use them for their HCAL in this model.

* 2017-04-13 Shaojun Lu ([PR#88](https://github.com/iLCSoft/lcgeo/pull/88))
  - Replaced the ILD_o1_v05(obsolete model) with ILD_l1_v01(optimisation model) in lcgeoTests.
     - ILD will focus on  the optimisation models ILD_l1/s1, and maintain/support for optimisation studies.

* 2017-04-07 Dan Protopopescu ([PR#85](https://github.com/ilcsoft/lcgeo/pull/85))
  - Added correct tracker region to SiD models as per #81
  - Added type flags in ECalBarrel_o2_v03_00.xml
  - Updated model version for SiD test

* 2017-04-07 Shaojun Lu ([PR#80](https://github.com/ilcsoft/lcgeo/pull/80))
  - Use MutliCollections configuration to split SiEcal into preShower inside ILD_o1_v05.
     - "CaloPreShowerSDAction" is not necessary in "ddsim_steer.py"
     - The same  "ddsim_steer.py" in ILDConfig will work for all ILD modules.

* 2017-04-07 Shaojun Lu ([PR#78](https://github.com/ilcsoft/lcgeo/pull/78))
  - Added tracker region in "ILD_o4_v01.xml" to fix the MC Truth link in CED display.

* 2017-04-04 Frank Gaede ([PR#79](https://github.com/ilcsoft/lcgeo/pull/79))
  - enforce definition of tracker volume in ddsim - models need to define:
          - tracker_region_zmax
          - tracker_region_rmax
  - this is needed to ensure that the MCTruth link for hits works correctly
  - allow to enable additional DebugDumpActions with *enableDetailedHitsAndParticleInfo*

* 2017-04-04 Dan Protopopescu ([PR#70](https://github.com/ilcsoft/lcgeo/pull/70))
  * SiD: Updated to-do list

* 2017-04-05 Marko Petric ([PR#82](https://github.com/ilcsoft/lcgeo/pull/82))
  - CLIC_o3_v09_NoGaps model is same as CLIC_o3_v09 but has not gaps in calorimeters between barrel and endcap, the calorimeters endcaps have an inner radius of 1cm, therefore the beampipe and lumi and beam call are removed. Model only meant for calo efficiency study.

* 2017-04-02 bogdan ([PR#72](https://github.com/ilcsoft/lcgeo/pull/72))
  - corrections for LHCal01.xml and LHCal_o1_v01_geo.cpp fixing overlaps.
   -  Inner cutout in LHCal envelope changed from octagon to tube, as for some reasons 
      PolyhedraRegular is not subtracted - causing lack of inner cutout in envelope and overlaps 
      with beam pipe.
  -  Fix of bug in geo driver which caused  internal overlaps

* 2017-04-03 luisaleperez ([PR#77](https://github.com/ilcsoft/lcgeo/pull/77))
  - Removed line in FieldMapXYZ that set B-field to zero if point was outside field-map range.

* 2017-04-03 StrahinjaLukic ([PR#75](https://github.com/ilcsoft/lcgeo/pull/75))
  - Corrected BeamCal envelope for the new dimensions of the graphite absorber.
  - Consolidated some BeamCal parameters. Removed a redundant parameter for the incoming beam pipe radius.
  - Added clearance between BeamCal and the incoming beam pipe.

* 2017-04-03 Marko Petric ([PR#74](https://github.com/ilcsoft/lcgeo/pull/74))
  - New model of CLIC detector CLIC_o3_v09, the only change to CLIC_o3_v08 is the Polystyrene -> G4_POLYSTYRENE in the HCal barrel and endcap. This changes the density of the scintilator from 1.032 to 1.06 g/cm3.

* 2017-03-27 Daniel Jeans ([PR#61](https://github.com/ilcsoft/lcgeo/pull/61))
  - Fix bugs in barrel modules size introduced when generalising barrel symmetry
  - Exit gracefully if too many sides requested for barrel

* 2017-03-29 Frank Gaede ([PR#63](https://github.com/ilcsoft/lcgeo/pull/63))
  - add info printout if sensitive action is changed from the default

* 2017-04-08 Andre Sailer ([PR#86](https://github.com/ilcsoft/lcgeo/pull/86))
  - DDSim: Added printout of execution time at the end of simulation

* 2017-03-24 Emilia Leogrande ([PR#60](https://github.com/ilcsoft/lcgeo/pull/60))
  - Replace ILDCellID0 with LCTrackerCellID

* 2017-03-30 Shaojun Lu ([PR#69](https://github.com/ilcsoft/lcgeo/pull/69))
  - Added one option to fill act list in case that we want to use MultiSegmentation and MultiCollection.

* 2017-03-30 Shaojun Lu ([PR#68](https://github.com/ilcsoft/lcgeo/pull/68))
  - Initialise a module "ILD_o4_v01" for a combined RPC/Scintillator HcalBarrel simulation. 
     - Both RPC and Scintillator have 48 layers.
     - MarlinProcessor may access RPC from slice:3 and access Scintillator from slice:6
     - Assign different segmentations and output collections to the different sensitive type.

* 2017-03-30 Frank Gaede ([PR#67](https://github.com/ilcsoft/lcgeo/pull/67))
  - fix a bug in the field maps  FieldMapBrBz.cpp  and FieldMapXYZ.cpp  that prevented overlay
  - (fixes https://github.com/iLCSoft/lcgeo/issues/65 )

* 2017-03-31 luisaleperez ([PR#73](https://github.com/ilcsoft/lcgeo/pull/73))
  - implemented FieldMapXYZ.cpp 
       - a 3D field map defined on a grid in (x,y,z) 
  - added field maps for the ILD solenoid and anti-DID
       -  fieldmaps/ild_fieldMap_antiDID_10cm_v1_20170223.root 
       -  fieldmaps/ild_fieldMap_Solenoid3.5T_StandardYoke_10cm_v1_20170223.root 
  - added example configuration for using these maps to ILD_o1_v05


* 2017-04-12 Daniel Jeans ([PR#87](https://github.com/ilcsoft/lcgeo/pull/87))
  various updates to beampipe description:
  - inner radius of central tube reduced from 14 to 13.5mm to match DBD
  - tube thickness in conical regions corrected for the cone angle
  - second cone: now made of Be-Cu mix, and given uniform thickness of 2.7mm. this is to simulate a 2mm-thick Be beampipe with 0.7mm of Cu cables around it.


* below are release notes before PR script was introduced

Andre Sailer 2017-03-17 
  - delete NULL is a no-op, fix misleading indentation and useless setting of member variable in d'tor, initialise member pointer to NULL

Andre Sailer 2017-02-28 
  - Add Werror to CI setup

TiborILD 2017-03-22 
  - bug fixes

Dan Protopopescu 2017-03-21 
  - Unchanged from o2_v02

Marko Petric 2017-03-21 
  - Mention LICENCE in README
  - Add LICENCE
  - fix README
  - Update README.md
  - Create README.md

StrahinjaLukic 2017-03-17 
  - Adding test of the BeamCal z location for ILD detectors.
  - Adjusting indentation to surrounding code
  - Removing unneeded sliceType graphiteShielding, correcting indentation.

StrahinjaLukic 2017-03-16 
  - Minor edit in BeamCal08.xml for consistence
  - Resolving incorrect whitespace by hand
  - Correcting code.
  - Restoring another lost line.
  - Lost closing bracketwq
  - Correcting a minor error in the code
  - Removing hardcoded check of BeamCal z.
  - Another minor error...
  - Correcting a minor error in the code
  - Removing hardcoded check of BeamCal z.
  - Changes requested by Andre.

StrahinjaLukic 2017-03-15 
  - Implementing the optional outer_radius attribute for the BeamCal graphite shield as suggested by Andre
  - Changes requested by Andre.
  - Minor edit to the BeamCal driver

StrahinjaLukic 2017-03-14 
  - Adding BeamCal_LHCal.xml for testing relative position of BeamCal and LHCal. Moving BeamCal towards LHCal.
  - Updated segmentation for BeamCal. BC sensor cutout 45 deg. BC phi spanning 360deg-cutout.
  - Adding BeamCal_LHCal.xml for testing relative position of BeamCal and LHCal. Moving BeamCal towards LHCal.
  - Minor fixes and debug messages.
  - Defining top_BCal_dGraphite in top_defs_common and adjusting BeamCal envelope.
  - Updating the BeamCal inner radius and adapting the segmentation.
  - Implementing the optional outer_radius attribute for the BeamCal graphite shield as suggested by Andre
  - BCal_rGraphite has 10mm clearance from LHCal bore. Defined in fcal_defs.xml, removed elsewhere.
  - Proportional segmentation for BeamCal
  - Adding BeamCal_LHCal.xml for testing relative position of BeamCal and LHCal. Moving BeamCal towards LHCal.

StrahinjaLukic 2017-03-13 
  - Adding new parameter BCal_rGraphite to allow graphite shielding to have different outer radius than BeamCal

bogdan 2017-03-18 
  - LHCal envelope correction

bogdan 2017-03-17 
  - asymetric inner cutout, placement fix -  no x-angle rotation
  - info print modified
  - top_Lcal_z_begin set to Ecal_endcap_zmin
  -  LCal aligned with begin ECalPlug and LHCal tied to ECalPug end
  - increased number of layer to 30

Dan Protopopescu 2017-03-15 
  - Update SiTrackerBarrel_o2_v02_00.xml
  - Update SiTrackerEndcap_o2_v02_01.xml
  - Update SiD_o2_v02.xml
  - Added updated Forward Tracker

Daniel Jeans 2017-03-15 
  - add a couple of comments
  - update ECAL: allow other polyhedral barrel shapes (endcaps still octagons); configurable dead region due to slab plug (default 0)

bogdan 2017-03-14 
  -  replace LHcal.xml with new LHCal01.xml in top XML
  -  -update LHcal geometry according to info from yuno@univ.kiev.ua  -new LHCal driver centering detector at outgoing beam pipe   and building LHcal with octagonal inner bore  - new extra utility, FCAL_l1_v01.xml to build FCAL components only

TiborILD 2017-03-08 
  - detecor type flags added/corrected

TiborILD 2017-03-02 
  - correct the cell size for Hcal_EndcapRing_SD in ILD_s2_v01.xml

Frank Gaede 2017-03-03 
  - switch back to constant B-field for ILD_o1_v05

Marko Petric 2017-03-02 
  - Add sourcing of env in SensThickness tests

Marko Petric 2017-03-01 
  - Replace the tube shaped cones with actual tubes

TiborILD 2017-02-28 
  - Hcal_Endcaps_SD_v02 box geometry (as of SHalSc04)
  - Hcal_Endcaps_SD_v02 box geometry (as of SHalSc04)

TiborILD 2017-02-24 
  - Hcal_Endcaps_SD_v01.xml : cleaning & typo correction
  - some typo corrections & cleaning
  - ILD_s2_v01.xml : typo correction for HcalEndcapsCollection
  - ILD_l2_v01.xml : typo correction for HcalEndcapsCollection
  - Hcal_Endcaps_SD_v01.xml: rot_z correction & cleaning

Andre Sailer 2017-02-24 
  - Add guineapig.particlesPerEvent option, including documentation
  - Revert " add readerParameters to input readers"

Frank Gaede 2017-03-01 
  - mv Geant4EventReaderGuineaPig.cpp to DD4hep

Marko Petric 2017-02-28 
  - New driver for simple tube support
  - Keep compilation going after error

Frank Gaede 2017-02-24 
  - implment skipping to event in  GuineaPig reader
  - add parameter ParticlesPerEvent to GuineaPig reader
  -  add readerParameters to input readers
  - add component DDParsers (needed on macos)

Andre Sailer 2017-02-23 
  - TubeX01: drop register, clang warning, useless anyway
  - SServices00: comment unsed private variable (clang warning)
  - Add new/correct spelling for ROOT_INCLUDE_DIRS, root headers properly ignored now

Andre Sailer 2017-02-22 
  - LcgeoExceptions: Fix warning for llvm
  - TPCSDAction: fix shadow and unused variable warnings
  - VXD04: fix shadow and unused variable warnings
  - ECal05, plug: comment unused (Placed)Volumes, variables
  - SECal04_*: comment unsed placedvolume variable
  - SECal05_helper: delete NULL is a no-op we don't have to test before deleting
  - LumiCal_o1_v02: fix warnings for unused variables
  - ODH: static --> inline to avoid unused-function warnings
  - Fix warnings: set but not used and shadowing in detectors used by CLIC

TiborILD 2017-02-21 
  - added  Hcal_Endcaps_SD_v01.xml

Tibor Kurca 2017-02-21 
  - added  Hcal_Endcaps_SD_v01.cpp
  - added  Hcal_EndcapRing_SD_v01.cpp
  - Hcal_Barrel_SD_v01.cpp  5 modules & updates
  - ILD_s2_v01.xml modified for SDHcal
  - ILD_l2_v01.xml modified for SDHcal
  - Hcal_EndcapRing_SD_v01.xml added
  - graphite material added
  - added absorber slice, swapped epoxy/PCB slices
  - HcalSD_ extra parameters names for SDHcal

Dan Protopopescu 2017-02-21 
  - Rename detector/SiTrackerEndcap_o2_v02ext_geo.cpp to detector/tracker/SiTrackerEndcap_o2_v02ext_geo.cpp
  - Delete SiTrackerEndcap_o2_v02ext_geo.cpp

bogdanmishchenko 2017-02-17 
  - Add files via upload

bogdanmishchenko 2017-02-16 
  - Add files via upload

Dan Protopopescu 2017-02-14 
  - Updated driver after closing previous pull request

Marko Petric 2017-02-15 
  - Add the copper cables to the shell of the ourter tracker
  - Move airshell backwards
  - Remove cable from vertex

Marko Petric 2017-02-13 
  - Remove obsolete gitlab CI badge
  - Change surfaces to form a right handed system

Frank Gaede 2017-02-14 
  - split pair files into events in Geant4EventReaderGuineaPig

Frank Gaede 2017-02-13 
  - add vertex info in  Geant4EventReaderGuineaPig

Dan Protopopescu 2017-02-08 
  - Update ECalBarrel_o1_v03_geo.cpp
  - Update SiD_o2_v02.xml

Dan Protopopescu 2017-02-07 
  - Update SiD_o2_v02.xml
  - Update Solenoid_o2_v02_00.xml
  - Create ECalBarrel_o2_v03_00.xml

simoniel 2017-02-08 
  - Added region definitions for G4 material scan.

Frank Gaede 2017-02-08 
  - add ILD_s2_v01 with SDHcal from ILD_o2_v01
  - add ILD_l2_v01 with SDHcal from ILD_o2_v01

Shaojun Lu 2017-02-08 
  -  Following the update of small ILD, the HCAL radius has been reduced by ~34cm (HBU:36cm). Added one version 'SHcalSc04_Endcaps_sv01.xml' for the small AHCAL Endcaps, which has reduced one HUB at x and y, and 14 towers at x direction.(large ILD has 16 towers in endcaps from HBU.)

Andre Sailer 2017-01-27 
  - DDSim: Change default tracker action to Geant4TrackerWeightedAction

Dan Protopopescu 2017-02-01 
  - Delete SiTrackerEndcap_o2_v02_00.xml
  - Rename SiTrackerEndcap_o2_v01_01.xml to SiTrackerEndcap_o2_v02_01.xml
  - Add files via upload
  - Update SiVertexBarrel.xml
  - Rename Solenoid_o2_v01_00.xml to Solenoid_o2_v02_00.xml
  - Rename MuonEndcap_o2_v01_01.xml to MuonEndcap_o2_v02_01.xml
  - Rename MuonBarrel_o2_v01_01.xml to MuonBarrel_o2_v02_01.xml
  - Rename LumiCal_o2_v01_00.xml to LumiCal_o2_v02_00.xml
  - Rename HCalEndcap_o2_v01_01.xml to HCalEndcap_o2_v02_01.xml
  - Rename HCalBarrel_o2_v01_00.xml to HCalBarrel_o2_v02_00.xml
  - More files copied over from SiD_o2_v01_prelim
  - Delete test
  - Rename ECalEndcap_o2_v01_00.xml to ECalEndcap_o2_v02_00.xml
  - Rename ECalBarrel_o2_v01_00.xml to ECalBarrel_o2_v02_00.xml
  - Add files via upload
  - Rename SiVertexEndcap_o2_v01.xml to SiVertexEndcap_o2_v02.xml
  - Rename SiTrackerEndcap_o2_v01_00.xml to SiTrackerEndcap_o2_v02_00.xml
  - Rename SiTrackerBarrel_o2_v01_00.xml to SiTrackerBarrel_o2_v02_00.xml
  - Add files via upload
  - Add files via upload
  - Add files via upload
  - Create BeamPipe.xml
  - Create BeamCal_o2_v02_00.xml
  - Create test
  - Update SiD_o2_v02.xml
  - Update SiD_o2_v02.xml
  - Create SiD_o2_v02.xml

Frank Gaede 2017-02-02 
  - fix orientation of surface vectors in SIT, FTD

Marko Petric 2017-02-01 
  - Remove deprecated drivers to new ones
  - Add test for CLIC_o3_v08
  - Unify all tracker readout IDs with global parameter
  - Add missing materials
  - Add DD4hepVolumeManager
  - Fix the possition of the ECal surface
  - Remove regions since we do not need to go away from the global range cut
  - Fix the problem with missing outer radius for solenoid
  - The whole new tracker
  - Remove warnings from ConicalSupport
  - Fix conical support by removing a obsolete subtraction
  - This is the update of the Inner tracker
  - Update Vertex with cables and redefine materials
  - Add new drivers that read includes inside module block

Dan Protopopescu 2017-01-31 
  - Create README.md
  - Create README.md

Daniel Jeans 2017-01-27 
  - add rec geometry data to SEcal05_Endcaps.cpp driver

Daniel Jeans 2017-01-24 
  - remove SEcal04_ECRing.xml for previous driver
  - SEcal05_siw_ECRing for s1 model; fix FTD envelope definition (go to end of TPC) to avoid overlap
  - introduce SEcal05_ECRing driver: hole centered on outgoing beam; preshower status taken into account. Use this driver in ILD_(l,s)1_v01 models

Frank Gaede 2017-01-24 
  - mv sensitive  plugins; rm obsolete SDs

luisaleperez 2017-01-20 
  - Adding new guinea-pig event plugin
  - Adding new guinea-pig event plugin

luisaleperez 2017-01-19 
  - Using new field-map XYZ
  - Adding new field-map XYZ

Luis Alejandro Perez Perez 2017-01-19 
  - Adding new field-map XYZ

Frank Gaede 2017-01-18 
  - add example conversion:: guineapig_to_lcio.py

# v00-10

Shaojun Lu 2017-01-16
  - Added tube_defs.xml for ILD_s1_v01.

Frank Gaede 2017-01-16 
  - activate TPC_endplate_mix that was used before

Daniel Jeans 2017-01-16 
  - tidy up for new beampipe changes
  - fix lumical overlaps (internal + with beamtube)
  - fix overlaps between rewritten beamipe and ftd
  - fix vtx/beampipe overlaps
  - update beam pipe description

Frank Gaede 2017-01-12 
  - fix TPC cathode volumes (include in cathodeLog)

Shaojun Lu 2017-01-12 
  -  Included 'detector_types.xml' and added missed parameters 'tracker_region_rmax' and 'tracker_region_zmax' for reconstruction.

Shaojun Lu 2017-01-11 
  -  define Yoke relative to Hcal and Coil.

Frank Gaede 2017-01-05 
  - update version and release notes for v00-10

Daniel Jeans 2016-12-26 
  - define yoke barrel size with respect to outer radius of coil
  - tidy up ILD_?1_v01 compact description. Move all common definitions to ILD_common_v01. Keep only top-level steering file and mode-specific definitions in ILD_l1 and ILD_s1 directories. move unused files to extra_stuff directory.

simoniel 2016-12-15 
  - map cellid of surf <--> vector of cellid of neighbouring surf on the same layer (new NeighbourSurfacesStruct) filled in drivers currently used by CLIC

simoniel 2016-12-09 
  - fill of neighboring surfaces map done also for CLIC_o2_v04_p1 drivers and Vertex drivers in use
  - Depence from LCIO re-introduced. Compute neighbours of every surface and fill map of cellID of surface <--> vector of cellID of neighbouring surfaces for zPlanar and zPetalDisk data struct.

Shaojun Lu 2016-12-16 
  -  Updated Yoke05 parameters in compact XML files.
  -  Replaced old parameter 'Hcal_endcap_zmin' with the updated new naming convention parameter 'HcalEndcap_min_z'.

Andre Sailer 2016-12-16 
  - Declare Beampipe, Mask, and Solenoid driver deprecated, copied to DD4hep

Andre Sailer 2016-12-15 
  - Add badge to README, fix typos, modify readme for github
  - Add Travis Configuration

Frank Gaede 2016-12-15 
  - link ${Geant4_LIBRARIES} to TestSensThickness

Shaojun Lu 2016-12-15 
  -  Added a new parameter Coil_Yoke_gap and value 249.0*mm, to scale 'Yoke_barrel_inner_radius' following 'Coil_outer_radius'.

Frank Gaede 2016-12-13 
  - link ROOT libraries to TestSensThickness

Andre Sailer 2016-12-06 
  - Add tests to check that sensitiveThickess is correctly set for the different drivers used in latest clic models
  - set thicknessSensitive for VertexEndcap_o1_v[456] datastructs
  - set thicknessSensitive for TrackerEndcap_o2 datastructs

Marko Petric 2016-12-06 
  - Add comment to TrackerBarrel_o1_v03 and TrackerBarrel_o1_v04 to explain the difference

Daniel Jeans 2016-12-06 
  - enabled service drivers in ILD_l1 and ILD_s1 models

Marko Petric 2016-12-05 
  - Add correct version number

Marko Petric 2016-12-02 
  - Revert the Vertex endcap geometry back to previous one and use new one in the new model
  - Add new model for the updated material budget of the tracker
  - Inherit compiler flags form DD4hep so that TLS is inherited and drop ROOT5

Marko Petric 2016-11-30 
  - Update design to use the new driver
  - Updated vertex endcap design with corrected airflow

Frank Gaede 2016-11-30 
  -  - update release notes and version for v00-09-01

Daniel Jeans 2016-11-30 
  - bug fixes (character to int conversion; strip layer configuration)

Andre Sailer 2016-11-28 
  - Add sensitive thickness to DataStructs for TrackerEndcap drivers

Bogdan Pawlik 2016-11-24 
  -  updated LumiCal for ILD_s1_v01
  -  store layer phi stagger in LayeredCalorimeterData

Shaojun Lu 2016-11-24 
  -  - follow the DD4hep update of DDRec::LayeredCalorimeterData::Layer  - remove usage of DDRec::LayeredCalorimeterData::Layer.thickness    - replace where needed with inner/outer_thickness

Dan Protopopescu 2016-11-23 
  - Files from Amanda
  - Updates by Ross
  - GitHub copy updated by Ross Gordon McCoy

Shaojun Lu 2016-11-23 
  -  Fix the missing 'EcalEndcapParameters', which is needed in Marlin reconstruction. And the 'layoutType' of 'DDRec::LayeredCalorimeterData' will be used by 'convertToGear' to create the correct layout type 'Endcap' for EcalEndcap.

Marko Petric 2016-11-23 
  - Unify ReadoutID with the same schema as in o2_v04_p1

Marko Petric 2016-11-22 
  - Make if forward compatible with CLIC_o3_v07
  - The new pached CLIC_o2_v04 with patched readouts
  - Add tracker version for patch
  - Remove obsolete tracker models
  - Make a patch to CLIC_o2_v04

Daniel Jeans 2016-11-22 
  - fix offset of magic megatiles to give consisent cell indices
  - added forgotten xml files to svn...

Andre Sailer 2016-11-21 
  - SHcalSc04_Barrel_v01: add DetElements for staves, needed for reconstruction

Daniel Jeans 2016-11-21 
  - bug fix (identified by overlap check)

Shaojun Lu 2016-11-21 
  -  Fix the 'SEcal05_Helpers' to follow the update in  the 'DD4hep::DDRec::LayeredCalorimeterStruct::Layer' by removing thickness.
  - remove ECAL preshower layer (barrel&endcap) from ILD_l1_v01 and ILD_s1_v01 models
  - use SEcal05 model for ILD_o3_v05 (scECAL)
  - update ILD_l1_v01 and ILD_s1_v01 to use new SEcal05 driver
  - remove SEcal04Hybrid* (second try...)

Daniel Jeans 2016-11-18 
  - removed SEcal04Hybrid (not compatible with updated Megatile segmentation class); added new SEcal05* drivers

Dan Protopopescu 2016-11-17 
  - Included end plate, air gaps, and ajusted total layer thickness to 19mm

Shaojun Lu 2016-11-17 
  -  Updated 'CaloPrototype_v01.cpp' to be a more generic driver.

Daniel Jeans 2016-11-17 
  - remove new SEcal05 drivers until new segmentation class is available
  - new SEcal05 drivers

Bogdan Pawlik 2016-11-15 
  -  printout modification
  -  Set non-zero phi-offset to fix improper number of phi sectors reported by ddsim, 49 instead of 48

Shaojun Lu 2016-11-15 
  -  Clean up the unused variable 'motherVol'. When DD4hep/lcgeo moved to 'mandatory envelope volume', each sub-detector has been placed into 'envelope' volume. The envelope has been implemented into 'DD4hep/DDCore/src/XML/Utilities.cpp', and placed into the 'mother' volume (world volume) there.
  -  We do not need this user defined envelope shape. The envelope has been implemented into 'DD4hep/DDCore/src/XML/Utilities.cpp'
  -  fix warning unused variable, and commented out the lines. If the users have further implementation with the unused variable, the lines could be used again.
  -  fix warning: declaration of 'RailSeparation' shadows a member of 'this'.
  -  fix many warning for shadowed declaration, and comment out or remove unused variable.
  -  fix warning for shadowed declaration, and remove unused variable.

Aidan Robson 2016-11-14 
  - systemID must be 29 (defined by UTIL::ILDDetID) in CaloFaceEndcapSurfacePlugIn for MarlinTrk extrapolation to calo surface to work.  Not the same as detector id.

Andre Sailer 2016-11-14 
  - CLIC_o3_v07,  LumiCal: Correct segmentation to start first bin at phi=0, not have phi=0 be the center of the 'first' bin which causes problems at phi=+/-pi



# v00-09

  - developers release wih many changes and updates to CLIC, ILD and SiD simulation models, e.g.:
    - made compatible w/ DD4hep v00-18 ( removed LayeredCalorimterData::Layer::thickness )
    - ILD:
      - add models ILD_l1_v01 and ILD_s1_v01
      - new LCal driver
      - new Hcal barrel w/ staircase layout
      - ...
    - CLIC:
      - new mdoels CLIC_o3_v06 and CLIC_o3_v07
      - ...
    - SiD:
      - new model: SiD_o2_v01
      - make compatible w/ reconstruction
      -...


# v00-08

S. Lu
   - Adapted to Hcal endcap Ring to DD4hep envelope, and improve the code to read the envelope information clearly from compact file directly.
   - Updated to read the 'HcalEndcapRing_inner_radius' directly from compact file, and derivative dependence on 'Ecal_endcap_outer_radius' implemented in compact file.
   - added 'Hcal_endcap_thickness' the value came from 26.5*mm*48+15.0*mm. Fixed the 'Hcal_barrel_thickness' the value came from engineer, it is about 26.5*mm*48 and 0.3*mm of tolerance for the 'Hcal_outer_radius' maximum 3395.5*mm
   - Added two parameters 'Ecal_Barrel_thickness', and 'Hcal_Barrel_thickness', and use them to increase the ILD calorimeters (ECAL, HCAL) radius automaticly.
   - Update material for the Ecal ECRing module as Mokka used
   - Change the thickness of the scintillator in the active layer.
   - Updated 'detector/CaloTB/CaloPrototype_v02.cpp' and remove the hardcode layer/ identifier 'layer/K' in the testbeam geometry driver, get this string via segmentation, and use the defination in the compact xml file, whatever the user wanted, 'layer' or, 'K' ... .
   - added ID 'slice' in the geometry driver, and to be used in the cellEncoding string in compact file.
   - Updated SEcal04_Endcaps.cpp to set the Magic Wafer size in group tower, and set the group indentifier in the compact file for both Barrel (in layer group) and Endcaps (in tower group).
   - Use WaferGridXY segmentation for the EcalEndcapsCollection digitization in the DD4hep ILD_o1_v05, to improve the Magic wafer part in each layer too, as EcalBarrelCollection

D. Protopopescu
   - AHCal Barrel implementation of the layout from 70th SiD Optimization Meeting presentation, from Ross Gordon McCoy (ross.mccoy@mavs.uta.edu),
HCalEndcap with the same layer structure, and SiD_o1_v03 including the two new XMLs. Perhaps the Scintillator HCal should be o2 (option 2)?
   - Added new materials, and terminator absorber layer to the Scintillator HCal		   
   - Fix for SiD 'make test'
   - Updated Muon endcap and barrel XMLs with latest dimensions, and using the Generic drivers (will have to change this!). Added envelope (Aidan).
Modified Fe slice thicknesses (20->19.6cm) to fit within the new dimensions.
   - Added test for SiD_o1_v1

M. Petric
   - Make the sum of readout bits 32 otherwise we have a problem with the encoding of the surfaces
   - Fix problem with replication of color name
   - Add interlink to tracker
   - Change CI to new Geant4 10.2.2
   - converged design for tracker
   - Add diagnostic color to gcc 4.9 and move CI build to Ninja and add output on failure to tests
   - New tracker layot as requested by the tracking group
   - Addopt to 40 layers ECal in the endscaps 
   - Add change ECal to 40 layes and move everything in barrel after ECal 27mm outwards
   - Converted readme to markdown and added badge
   - change dd4hep init
   - Make lcgeo ROOT6 compatible. There is no Reflex in 6.

N. Nikiforou
   - Modified CLIC_o2_v04 BeamCal to avoid dummy layer

F. Gaede
   - updated enevelope parameters of ILD_o1_v05 simulation model
   - updated tex file
   - added picture for barrel enevelopes 

K. Coterra
   - Sc_Si_hybrid Ecal drivers for barrel and endcaps were created. So far, they dedicated for Sc-strip ECAL.
   - parameter files for ILD_o3_v05 model, Sc_Si_hybridECAL with AHCAL, were created. So far, for only Sc strip ECal, not for Hybrid ECAL.

Y. Voutsinas
   - updates in VXD material/surface description
   - external cabling and internal strip lines for the innermost VXD layer added
   - adding surfaces for beryllium annulus blocks
   - adding surfaces for the electronics at the end of the ladders
   - adding surface for the beryllium shell cone 

A.Sailer 05/07/2016

  - DDSim: * If no random seed is defined we get a random random seed.
           * new option random.enableEventSeed to calculate reproducible random seeds using the same method as the EventSeeder in Marlin

A.Sailer 08/03/2016

  - DDSim: Add enableG4Gun and enableG4GPS flags to enable the Geant4 Gun or GeneralParticleSource
    see examples/gun.mac or examples/gps.mac
    use with ddsim --enablgeG4Gun --macroFile gun.mac  --compactFile ...

# v00-07

  S.Lu
   - added example for a test beam calorimeter:
      - CaloTB/compact/MainTestBeamSetup.xml
      - ./detector/CaloTB/CaloPrototype_v01.cpp
   - run simulation with: 
      ./CaloTB/run_sim/run_sim.sh
 

# v00-06

F.Gaede
 - added detector type flags to xml files for ILD_o1_v05 and Simplified_ILD_o1_v05
 - updated (ILD) geometry constructors to call XML::setDetectorTypeFlag(e,det)
 - fixed thicknesses for surfaces in FTD (FTD_Simple_Staggered_geo.cpp)
 - added surface for Be endplate of VXD shell (VXD04_geo.cpp)
 - added surface for support shell barrel  (VXD04_geo.cpp)
 - fixed layer structure BeamCal08.xml: include graphite in first layer
 - ILD_o1_v05.xml
   - fixed length of B-field
   - fixed z-position of EcalEndcap caloface surface
   - changed calo face layers to be constructed by plugin now

 - added example reconstruction steering file for SiD: ./example/run_sid_reco.xml
 - fixed material for encap helper surfaces 
   ( moved origin away from centerwhich is in the beam pipe)
    -./other/TrackerEndcapSupport_o1_v01.cpp
    -./tracker/TPC10_geo.cpp
    -./tracker/VXD04_geo.cpp
 - added lcgeo.h defining version macros and
   method versionString()
 - renamed DDSimExceptions.h to LcgeoExceptions.h
 - added RunHeader to lcio file in lcio_particle_gun.py
 - added first example implementation of SiDloi3:
    - ./SiD/compact/sidloi3/sidloi3_v00.xml
 - added connical surfaces to beam pipe:
   ./detector/other/Beampipe_o1_v01_geo.cpp
 - added sensitive action CaloPreShowerSDAction
   - creates pre-shower collections for hits from first layer
   - to be used for (ILD) Ecals
 - protect usage of new member variables in LayeredCalorimeterData
   by usage of #if DD4HEP_VERSION_GE( 0, 15 )
   ( added workaround for this to work with ilcsoft v01-17-08/DD4hep v00-14 )
    - ./detector/calorimeter/GenericCalBarrel_o1_v01_geo.cpp
    - ./detector/calorimeter/GenericCalEndcap_o1_v01_geo.cpp
    - ./detector/fcal/BeamCal_o1_v01_geo.cpp
    - ./detector/fcal/LumiCal_o1_v01_geo.cpp

S.Lu
 - ILD_o1_v05.xml
   - Remove unused 'slice' from ILD calorimeters

  - Added two simplified SEcal drivers by removing detail Si wafer structure in the Si layer, which could 
    improve the performance by saving 9.4 s/event in simulation. Todo: thinking about a sensitive driver (segmentation) 
    for virtual cell and virtual gap from guard_ring and so on:
     - ./detector/calorimeter/SEcal04_Barrel_v01.cpp
     - ./detector/calorimeter/SEcal04_Endcaps_v01.cpp
  - ./detector/calorimeter/SEcal04_Barrel.cpp:
    - Added 'Ecal_fiber_thickness' and 'Ecal_Slab_shielding' thickness to the layer placement position. 
    - Update the Ecal wafer placement position 'wafer_pos_z' with correct dimension parameter to work together with 
      'Ecal_guard_ring_size'

  - Remove 'slice' from the following calorimeter drivers
    - ./detector/calorimeter/Yoke05_Barrel.cpp
    - ./detector/calorimeter/Yoke05_Endcaps.cpp
    - ./detector/calorimeter/SHcalSc04_EndcapRing.cpp
    - ./detector/calorimeter/SHcalSc04_Endcaps.cpp
    - ./detector/calorimeter/SHcalSc04_Barrel_v02.cpp
    - ./detector/calorimeter/SHcalSc04_Barrel_v01.cpp
    - ./detector/calorimeter/SHcalSc04_Barrel.cpp




A.Sailer 02/02/2016
-------------------
  - ddsim
    * Add ExtraParticle tables
    * Add flag to define distance for vertexIsNotEndpointOfParent


N.Nikiforou 12/02/2016
----------------------
  - Changed Solenoid_o1_v01_geo.cpp and SCoil02_geo.cpp for proper filling of 
    solenoid extent in the LayeredCalorimeterStructure for use by DDMarlinPandora

N.Nikiforou 11/02/2016
----------------------
  - Fixed errors due to change of namespace of "_toDouble"
  - Added using of "DetType_AUXILIARY" for CLIC_o2_v04  ECalPlug and HCalRing


N.Nikiforou 10/02/2016
----------------------
  - Modified Tracker*_geo.cpp, Vertex*_geo.cpp, ZPlanarTracker_geo.cpp, GenericCal*_geo.cpp  to add setting detector type via xml tag
  - Moved CLIC_o2_v04 to using TrackerEndcap_o1_v04 drivers that set ZPetalDisksData according to minimum/maximum R and average Z
    of all rings in one endcap tracker disk (for DDMarlinPandora)
  - Added "type_flags" tags to all CLIC_o2_v04 subdetector xml files and added include of detector_types.xml

N.Nikiforou 03/02/2016
----------------------
  - Added deprecation warning for GenericSurfaceInstaller plugin which is now part of DD4hep. 
  - Moved CLIC_o2_v04 to using DDhep_GenericSurfaceInstallerPlugin rather than the local lcgeo implementation which is deprecated

A.Sailer 18/01/2016
-------------------
  - ddsim
    * Add RangeCut parameter
    * Add LargestAcceptableStep for magneticField (10*m)

N.Nikiforou 15/12/2015
----------------------
  - Switched CLIC_o2_v04 to using the GenericSurfaceInstaller. Changed therefore to using the following
    includes (notice _02 at the end of the name):
    * Vertex_o2_v04_02.xml
    * InnerTracker_o2_v03_02.xml
    * OuterTracker_o2_v03_02.xml
  - Changed the orientation of the endcap uvn vectors for Vertex and Tracker to 
    u(-1.,0.,0.), v(0.,0.,1.), n(0.,1.,0.), o(0.,0.,0.) using the GenericSurfaceInstaller rather than
    the individual surface installers (for which the orientation of uvn remains as before). The 
    vectors were set using plugin arguments in the xml. The resulting orientation puts v along R
    i.e. the "bad" measurement direction and u along -phi (clockwise). 
  - GenericSurfaceInstaller: by default all u,v,n,o vector components are now set to 0 which means 
    that the user has to define the non-zero components as arguments, otherwise the plugin fails
    
A.Sailer 15/12/2016
-------------------
  - ddsim
    * Add default filter with Edep > 0

N.Nikiforou 10/12/2015
----------------------
  - Added GenericSurfaceInstaller in detector/other to install surfaces to a given
    box-like volume. The u,v,n,o vectors are provided component-wise as arguments to the plugin.
    If the volumes are trapezoids (as in the case of the endcap) the thickness dimension is supposed to be 
    the same for both sides of the trapezoid (for example, dy1=dy2=dy) which is supposed to be handled
    when "casting" the Trapezoid shape into a Box
  - Added InnerTracker_o2_v03_02.xml and OuterTracker_o2_v03_02.xml in CLIC_o2_v04 (not used by default)
    as examples using the new GenericSurfaceInstaller
  
N.Nikiforou 04/12/2015
----------------------
  - Modified CLIC_o2_v04 to use proper convention for listing slices in module composition.
    In this convention, when you list a layer or a slice/component in a module, the top-most
    element in the xml is the component that is closest to the IP (top-most is inner-most).
  - Created new drivers and new xml files to use them
     * TrackerBarrel_o1_v03 and TrackerEndcap_o1_v03 used by Inner/OuterTracker_o2_v03.xml
     * VertexEndcap_o2_v03 used by Vertex_o2_v04_01
  - Created new surface installer plugins copied from DD4hep that follow the convention discussed above
    and used by the new drivers. They are listed below with the u,v,n and origin vectors as defined in 
    the driver coordinate system :
     * TrackerBarrelSurfacePlugin: u(-1.,0.,0.), v(0.,-1.,0.), n(0.,0.,1.), o(0.,0.,0.)
        + Assumes the driver builds the module along the original +z axis, with the top-most slice
          stacked perpendicular to z, closest to (0,0,0) with eventual rotations to place the module
          in the barrel so that the slice is inner-most. This is the case now in TrackerBarrel_o1_v03 
     * TrackerEndcapSurfacePlugin: u(0.,0.,1.), v(1.,0.,0.), n(0.,1.,0.), o(0.,0.,0.);
        + Assumes the driver builds the module along the original +y axis, with the top-most slice
          stacked perpendicular to y, closest to (0,0,0) with eventual rotations to place the module
          in the barrel so that the slice is inner-most. This is the case now in TrackerEndcap_o1_v03 

A.Sailer 30/11/2015
----------------------
  - Added steering of filter to ddsim
    * assign filters to sensitive detectors
    * create your own filters
  - Added dumpSteeringFile option to ddsim
    * prints steeringFile to stdout with current values as default

A.Sailer 23/11/2015
----------------------
 - Added information to the lcio runHeader in ddsim. This needs the DD4hep revision 1984!
   * The steeringFile content, commandline and resulting ddsim parameters are written
   * geant4 and dd4hep version is written directly in dd4hep.
   * user, workingdirectory and date are written as well

A.Sailer 06/10/2015
----------------------
 - ddsim python simulator
 * "ddsim --help" to print all the command line parameters currently available
 * example/steeringFile.py to see an example steering file with some options
 * to enable tab-completion for ddsim: requires bash (or zsh),bash-completion, and python argcomplete
   (e.g. install with: [sudo] pip install argcomplete)
   in the shell run: eval "$(register-python-argcomplete ddsim)"


N.Nikiforou 24/08/2015
----------------------
 - Modified CLIC_o2_v03/ECalBarrel_o2_v01_01.xml, ECalEndcap_02_v01_01.xml and
   ECalPlug_o1_v01_01.xml to avoid having layers without active element
   which complicates implementation for reconstruction.
   This was the case for the last (outermost) dead layer (absorber only).
   Instead, added one more type of layer with only one repeat which includes the extra
   absorber plate
 - Fixed bug (?) in CLIC_o2_v03/ECalEndcap_02_v01_01.xml and ECalPlug_o1_v01_01.xml
   which had different siPCBMix thickness compared to the ECalBarrel (1.20 mm -> 0.80 mm). 


N. Nikiforou 19/08/2015 
-----------------------

   M.Petric
   - Made a new version for CLIC named CLIC_o2_v03 which includes better representations
     of X0 in the Trackers
   N. Nikiforou
   - Added new Generic Calorimeter drivers based on Polyhedra Calorimeters named
     GenericCalBarrel and GenericCalEndcap. The Generic drivers
     were obtained by svn cp from HCalBarrel and HCalEndcap respectively and are used
     in the CLIC_o2_v03 model by the ECal, HCal and Yoke Endcaps, Barrels and Rings.
     These new drivers (almost identical to the old ones) will be used to avoid making 
     changes in all the drivers while developing the reconstruction framework.
     Changes include:
     - Cleaning up of old comments, old code snippets and featurs no longer used (like
       multiple sensitive elements per layer which will not be supported by the reco).
     - Introduction of the new LayeredCalorimeterStruct variables defined in DDRec:
       - inner/outer thicknesses, nRad/nIntLengths and sensitive element thicknesses
     - Avoidance of cloning of element volumes in the barrel case
     - Start of counters (e.g. layer, module, stave, etc) from 0 rather than 1

# v00-05

   N. Nikiforou
   - Added TrackerEndcap_o1_v02_geo.cpp (svn cp'd from TrackerEndcap_o1_v01_geo.cpp) which supports 
     ILD ID encoding and implements ring modules. Elements (enumerated as "sensors") in one ring 
     (constant R) have the same "module" number but different sensor number. 
     Also modified the CLIC_o2_v02 model compact (InnerTracker_o2_v01_01.xml 
     and OuterTracker_o2_v01_01.xml) to use the new driver. Also now implements the "side" bitfield
     side=0 is positive, side=1 is negative. 
   - Introduced new driver VertexEndcap_o1_v02 and necessary changes to compact file for CLIC_o2_v02.
     The new xml also implements double layers in the endcap but with each layer having its own ID. 
     Also, in the new VertexEndcap driver similar changes as above VertexEndcap.


   S.Lu 
   - updated Ecal and Hcal in Share_ILD_o1_CLIC

   T.Quast
   - fixed DetectorData structure for several drivers
     ( as needed for drawing the detector w/ DD4hep in CED )


   M. Petric
   - added CLIC_o2_v02 model w/ simplified Ecal barrel
     ECalBarrel_o2_v01_01.xml in order to set correct W thicknesses

   F.Gaede

   - added TPCSDAction.cpp 
     - ported from Mokka/TPCSD04.cc
   - introduces Geant4 dependency 
      -> to be addressed ...
   - activated in TPC10_geo.cpp
     -> hits should be exactly on pad row centers ...

   - updated ddsim.py :
      - use Geant4ScintillatorCalorimeterAction
        as defaults for all calorimeters
      - use TPCSDAction for TPC 


   A.Sailer
     - Implemented FieldMapBrBz.cpp, based on the 2D Fieldmap of Mokka FieldX03
     - Example XML for the fields section of the compact XML

            <field name="DetectorMap" type="FieldBrBz"
                   filename="${lcgeo_DIR}/fieldmaps/ILDMap_KB_20150204_BRhoZ.root"
               tree="fieldmap:rho:z:Brho:Bz"
               rScale = "1.0"
               zScale = "1.0"
               bScale = "1.0"
               rhoMin = "5*mm"
               zMin = "5*mm"
               rhoMax = "7005*mm"
               zMax = "7005*mm"
               nRho = "701"
               nZ = "701"
               >
            </field>

     - added ./fieldmaps/ILDMap_KB_20150204_BRhoZ.root 
       latest field simulation for ILD by K.Buesser

     - requires BOOST

   N. Nikiforou
     - added  EcalBarrelFace_v00.xml/EcalEndcapFace_v00.xml to CLIC_o2_v01
     - Enabled detailed shower mode for calorimeters by default in ddsim.py/DD4hepSimulation.py

   F.Gaede
     - added PolyhedralBarrelSurfaces_geo.cpp/PolyhedralBarrelSurfaces_geo.cpp
     - to be used for track states at the calorimeters
     - added as EcalBarrelFace_v00.xml/EcalEndcapFace_v00.xml to ILD_o1_v05


# v00-04-01: patch release


   F.Gaede

    - added a new starting point for SDHcal development: ILD_o2_v01
       - based on ILD_o1_v05 as released in (v00-04/ilcsoft v01-17-07)
       - copied SHcalSc04_Barrel_v01.cpp to Hcal_Barrel_SD_v01.cpp
       - copied SHcalSc04_Barrel_v01.xml to Hcal_Barrel_SD_v01.xml

    - removed outdated ILD_o2_v00


     N. Nikiforou
     - Changed CLIC calorimeters to use new isRadiator() helper function and 
       new extended LayeredCalorimeterStruct

# v00-04:  ****** renamed to lcgeo ******************

     - many changes and developements:

     - first 'complete' prototypes of ILD_o1_v05 
       and CLIC_o2_v01 simulation models

     - introduced envelopes for all sub detector

     - added examples/ddsim.py for running the 
       simulation

     - plugins for adding surfaces to the CLIC model
     - ... 

//======================================================================

# v00-03: third beta release of DDSim

       - ...


# v00-02: second beta release of DDSim

       M.Frank:  
       - add initial version of DDEve.xml
       - create header file with class XMLHandlerDB

       Sh.Lu:
       - implement EndcapRings for Ecal and Hcal
       - many fixes to calorimter drivers


       A.Sailer:
       - implement 'canonical' SD for creating lcio::SimCalorimeterHits
       	 using the readout/segmentation to compute cellID and position  
       - fixes for BeamCal driver


# v00-01: first beta release of DDSim

	     - ILD detector (ported from Mokka ILD_o1_v05) with :
                 VXD, SIT, FTD, SET, TPC,
                 EcalBarrel, HcalBarrel, EcalEndcap, HcalEndcap, EcalRing, beamcal
             - sensitive detectors still experimental
