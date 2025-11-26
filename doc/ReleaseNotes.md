# v00-24

* 2025-11-26 Lorenzo Pezzotti ([PR#545](https://github.com/key4hep/k4geo/pull/545))
  - Set optical fiber radius identical for IDEA_o2 endcap and barrel hadronic calorimeters
  - Do not include IDEA_o2 hadronic barrel calorimeter casing by default

* 2025-11-26 Daniel Jeans ([PR#540](https://github.com/key4hep/k4geo/pull/540))
  - Updated ILD_FCCee_v0{1,2} to use MDI description MDI_o1_CADBased_v01
  - world volume material changed to vacuum (to ensure vacuum within beampipe)

* 2025-11-21 sss ([PR#535](https://github.com/key4hep/k4geo/pull/535))
  - Updated HCal segmentation classes to more efficiently store cell edge information as vectors rather than as hash maps.

* 2025-11-19 Juan Miguel Carceller ([PR#546](https://github.com/key4hep/k4geo/pull/546))
  - Set the env variable K4GEO in thisk4geo.sh to the current source directory, since in the installation there may be or not a folder with the compact files (depending on the cmake option `INSTALL_COMPACT_FILES`)

* 2025-11-17 Lorenzo Pezzotti ([PR#539](https://github.com/key4hep/k4geo/pull/539))
  - Avoid energy deposition filters for the scepcal in the IDEA_o2 steering file for correct usage of its SDAction
  - Reproduce the expected experimental light yields and light fluctuations for the scepcal main crystals and timing layer
  - Do not simulate Cerenkov light in the scepcal timing layer

* 2025-11-17 sss ([PR#536](https://github.com/key4hep/k4geo/pull/536))
  - Fixed cmake development warnings.

* 2025-11-13 Giovanni Marchiori ([PR#543](https://github.com/key4hep/k4geo/pull/543))
  - fix comments in turbine segmentation class

* 2025-11-05 Juan Miguel Carceller ([PR#531](https://github.com/key4hep/k4geo/pull/531))
  - Remove checks for old versions of DD4hep (10+ years old)

* 2025-11-05 SanghyunKo ([PR#529](https://github.com/key4hep/k4geo/pull/529))
  - Allow negative time optical photons (from BIB overlay) in the DRCaloSDAction

* 2025-10-31 aciarma ([PR#534](https://github.com/key4hep/k4geo/pull/534))
  - [FCCee] MDI_o1_ShapeBased_v01 started from MDI_o1_v00 
  - [FCCee] added antisolenoid field map to ShapeBased
  - [FCCee] added final focus quadrupoles and their field to ShapeBased. Dimension and magnetic field defined via parameter xml based on different accelerator lattice versions and beam energy
  - [FCCee] add quadrupoles also to MDI_CADbased_o1_v01"

* 2025-10-30 BrieucF ([PR#532](https://github.com/key4hep/k4geo/pull/532))
  - [FCCee] Remove the detector solenoid field from the CAD imports (it is defined by the detector xml's)
  - [FCCee] Remove the compensating solenoid analytical field map from the detector xml
  - [FCCee] Switch the commented CAD beampipe to a newer version and add a workaround to have vacuum inside the pipe

* 2025-10-29 aciarma ([PR#530](https://github.com/key4hep/k4geo/pull/530))
  Added field map from Mike Koratzinos for the magnetic field of IR antisolenoids.
  The fieldmap accounts only for the antisolenoids, so it must be overlapped with the detector field.
  For the moment I used an example description for the 2T detector solenoid field (fields_solenoid.xml) but this part may be different for each experiment and therefore in each dedicated xml, not in /MDI.
  The antisolenoids are common instead.
  
  Example of the field on the ZX plane from `dumpBfield`:
  <img width="927" height="587" alt="image" src="https://github.com/user-attachments/assets/435946c8-d8ab-416b-8576-9f205401613d" />

* 2025-10-24 sss ([PR#524](https://github.com/key4hep/k4geo/pull/524))
  - Updated FCCSWHCalPhiTheta_k4geo and FCCSWHCalPhiRow_k4geo to avoid the thread-unsafe use of mutable members.

* 2025-10-22 Stefano Fanchellucci ([PR#497](https://github.com/key4hep/k4geo/pull/497))
  - Harmonize object ID for Vertex Barrel

* 2025-10-21 aciarma ([PR#528](https://github.com/key4hep/k4geo/pull/528))
  * FCCee/MDI: Moved the `stl_files` folder up one level:
    from `FCCee/MDI/compact/MDI_o1_v01/stl_files` to `FCCee/MDI/compact/stl_files`
    Provide a common space for CAD versions sharing the same files.

* 2025-10-17 aciarma ([PR#527](https://github.com/key4hep/k4geo/pull/527))
  CAD version of MDI had the crotch rotated 180deg, so that the notch for the SR mask was on the wrong side.
  The SR mask is only on the lower part of the upstream beampipes - negative x coordinates.
  
  Created new folder MDI_o1_CADBased_v01 with the following fixes:
  - fixed orientation of crotch beam pipe
  - added MDI_standalone_v01.xml files to ease debug of MDI geometries
  
  New naming convention decided:
  - MDI_o1_ShapeBased_vXX for shape based geometries
  - MDI_o1_CADBased_vXX    for CAD based geometries

* 2025-10-14 mahmoudali2 ([PR#526](https://github.com/key4hep/k4geo/pull/526))
  This PR introduces improvements to the muon system of the IDEA detector:
  1. Automated Endcap Rmax Calculation:
  - The maximum radius (Rmax) of the endcap is now calculated automatically based on barrel detector geometry, which its last layer encloses the endcap within it. Resulting in reducing the need for parameters manual adjustments.
  - It also simplifies changes to the number of layers, as Rmax no longer needs to be adapted manually for each design.
  
  2. Fixed Negative Endcap Layering Order:
  
  - Corrected the ordering of layers in the negative endcap to ensure proper reconstruction and simulation consistency.
  - This fix addresses issues where layer sequences in negative endcap were previously inverted, improving track propagation and hit association.

* 2025-10-13 Juan Miguel Carceller ([PR#525](https://github.com/key4hep/k4geo/pull/525))
  - Do not run tests that require the beampipe files when not downloading them, since they will fail

* 2025-10-13 Andrea Sciandra ([PR#523](https://github.com/key4hep/k4geo/pull/523))
  - FCC-ee ALLEGRO detector concept: include plugins attaching additional info to the geometry, needed to generate `edm4hep::TrackerHitPlane`.

* 2025-10-06 Alvaro Tolosa Delgado ([PR#522](https://github.com/key4hep/k4geo/pull/522))
  - Refactor of the ARC DDSim configuration example, used in CI. The main change is that hit cellID assignment now works correctly, thanks to consistent handling between the optical physics configuration (SD not triggered at boundaries) and the DD4hep SD action for the optical tracker.

# v00-23

* 2025-10-01 sss ([PR#517](https://github.com/key4hep/k4geo/pull/517))
  - Code cleanups in the segmentation classes.

* 2025-09-30 Paolo Andreetto ([PR#508](https://github.com/key4hep/k4geo/pull/508))
  - Add the MuSIC detector concept geometry. (Ported from muon collider fork)

* 2025-09-25 sss ([PR#516](https://github.com/key4hep/k4geo/pull/516))
  - Small cleanup of HCalTileBarrel_01_v02_geo: take the number of segments from the XML file rather than recalculating in C++.

* 2025-09-19 sss ([PR#514](https://github.com/key4hep/k4geo/pull/514))
  - Fix overlaps in the tile hcal barrel.

* 2025-09-18 Juan Miguel Carceller ([PR#515](https://github.com/key4hep/k4geo/pull/515))
  - Add override to fix warnings in Clang, fix variable names in docs and add missing includes to the FCCSWHCalPhi segmentations

* 2025-09-18 sss ([PR#513](https://github.com/key4hep/k4geo/pull/513))
  - Update outdated comment in FCCSWGridModuleThetaMerged_k4geo::position.

* 2025-09-18 sss ([PR#511](https://github.com/key4hep/k4geo/pull/511))
  - Make sure segmentation classes properly override the virtual neighbours() method for the FCCSWHcalPhi segmentations.

* 2025-09-15 sss ([PR#512](https://github.com/key4hep/k4geo/pull/512))
  - Small code cleanup in a few files to remove the use of exceptions for normal control flow.

* 2025-09-10 Juan Miguel Carceller ([PR#510](https://github.com/key4hep/k4geo/pull/510))
  - Remove a default destructor that is not needed in `DRparamBase_k4geo.h`

* 2025-09-10 Sanghyun Ko ([PR#486](https://github.com/key4hep/k4geo/pull/486))
  - Use the front end instead of the rear end for the default DRC calo hit position
  - Added a dedicated function to calculate the rear end position of the fiber DRC cells (only needed for the DRC digitizer)
  - Made the neighborhood definition of the fiber DRC configurable via a compact file
  - Fixed a bug that the fiber DRC SD action did not take the full active volume of the scintillation fiber for the energy deposit

* 2025-09-09 Robyn Evren ([PR#499](https://github.com/key4hep/k4geo/pull/499))
  - New variables added for partial construction of SCEPCal

* 2025-09-09 Robyn Evren ([PR#488](https://github.com/key4hep/k4geo/pull/488))
  - added variables to allow for choice of phi segment and number of towers loaded.
  - altered for-loops for continuous phi segments and forward only tower placement

* 2025-08-27 Leonhard Reichenbach ([PR#507](https://github.com/key4hep/k4geo/pull/507))
  - fix: [CLD_o2_v08] set tracking_volume parallelworld to not be connected

* 2025-08-18 Thomas Madlener ([PR#505](https://github.com/key4hep/k4geo/pull/505))
  - Set the helicity for the MCParticle once it becomes available. See [EDM4hep#404](https://github.com/key4hep/EDM4hep/pull/404) and [DD4hep#1488](https://github.com/AIDASoft/DD4hep/pull/1488) for more details.

* 2025-08-18 Andre Sailer ([PR#504](https://github.com/key4hep/k4geo/pull/504))
  - Lumical_o1_v01: only set calorimeter type flags for detectors with sensitive elements, otherwise use SUPPORT and AUXILIARY. Fixes #495

* 2025-08-11 Federico Meloni ([PR#503](https://github.com/key4hep/k4geo/pull/503))
  - Muon collider compact files are installed in k4geo/share

* 2025-08-05 jmcarcell ([PR#502](https://github.com/key4hep/k4geo/pull/502))
  - Update tests to allow running them without installing (i.e. we don't need to run `make install` or `ninja install` before running the tests
  - Remove "test_IDEA_with_DRC_o1_v03" that was duplicated, keep only the one that simulates electrons

* 2025-07-31 Kevin Nelson ([PR#480](https://github.com/key4hep/k4geo/pull/480))
  Add a straw tube tracker for ALLEGRO concept.  The implementation was discussed at the FCC Full Sim Meeting on May 28: https://indico.cern.ch/event/1551837/#71-straw-tube-tracker-implemen.  This discussion included the memory load of the straw concept, which is on the order of 200MB.  The logical volume count is minimized by sharing logical tube volumes for each multilayer, however due to the use of cylindrical tube geometry the number of faces is large.
  
  The straw tracker concept uses thin-wall straws (12um) to achieve high spatial resolution.  The main possible advantage over a drift chamber concept is a modular construction with each readout channel corresponding to a single unit.  Additionally, because each tube is a single unit there are no B field correction factors needed to compute spatial position, and it is expected that the straws may achieve a slightly higher single hit position resolution compared to a drift chamber.
  
  Changes:
   - Create new ALLEGRO version o2_v01 with straw tracker
   - Corresponding READMEs in ALLEGRO and detector/tracker
   - Straw tracker C++ implementation in detector/tracker folder
   - Uses same dimensions as drift chamber for easy interchange of the two tracker concepts
   - Both .xml and .cpp code for straw tracker have readable, useful comments
   - For ALLEGRO_o2_v01 ONLY (the straw tube implementation), extend the silicon wrapper z extent to +/-3m to make full use of the available space
   - Default straw tracker implementation has +/-2 degree stereo angle alternating for even/odd multilayers
   - Straw tracker tube size matches IDEA for innermost cells (8mm) and increases linearly with radial distance (saturating at 2.5cm diameter tubes) from interaction point to minimize total N tubes, which minimizes multiple scattering and improves resolution.

* 2025-07-29 jmcarcell ([PR#500](https://github.com/key4hep/k4geo/pull/500))
  - Use . instead of -> for objects that are not pointers

* 2025-07-29 Daniel Jeans ([PR#483](https://github.com/key4hep/k4geo/pull/483))
  - FieldMapBrBz: Correct bug in the conversion of a radial field component to x-y components. Previously x and y were interchanged.

* 2025-07-25 Andreas Loeschcke Centeno ([PR#494](https://github.com/key4hep/k4geo/pull/494))
  - Updated the `stereosign` bitfield in the IDEA DriftChamber to host two bits

* 2025-07-18 varnes ([PR#492](https://github.com/key4hep/k4geo/pull/492))
  This PR incorporates two bug fixes and one new feature in ECalEndcap_Turbine_o1_v03:
  
  Bug fixes: 
  
  1) Fix a 1 cm overlap where outer support ring extended beyond the LAr bath
  2) Fix duplicated volume errors when running in calibration mode (i.e. with electrode, absorber, glue, and cladding all set to sensitive
  
  New feature: add mechanical support structure behind calorimeter.  Note that this is a placeholder at the moment, with the final shape/size pending engineering studies.

* 2025-07-14 Sanghyun Ko ([PR#493](https://github.com/key4hep/k4geo/pull/493))
  - Fix build failure when `K4GEO_USE_LCIO=OFF`

* 2025-07-08 Giovanni Marchiori ([PR#490](https://github.com/key4hep/k4geo/pull/490))
  - Change ALLEGRO calorimeter types from BirksLawCalorimeterSD and SimpleCalorimeterSD to Geant4ScintillatorCalorimeter. This does not affect the behaviour (the same action is used) but  removes the "Sensitive type: Unknown" warning in the output.

* 2025-07-07 Giovanni Marchiori ([PR#489](https://github.com/key4hep/k4geo/pull/489))
  - Fix Rmax of muon tagger barrel and response type of muon tagger endcap for ALLEGRO

* 2025-07-03 Andre Sailer ([PR#478](https://github.com/key4hep/k4geo/pull/478))
  - CLD_o2_v08: new detector model taking into account Birks Law for the scintillator hcal.

* 2025-06-30 BrieucF ([PR#487](https://github.com/key4hep/k4geo/pull/487))
  - Improve ALLEGRO_o1_v03 ECAL visual rendering

* 2025-06-25 wonyongc ([PR#443](https://github.com/key4hep/k4geo/pull/443))
  - Added Segmented Crystal EM Precision Calorimeter (SCEPCal)
  - Separate Main/Timing layer subdetectors
  - Separate readout collections for edep, S/C photon counts using Geant4CalorimeterHit/edm4hep::SimCalorimeterHit
  - Support different F/R crystal granularities in main layer
  - Single, pure-projective timing layer
  - New main layer projective crack migitation scheme with angular offset
  - Added OpticalSurfaces

* 2025-06-23 Nazar Bartosik ([PR#469](https://github.com/key4hep/k4geo/pull/469))
  - Add the "historic" (i.e. not MAIA or MuSIC) muon collider geometries to `MuColl/MuColl/compact`

* 2025-06-19 Federico Meloni ([PR#468](https://github.com/key4hep/k4geo/pull/468))
  - Add the `MAIA` muon collider detector concept

* 2025-06-16 Sanghyun Ko ([PR#485](https://github.com/key4hep/k4geo/pull/485))
  - Add detector type flag for IDEA o1 DRC

* 2025-06-12 jmcarcell ([PR#482](https://github.com/key4hep/k4geo/pull/482))
  - Fix pre-commit after updating to LLVM 20

# v00-22

* 2025-05-27 Archil Durglishvili ([PR#479](https://github.com/key4hep/k4geo/pull/479))
  - ALLEGRO_o1_v03: Birks constant value is set for Polystyrene scintillator used by HCal. This fixed the abnormal response to hadrons.

* 2025-05-27 Zhibo Wu ([PR#458](https://github.com/key4hep/k4geo/pull/458))
  - Crosstalk coefficients can be provided as input argument in the `det::crosstalk::getNeighboursModuleThetaMerged` function
  - No impact on other parts of k4geo, but [k4RecCalorimeter](https://github.com/HEP-FCC/k4RecCalorimeter/blob/main/RecFCCeeCalorimeter/src/components/CreateFCCeeCaloXTalkNeighbours.cpp) needs to be updated accordingly.

* 2025-05-22 Sanghyun Ko ([PR#477](https://github.com/key4hep/k4geo/pull/477))
  - Follow up fix of #471 & protecting from an exception by `std::map::at`

* 2025-05-14 Thomas Madlener ([PR#470](https://github.com/key4hep/k4geo/pull/470))
  - Add a `CODEOWNERS` file to designate ownership to different parts of the repository in order to facilitate review and approval processes.

* 2025-05-13 Sanghyun Ko ([PR#471](https://github.com/key4hep/k4geo/pull/471))
  - Developed neighboring CellID finding algorithm for IDEA_o1 DRC
  - Bug fix of the CellID position calculation of `GridDRcalo_k4geo`

* 2025-04-29 Leonhard Reichenbach ([PR#453](https://github.com/key4hep/k4geo/pull/453))
  - CLD_o2_v07: made tracking volume actually parallel

* 2025-04-25 Thomas Madlener ([PR#467](https://github.com/key4hep/k4geo/pull/467))
  - Add a `clang-format` hook to our pre-commit configuration to actually enforce the formatting that we expect via `.clang-format`.

* 2025-04-24 Erich Varnes ([PR#456](https://github.com/key4hep/k4geo/pull/456))
  This PR adds the ability to merge readout cells with the same z and rho indices in adjacent absorber blades.  The number of blades to merge in each wheel is controlled by the xml file.
  
  Also, some bugs in calculating the positions of the cells were fixed.

* 2025-04-17 Sanghyun Ko ([PR#466](https://github.com/key4hep/k4geo/pull/466))
  - Revert dual-readout version change of `README.md` in #465

* 2025-04-14 Sanghyun Ko ([PR#465](https://github.com/key4hep/k4geo/pull/465))
  - Added DDRec extension (`dd4hep::rec::LayeredCalorimeterData`) for DRC IDEA_o1 to provide necessary detector dimensions to the downstream

* 2025-04-14 sss ([PR#463](https://github.com/key4hep/k4geo/pull/463))
  - Fix a typo in comments.

* 2025-04-13 Sanghyun Ko ([PR#459](https://github.com/key4hep/k4geo/pull/459))
  - Switch monolithic DRC's calo hit type to edm4hep::SimCalorimeterHit to store positions

* 2025-04-12 sss ([PR#462](https://github.com/key4hep/k4geo/pull/462))
  - Adjust ALLEGRO HCal endcap (HCalThreePartsEndcap) cell positions so that they are symmetric under z-inversion.

* 2025-04-11 Giovanni Marchiori ([PR#461](https://github.com/key4hep/k4geo/pull/461))
  - Implement cellDimensions method for two segmentation classes used in the readout of ALLEGRO subdetectors
  - Properly calculate and save as caloData the cell sizes for the noble liquid inclined barrel EM calorimeter, for the HCAL barrel calorimeter and for the SimpleCylinder detector used for the muon system

* 2025-04-11 mahmoudali2 ([PR#444](https://github.com/key4hep/k4geo/pull/444))
  - Last touches on muon-system adoption for DDPlanarDigi( it finally works with MS).
  - Reduced number of segmentation.
  - Created a string for muon system readout to be called in DDPlanarDigi. (`${MuonSystemReadoutID}`)
  - Made the surfaces dimensions = 2 instead of 1.

* 2025-04-08 Victor Schwan ([PR#428](https://github.com/key4hep/k4geo/pull/428))
  - change the MDI version of the ILD FCCee models to `MDI_o1_v01`
  - minor: remove some trailing whitespaces in the ILD FCCee config files

* 2025-04-04 Sanghyun Ko ([PR#454](https://github.com/key4hep/k4geo/pull/454))
  - Migrate SiPM efficiency to HEP-FCC/k4RecCalorimeter#145

* 2025-03-26 jmcarcell ([PR#452](https://github.com/key4hep/k4geo/pull/452))
  - Remove the NeighbourSurfacesData extension in `TrackerBarrel_o1_v06_geo` since it is not used, so we don't need to fill it every time a detector using `TrackerBarrel_o1_v06_geo` is built

* 2025-03-25 jmcarcell ([PR#451](https://github.com/key4hep/k4geo/pull/451))
  - Fix compiler warnings caught with clang 19
  - Add a missing `#include <list>` that will be needed with https://github.com/AIDASoft/DD4hep/pull/1432

* 2025-03-25 jmcarcell ([PR#450](https://github.com/key4hep/k4geo/pull/450))
  - Exit a loop earlier and cache a map lookup in TrackerBarrel_o1_v06_geo

* 2025-03-25 armin.ilg ([PR#442](https://github.com/key4hep/k4geo/pull/442))
  - Fixing IDEA vertex detector dimensions
  - Covering gaps in detector coverage at theta=90 degrees

* 2025-03-20 Giovanni Marchiori ([PR#448](https://github.com/key4hep/k4geo/pull/448))
  - Fix envelope of endcap turbine saved in detector data extension for reconstruction

* 2025-03-19 Giovanni Marchiori ([PR#446](https://github.com/key4hep/k4geo/pull/446))
  - Update ALLEGRO tracker_region zmaz, rmax constants to also include the wrapper

* 2025-03-19 Giovanni Marchiori ([PR#445](https://github.com/key4hep/k4geo/pull/445))
  - Set detector type flag for the IDEA_o1_v03 drift chamber (this propagates to other detectors using the drift chamber)

* 2025-03-19 Giovanni Marchiori ([PR#441](https://github.com/key4hep/k4geo/pull/441))
  - port allegro-04 (minimal) changes back to v03 and remove v04

* 2025-03-12 Sanghyun Ko ([PR#440](https://github.com/key4hep/k4geo/pull/440))
  - Modified CTests for IDEA & Allegro to not to shoot particles to (0,0,1)

* 2025-03-12 BrieucF ([PR#351](https://github.com/key4hep/k4geo/pull/351))
  - Include common review comments in the PR template

* 2025-02-27 Giovanni Marchiori ([PR#439](https://github.com/key4hep/k4geo/pull/439))
  - Adds to ECalEndcap_Turbine_o1_v03_geo the code that fills the dd4rec calorimeter data (type and envelope)

* 2025-02-24 Giovanni Marchiori ([PR#438](https://github.com/key4hep/k4geo/pull/438))
  - Reverts change of HCal barrel readout introduced by mistake in https://github.com/key4hep/k4geo/pull/437

* 2025-02-23 Giovanni Marchiori ([PR#437](https://github.com/key4hep/k4geo/pull/437))
  - fixes a problem where compact files for two different detectors were using the same names for different constants

* 2025-02-23 mahmoudali2 ([PR#436](https://github.com/key4hep/k4geo/pull/436))
  - Adding general surface plugin for IDEA muon-system.
  - Fixing of the hierarchy of the detElemnt (the former one had errors when running with surfaces).
  - Now we have ~ 7256 surfaces, but 'ddPlanarDigi' still complaining that it doesn't see surfaces!

* 2025-02-22 Erich Varnes ([PR#426](https://github.com/key4hep/k4geo/pull/426))
  This PR is to support v3 of the ECal endcap turbine geometry. In this version there is flexibility to set parameters separately for the three wheels, to calibrate the response in both dimensions along the surface of a blade, and the segment the readout in both directions as well. The readout unit cells are allowed to be smaller than the calibration cells.
  
  (Note that the above changes were implemented in two steps, with the flexibility to set parameters separately introduced in "v2" of the geometry, and the remaining changes in v3. v2 was never included in a pull request.)

* 2025-02-20 Giovanni Marchiori ([PR#433](https://github.com/key4hep/k4geo/pull/433))
  - Fixes bug in calculation of radial distance between digitised hit position and inner radius of layer for inclined ecal barrel. The information is possibly used by reconstruction algorithms downstream (e.g. pandora)

* 2025-02-14 jmcarcell ([PR#434](https://github.com/key4hep/k4geo/pull/434))
  - Add LANGUAGES CXX to CMakeLists.txt to disable checks for a C compiler

* 2025-02-12 Alvaro Tolosa Delgado ([PR#431](https://github.com/key4hep/k4geo/pull/431))
  - Re-enable tests that depend on data extension of the drift chamber. This is a follow up of a previous PR, https://github.com/key4hep/k4geo/pull/400

* 2025-02-11 Giovanni Marchiori ([PR#419](https://github.com/key4hep/k4geo/pull/419))
  - Use theta rather than eta grid for muon tagger readout in ALLEGRO
  - Add detector type and calorimeter data to various detectors, for pandora PFA
  - Add simple cylinder geometry v02 with calo data and possibility to use it for the endcap with both endcaps created within the same volume
  - Add simple cylinder geometry v03 which introduces longitudinal segmentation

* 2025-02-07 jmcarcell ([PR#415](https://github.com/key4hep/k4geo/pull/415))
  - Change lcgeo_DIR to k4geo_DIR for all files in the repository

* 2025-02-04 sss ([PR#423](https://github.com/key4hep/k4geo/pull/423))
  - Fix cell position calculation for the Allegro ECal.

* 2025-01-30 Sanghyun Ko ([PR#367](https://github.com/key4hep/k4geo/pull/367))
  - Skip the scintillation optical generation for the (monolithic) DRC and migrate to the `Geant4RegexSensitivesConstruction`.

* 2025-01-29 Victor Schwan ([PR#424](https://github.com/key4hep/k4geo/pull/424))
  - collect the files of the individual MDI versions in a “global” MDI version configuration file in the folder of this version
  - avoid redundant specification of the MDI version number in each file name of each component of the MDI version (single source of truth)

* 2025-01-28 Victor Schwan ([PR#425](https://github.com/key4hep/k4geo/pull/425))
  - add Default Material to ILD as well
  - minor change to improve the readability in config of ILD_FCCee_v02.xml

* 2025-01-24 Alvaro Tolosa Delgado ([PR#422](https://github.com/key4hep/k4geo/pull/422))
  - Add Default material to CLD o4 v5 material file, so it can be used with CAD models. 
  - Add simulation test of 1 event for CLD o4 v5

* 2025-01-24 jmcarcell ([PR#400](https://github.com/key4hep/k4geo/pull/400))
  - Improve the check for DCH_info.h from DD4hep. It was added in version 1.29 so if the version is below that don't use it.

* 2025-01-23 Lorenzo Pezzotti ([PR#413](https://github.com/key4hep/k4geo/pull/413))
  - A new subdetector (DRBarrelTubes_o1_v01) is included. It is the dual-readout barrel calorimeter built exploiting the Hidra-like capillary-tubes-based technology. Together with the endcap just merged, it completes the hadronic calorimeter of the IDEA_o2 concept.
  - The IDEA_o2 muon system dimensions are adjusted to avoid overlaps with the barrel dual-readout calorimeter: (BarrelFirstLayerRadius=4630mm, BarrelLength=9260mm, EndcapFirstLayerZOffset=4630mm, EndcapLayersOuterRadius=5450mm)
  - The custom DRTubes sensitive detector action was extended to be used in both the endcap (DREndcapTubes) and the barrel (DRBarrelTubes) subdetectors.

* 2025-01-22 scott snyder ([PR#418](https://github.com/key4hep/k4geo/pull/418))
  - Fix potential uses of dangling temporaries.

* 2025-01-03 jmcarcell ([PR#417](https://github.com/key4hep/k4geo/pull/417))
  - Do not set the color flow from MCParticles, needed after https://github.com/key4hep/EDM4hep/pull/389

* 2024-12-19 Alvaro Tolosa Delgado ([PR#416](https://github.com/key4hep/k4geo/pull/416))
  - The drift chamber cell shape has been updated to a composite design, improving navigation speed in simulations by factor 10. To revert to the original twisted tube shape, add the optional `useG4TT` tag in the detector section.

* 2024-12-15 Jana ([PR#412](https://github.com/key4hep/k4geo/pull/412))
  - CLD_o4_v05: Fixing LAr_ECalBarrel.xml constants name (was changed by mistake in #361)

* 2024-12-06 Lorenzo Pezzotti ([PR#411](https://github.com/key4hep/k4geo/pull/411))
  - A new subdetector for the dual-readout capillary-tubes-based endcap calorimeter is included as discussed [here](https://indico.cern.ch/event/1439207/contributions/6056623/attachments/2903299/5092292/lopezzot_fccsim_2472024.pdf)
  - As INFN recently changed the baseline of the IDEA detector, a new IDEA_o2 XML file is created accordingly. It includes the new `DREndcapTubes` subdetector, changes the solenoid inner radius and removes the preshower. XML files for IDEA_o2 dimensions and materials are included as well. IDEA_o2 will be completed with the new crystal em calorimeter and the dual-readout capillary-tubes-based barrel calorimeter.
  - A new example `example/SteeringFile_IDEA_o2_v01.py` is included to run a simulation with the IDEA_o2 concept.
  - A new SDAction (`DRTubesSDAction`) is included and applied to the `DREndcapTubes` subdetector via `regexSD` in `example/SteeringFile_IDEA_o2_v01.py`, this speeds up the event rate w.r.t. previous dual-readout calorimetry simulations and reduces the memory fooprint as explained [here](https://indico.cern.ch/event/1457476/contributions/6136107/attachments/2929865/5144647/lopezzot_fccsim_1892024.pdf)

* 2024-12-03 Archil Durglishvili ([PR#409](https://github.com/key4hep/k4geo/pull/409))
  - HCalBarrel_TileCal_v03.xml and HCalEndcaps_ThreeParts_TileCal_v03.xml files have been created with new readouts: phi-theta and phi-row segmentation
  - Two new segmentation classes are created for HCal: FCCSWHCalPhiTheta_k4geo and FCCSWHCalPhiRow_k4geo
    - New segmentation classes can calculate the global positions of the cells and determine the cell neighbours for TopoClustering

* 2024-11-14 Swathi Sasikumar ([PR#361](https://github.com/key4hep/k4geo/pull/361))
  - This updates version ECalBarrel_NobelLiquid_InclinedTrapezoid_o1_v02_geo.cpp and ECalBarrel_NobelLiquid_InclinedTrapezoid_o1_v03_geo.cpp with layered calorimeter data for PandoraPFA as given in ECalBarrel_NobelLiquid_InclinedTrapezoid_o1_v01_geo.cpp
  - Update on the compact file LAr_ECalBarrel.xml based on the latest parameters and names...

* 2024-11-08 Giovanni Marchiori ([PR#407](https://github.com/key4hep/k4geo/pull/407))
  - Make material in 1st layer of noble-liquid ECAL absorber configurable. Set by default to G10 rather than LAr.

* 2024-10-25 Daniel Jeans ([PR#404](https://github.com/key4hep/k4geo/pull/404))
  - fix some inconsistent definitions of materials used in ILD simulation (sum of fractions != 1.000)
  - adjust composition of RPC gas to R134a:0.93 + CO2:0.05 + SF6:0.02 [exchange 5% isoButane for CO2, as suggested by G. Grenier]

* 2024-10-25 Victor Schwan ([PR#392](https://github.com/key4hep/k4geo/pull/392))
  - formatting and refactorization
  - `material_plots_2D.py`:
    - `angleDef` is a choice parameter + `thetaRad` choice added
    - `outputDir` option added
  - `material_scan_2d.py`:
    - argParser added (`compactFile`,`outputFileBase`, `angleDef`)
    - `outputDir` option added

* 2024-10-17 tmadlener ([PR#402](https://github.com/key4hep/k4geo/pull/402))
  - Fix the CI badge and add a zenodo badge to the README

* 2024-10-17 Daniel Jeans ([PR#390](https://github.com/key4hep/k4geo/pull/390))
  - Fix severe overlaps in ILD_l5_v11 inner tracker (outer barrel layer)
  - Introduce new experimental models ILD_FCCee_v01 and ILD_FCCee_v02 under FCCee/ILD_FCCee directory
  -- MDI from FCCee common MDI area
  -- vertex & lumi cal identical to CLD_o2_v07
  -- for ILD_FCCee_v01: inner tracker now based on TrackerBarrel_o1_v06 (adapted from CLD_o2_v07), TPC inner radius slightly increased 329 -> 365 mm to avoid MDI envelope
  -- for ILD_FCCee_v02: inner tracker identical to CLD_o2_v07, TPC inner radius significantly increased 329 -> 701 mm to accommodate larger inner tracker
  -- no ECAL ring, LHCAL, beamcal
  -- other subdetectors same as ILD_l5_v02
  -- some moving/renaming of definitions in ILD_common_v02 was necessary

* 2024-10-16 BrieucF ([PR#403](https://github.com/key4hep/k4geo/pull/403))
  - Lower IDEA drift chamber verbosity

# v00-21-00

* 2024-10-02 mahmoudali2 ([PR#401](https://github.com/key4hep/k4geo/pull/401))
  - IDEA_o1_v03: Moving IDEA muon system starting point in both barrel (rmin) and endcap (z-offset) 30 mm, backwards to avoid overlapping with Dual readout calorimeter.
    - The change including also updating chamber names to distinguish between muon chambers and pre shower chamber, since they are using the same uURWELL chamber.
    - Change dimensions of solenoid/endplate and move the pre-shower to avoid overlaps with dual readout calorimeter

* 2024-10-02 michaela mlynarikova ([PR#395](https://github.com/key4hep/k4geo/pull/395))
  - fix printout messages in HCal Tile barrel and Endcap three parts detector builder

* 2024-10-01 tmadlener ([PR#398](https://github.com/key4hep/k4geo/pull/398))
  - Switch the `pre-commit` action to run in a Key4hep environment
  - Add `ruff` formatting to `pre-commit`
  - Fix a few python2 style `print` statements
  - Fix format of all python files

* 2024-10-01 mahmoudali2 ([PR#397](https://github.com/key4hep/k4geo/pull/397))
  - Generalizing the muon system builder to adopt pre-shower description, the changes include:
       - Making the variable names more general, not specific only for muon system.
       - Changing the detector side's volume thicknesses in case there is only one chambers row in the side (like the pre-shower case), and it's in general it is the case for any -almost- circular shape for the detector. 
       - Disallowing the overlap rotation in case of single chamber side.
  - Overall, the changes make the builder more general to adopt different cases with different structures (polyhedron & cylinder of chambers).

* 2024-10-01 Thomas Madlener ([PR#394](https://github.com/key4hep/k4geo/pull/394))
  - Improve readability of README for FCCee MDI

* 2024-09-25 armin.ilg ([PR#396](https://github.com/key4hep/k4geo/pull/396))
  - No more warnings in silicon wrapper
  - Improvements in vertex builder printouts
  - Adding all materials of beam pipe also to material_plots_2D.py, as without having the beam pipe enabled, the vertex material budget estimation will fail.
  - Changed paths to .stl files in vertex to use https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/IDEA/IDEA_o1_v03/STL_FILES/ (still commented out due to overlaps)

* 2024-09-19 Armin Fehr ([PR#363](https://github.com/key4hep/k4geo/pull/363))
  - Update of IDEA vertex, with the ability to use the ultra-light vertex concept in-situ.
  - No overlaps in all of vertex and silicon wrapper (not including the DDCAD imported vertex support and cooling cones yet), more performant silicon wrapper (only silicon wrapper barrel sensors are simplified)

* 2024-09-18 Erich Varnes ([PR#379](https://github.com/key4hep/k4geo/pull/379))
  * ECalEndcap_Turbine_o1_v01_geo: Fix issues with printout (to allow verbosity to be controlled from run script).  
  * Add ECalEndcap_Turbine_o1_v02_geo of the "turbine" endcap geometry: which allows for more flexibility than v01 (for example, one can set different blade angles for the three wheels in v02).  As v02 is still a work in progress, the default xml points to v01.

* 2024-09-16 JEANS Daniel Thomelin Dietrich ([PR#388](https://github.com/key4hep/k4geo/pull/388))
  - For ILD models only: apply the same step limits as defined for the tracker ("Tracker_limits", currently 5mm) inside the beampipe volume and MDI region. This is important for tracking of low momentum particles  (eg beamstrahlung pairs) especially in non-uniform fields. Should have no noticeable effect in other situations.

* 2024-09-12 Andre Sailer ([PR#391](https://github.com/key4hep/k4geo/pull/391))
  - CLD_o2_v07: change LumiCal envelopes from boolean of boolean to assembly, fixes #306, speeds up overlap check (of LumiCal only) with /geometry/test/resolution 300000 down to 13s instead of 3m10s

* 2024-09-10 jmcarcell ([PR#387](https://github.com/key4hep/k4geo/pull/387))
  - Use the Key4hepConfig flag to set the standard, compiler flags and rpath magic.

* 2024-09-03 jmcarcell ([PR#386](https://github.com/key4hep/k4geo/pull/386))
  - Do not link against podio and EDM4hep dictionaries. Introduced in https://github.com/key4hep/k4geo/pull/346, I think it's never necessary to link to the dictionaries.

* 2024-09-02 Andre Sailer ([PR#385](https://github.com/key4hep/k4geo/pull/385))
  - FieldMapXYZ, FieldMapBrBz: adapt to variable rename from DD4hep, fix "OverlayedField   ERROR add: Attempt to add an unknown field type.", fixes #384

* 2024-08-29 Andre Sailer ([PR#383](https://github.com/key4hep/k4geo/pull/383))
  - CLD_o2_v07: fix overlaps related to the LumiCal, slight correction in the position of the envelopes and passive material. Fixes #376

* 2024-08-28 Leonhard Reichenbach ([PR#369](https://github.com/key4hep/k4geo/pull/369))
  - Added TrackerBarrel_o1_v06 using a stave assembly instead of directly placing the sensors into the layers
  - Added CLD_o2_v07 using the new TrackerBarrel_o1_v06

* 2024-08-28 michaela mlynarikova ([PR#350](https://github.com/key4hep/k4geo/pull/350))
  - added new HCalEndcaps_ThreeParts_TileCal_v02.xml: migrated to use FCCSWGridPhiTheta_k4geo; fixed radial dimensions, so the outer radius of all three parts is the same; renamed nModules to nsegments for number of layers in the second cylinder; uses geometry CaloThreePartsEndcap_o1_v02
  
  - added new HCalBarrel_TileCal_v02.xml which uses geometry HCalTileBarrel_o1_v01
  
  - updated ALLEGRO_o1_v03.xml to include HCalBarrel_TileCal_v02.xml and HCalEndcaps_ThreeParts_TileCal_v02.xml
  
  - added new HCalThreePartsEndcap_o1_v02_geo.cpp: added extension to store the radii of each radial layer as well as dimensions of cells. These will be used by the CellPositionsHCalPhiThetaSegTool to calculate the radii of each layer. Improved code readability and variables naming 
  
  - updated HCalTileBarrel_o1_v01_geo.cpp: added extension to store the radii of each radial layer as well as dimensions of cells. These will be used by the CellPositionsHCalPhiThetaSegTool to calculate the radii of each layer. Improved code readability and variables naming

* 2024-08-22 Victor Schwan ([PR#378](https://github.com/key4hep/k4geo/pull/378))
  - 2nd SIT barrel layer ID was corrected for `ILD_l5_v11`; the error stemmed from out-commenting 2nd out of 3 layers without adjusting hard-coded layer IDs

* 2024-08-20 BrieucF ([PR#372](https://github.com/key4hep/k4geo/pull/372))
  - [FCCeeMDI] Use absolute path to import CAD files

* 2024-08-09 jmcarcell ([PR#374](https://github.com/key4hep/k4geo/pull/374))
  - Fix a few compiler warnings

* 2024-08-09 Erich Varnes ([PR#373](https://github.com/key4hep/k4geo/pull/373))
  FCCSWEndcapTurbine_k4geo segmentation: Correct the y position for cells in the endcap on the -z side of the detector (a minus sign is needed since this detector is a mirrored copy of the +z side).

* 2024-08-09 jmcarcell ([PR#368](https://github.com/key4hep/k4geo/pull/368))
  - Clean up includes

* 2024-08-09 Alvaro Tolosa Delgado ([PR#365](https://github.com/key4hep/k4geo/pull/365))
  - IDEA_o1_v03: Dual Readout Calorimeter (DRC) is not loaded by default
     - Added Test for IDEA with DRC

* 2024-08-09 jmcarcell ([PR#353](https://github.com/key4hep/k4geo/pull/353))
  - muonSystemMuRWELL_o1_v01.cpp: Use + std::to_string to append to a string, instead of adding an integer to a string (introduced in https://github.com/key4hep/k4geo/pull/322). Adding a string and an integer cuts the string by as many characters as the value of the integer.

* 2024-08-08 BrieucF ([PR#371](https://github.com/key4hep/k4geo/pull/371))
  - Put the stl files for CAD beampipe, downloaded with cmake, at the right place

* 2024-08-06 Alvaro Tolosa Delgado ([PR#359](https://github.com/key4hep/k4geo/pull/359))
  - New CMake option `INSTALL_BEAMPIPE_STL_FILES` can be used to download the STL (CAD model) beam pipe files from the web EOS

* 2024-08-06 Sungwon Kim ([PR#346](https://github.com/key4hep/k4geo/pull/346))
  - Add DRC geometry construction code under `detector/calorimeter/dual-readout` directory
  - Add .xml compact files under `FCCee/IDEA/compact/IDEA_o1_v03` directory
  - Add custom SD action, output file, fast simulation (boosting optical photon transportation) for Monolithic fiber DRC under `plugin` directory
  - Fixed CMakeLists to compile all above

* 2024-07-30 Leonhard Reichenbach ([PR#362](https://github.com/key4hep/k4geo/pull/362))
  - TrackerEndcap_o2_v06_geo: Fixed endcap radius calculation for the event display (CED), and only the event display, fixes #355

* 2024-07-30 jmcarcell ([PR#360](https://github.com/key4hep/k4geo/pull/360))
  - Add aliases for the detectorCommon and detectorSegmentations libraries

* 2024-07-22 Giovanni Marchiori ([PR#357](https://github.com/key4hep/k4geo/pull/357))
  - [ALLEGRO_o1_v03 ECAL barrel] Get number of modules passed to readout from constant defined before in xml

* 2024-07-22 Erich Varnes ([PR#347](https://github.com/key4hep/k4geo/pull/347))
  - Added a new driver, `ECalEndcap_Turbine_o1_v01`, to build a Noble Liquid ECAL endcap with inclined blades (aka turbine geometry)
  - Added a new segmentation (`FCCSWEndcapTurbine_k4geo`) for the Noble Liquid ECAL endcap turbine geometry
  - Replaced the ALLEGRO_o1_v03 ECAL endcap made of disks perpendicular to the z axis by the turbine geometry built with `ECalEndcap_Turbine_o1_v01`

* 2024-07-19 aciarma ([PR#344](https://github.com/key4hep/k4geo/pull/344))
  - added `k4geo/FCCee/MDI` folder 
  - put the shape based beampipe in `MDI_o1_v00`
  - prepared `k4geo/FCCee/MDI/compact/MDI_o1_v01/` which will contain the CAD beampipe
  - modified the main compact files of `ALLEGRO_o1_v03` and `IDEA_o1_v03` to include the centralized beampipe and prepare them for the CAD ones
  - removed `HOMAbsorbers.xml` from ALLEGRO_o1_v03 since they are not needed with the low impedance beam pipe.

* 2024-07-16 jmcarcell ([PR#354](https://github.com/key4hep/k4geo/pull/354))
  - Rename the lcgeoTests folder to test

* 2024-07-08 mahmoudali2 ([PR#322](https://github.com/key4hep/k4geo/pull/322))
  - Define the first draft of the detailed muon system, which depend on mosaics of 50 * 50 cm^2 mRWELL chambers.
  - Define a suitable XML for the new detailed version.
  - Describe µRWELL materials.
  - Add the parameters of the muon system into the full IDEA implementation.

* 2024-07-04 Giovanni Marchiori ([PR#349](https://github.com/key4hep/k4geo/pull/349))
  - fix detector type in ALLEGRO v03 ECAL calibration scripts

* 2024-07-04 jmcarcell ([PR#348](https://github.com/key4hep/k4geo/pull/348))
  - CMake: fix printout for missing header file, by printing the actual missing file

* 2024-06-28 tmadlener ([PR#343](https://github.com/key4hep/k4geo/pull/343))
  - Make the TPC have detector ID 4 and use a consistent cellID for all tracking detectors in `ILD_l5_v11` in order to make tracking code run again

* 2024-06-19 jmcarcell ([PR#342](https://github.com/key4hep/k4geo/pull/342))
  - Remove a warning by deleting an unused string

* 2024-06-12 jmcarcell ([PR#340](https://github.com/key4hep/k4geo/pull/340))
  - Remove a CentOS7 workflow using the CLIC nightlies

* 2024-06-10 Frank Gaede ([PR#333](https://github.com/key4hep/k4geo/pull/333))
  - fix CED event display for CLIC like detectors using TrackerEndcap_o2_v0x_geo 
         -  fix nPetals in ZDiskPetalsData (for CEDViewer) to use nmodules (e.g. 48 ) rather than nrings
         -  store the number of rings in  `ZDiskPetalsData::sensorsPerPetal`

* 2024-06-06 BrieucF ([PR#339](https://github.com/key4hep/k4geo/pull/339))
  - [FCCee-ALLEGRO_o1_v03] Vertex detector and drift chamber updated to the last IDEA version
  - [FCCee-ALLEGRO_o1_v03] Added solenoidal and MDI magnetic fields
  - [FCCee-ALLEGRO_o1_v03]  Removed ALLEGRO_o1_v03_ecalonly.xml and ALLEGRO_o1_v03_trackeronly.xml to ease maintenance (they can be obtained by commenting out sub-detectors)

* 2024-05-30 BrieucF ([PR#335](https://github.com/key4hep/k4geo/pull/335))
  - Fix for the IDEA_o1_v03 solenoid position

* 2024-05-15 Alvaro Tolosa Delgado ([PR#330](https://github.com/key4hep/k4geo/pull/330))
  - Implementation of IDEA drift chamber, o1_v02. It is based on an original description. The code was redesign to be light. Standalone overlap ctest added.

* 2024-05-06 Zhibo WU ([PR#334](https://github.com/key4hep/k4geo/pull/334))
  - xtalk_neighbors_moduleThetaMergedSegmentation.cpp: change the loop variable itField from int to size_t, in order to remove a compilation warning.

* 2024-04-18 Zhibo WU ([PR#331](https://github.com/key4hep/k4geo/pull/331))
  - Add new functions related to the crosstalk neighbour finding for ALLEGRO ECAL barrel cells.

* 2024-04-16 Giovanni Marchiori ([PR#332](https://github.com/key4hep/k4geo/pull/332))
  - New version of ALLEGRO detector o1_v03 with ecal barrel with 11 layers with cell corners projective along phi. No changes to the other sub detectors.

* 2024-03-24 BrieucF ([PR#326](https://github.com/key4hep/k4geo/pull/326))
  - CLD_o2_v06: add new detector model

* 2024-03-21 Brieuc Francois ([PR#329](https://github.com/key4hep/k4geo/pull/329))
  - Add a description of CLD_o4_v05 in the FCCee README

* 2024-03-21 BrieucF ([PR#327](https://github.com/key4hep/k4geo/pull/327))
  - Add a solenoid and magnetic fields for IDEA. The solenoid has the right material budget (0.75 X0) and spacial extent but its internals should be revised by ultra-thin solenoid designers.
  - Add the endplate lead absorbers for IDEA

* 2024-03-06 Anna Zaborowska ([PR#328](https://github.com/key4hep/k4geo/pull/328))
  - Make LCIO an optional dependency. If LCIO is not found, some detectors (trackers) will not be built.

* 2024-02-25 jmcarcell ([PR#324](https://github.com/key4hep/k4geo/pull/324))
  - Remove the old key4hep build workflow since there is a newer one that builds for all the supported OSes

* 2024-02-25 jmcarcell ([PR#302](https://github.com/key4hep/k4geo/pull/302))
  - Clean up unused variables

# v00-20-00

* 2024-02-22 BrieucF ([PR#323](https://github.com/key4hep/k4geo/pull/323))
  - Increase IDEA test timeout

* 2024-02-21 Giovanni Marchiori ([PR#318](https://github.com/key4hep/k4geo/pull/318))
  * ALLEGRO: Update some xml config files

* 2024-02-16 Andre Sailer ([PR#320](https://github.com/key4hep/k4geo/pull/320))
  - CMake: Add shim for older LCIO versions that do not have LCIO::lcio yet
  - LinearSortingPolicy: adapt check for existing extension to silence error message from DD4hep

* 2024-02-09 Andre Sailer ([PR#317](https://github.com/key4hep/k4geo/pull/317))
  - ALLEGRO: EcalBarrel: fix the length of the Z extent for the calorimeter data extension for reconstruction

* 2024-02-06 BrieucF ([PR#315](https://github.com/key4hep/k4geo/pull/315))
  - [ALLEGRO_o1_v02] Replace the simplified drift chamber with the detailed one
  - [TESTS] Add a test for ALLEGRO_o1_v02
  - [ALLEGRO_o1_v02] File and variable renaming plus propagation of overlap fix to upstream and calibration xml's

* 2024-02-06 Juraj Smiesko ([PR#309](https://github.com/key4hep/k4geo/pull/309))
  - CMake: Add printing of "Found k4geo" message in downstream projects picking up k4geoConfig.cmake

* 2024-01-18 jmcarcell ([PR#310](https://github.com/key4hep/k4geo/pull/310))
  - Change LCIO::LCIO to LCIO::lcio. The target provided by LCIO is LCIO::lcio; DD4hep provides LCIO::LCIO so building without DD4hep and using LCIO::LCIO doesn't seem to work.

* 2024-01-09 Jana ([PR#308](https://github.com/key4hep/k4geo/pull/308))
  - fixing the overlaps in ALLEGRO Ecal barrel geometry by increasing the LAr bath volume.

* 2023-12-18 Swathi Sasikumar ([PR#312](https://github.com/key4hep/k4geo/pull/312))
  - CLD_o4_v05: The name of the inclined calorimeter was given as `EmCaloBarrelInclined` before. This is now changed to the latest naming convention as `ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v01`
  - ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v01: The units for cell sizes has been added to the CaloLayerData ensure that they are given correctly and not mistaken due to confusion in units. now they are represented in mm.

* 2023-12-14 BrieucF ([PR#307](https://github.com/key4hep/k4geo/pull/307))
  - Added DriftChamber_o1_v01, a first version of the detailed IDEA drift chamber to enable further technical developments (digitization, tracking, PFlow, ...)

* 2023-12-13 Giovanni Marchiori ([PR#304](https://github.com/key4hep/k4geo/pull/304))
  - Code changes related to k4geo migration of FCC segmentation classes and related utilities
  - Also moves HCal readout for Allegro from eta-based to theta-based
  - Basically a rebasing of https://github.com/key4hep/k4geo/pull/296 after https://github.com/key4hep/k4geo/pull/298 was merged

* 2023-12-12 BrieucF ([PR#311](https://github.com/key4hep/k4geo/pull/311))
  - Changed the sensitive detector type of the ALLEGO simplified drift chamber to be able to run with ddsim

# v00-19-00

* 2023-11-22 alvarotd ([PR#289](https://github.com/key4hep/k4geo/pull/289))
  * CLD_o3_v05 is created. This option evolves from option 2 version 05. The main characteristics with respect the mother version are:
  	- ARC subdetectors included between the trackers and the ECal. 
  	- Trackers are shrunk by 20 cm in radial and Z direction to accommodate the ARC. The space between layers is reduced proportionally, and the angular coverage is preserved.

* 2023-11-16 Leonhard Reichenbach ([PR#300](https://github.com/key4hep/k4geo/pull/300))
  - Add pre-commit check-xml test to validate compact files

* 2023-11-16 Daniel Jeans ([PR#299](https://github.com/key4hep/k4geo/pull/299))
  - fix the xml issues relating to ILD models from https://github.com/key4hep/k4geo/pull/297#issuecomment-1803985971

* 2023-11-16 Andre Sailer ([PR#298](https://github.com/key4hep/k4geo/pull/298))
  - CMake: The lcgeo project is renamed to k4geo, this also changes library names
  - CMake: add the creation of a k4geoConfig.cmake file

* 2023-11-15 SwathiSasikumar ([PR#286](https://github.com/key4hep/k4geo/pull/286))
  ECalBarrelInclined_geo: Implemented filling of the LayeredCalorimeterData extension, which now allows reconstruction of clusters by PandoraPFA. The detector geometry compact files used for the reconstruction are added in```CLD_o4_v05```.

* 2023-11-08 Leonhard Reichenbach ([PR#297](https://github.com/key4hep/k4geo/pull/297))
  - Thank you for writing the text to appear in the release notes. It will show up
    exactly as it appears between the two bold lines
  - ...

* 2023-11-03 Alvaro Tolosa Delgado ([PR#293](https://github.com/key4hep/k4geo/pull/293))
  - This PR is the equivalent to [PR56](https://github.com/HEP-FCC/FCCDetectors/pull/56/) in FCCDetectors
    - New segmentation imported from FCCDetectors, refactored as `FCCSWGridModuleThetaMerged_k4geo`
    - New version of ALLEGRO_o1_v02, evolving the ECAL barrel  to include the new segmentation

* 2023-10-04 Daniel Jeans ([PR#288](https://github.com/key4hep/k4geo/pull/288))
  - Update to SServices00 to allow individual components to be included/excluded
  -- default behaviour unchanged
  
  - add new model ILD_l5_v10 : ILD model for ILC
  -- inner silicon (VTX, SIT, FTD) replaced by FCCee_o2_v02-inspired design
  
  - add new model ILD_l5_v11 : ILD model for FCCee
  -- inner silicon replaced by FCCee_o2_v02-inspired design
  -- beampipe, MDI from FCCee_o2_v02

* 2023-10-04 Alvaro Tolosa Delgado ([PR#282](https://github.com/key4hep/k4geo/pull/282))
  * Migration of IDEA and ALLEGRO detectors from FCCDetectors repository, becoming each one the option 1 version 01 of the corresponding detector
  	- The compact files and detector constructors files, and the detector name therein, are renamed to avoid collision with FCCDetector components by appending the suffix `_o1_v01`
  	- IDEA compact files are copied into `./FCCee/IDEA/compact/IDEA_o1_v01`
  	- ALLEGRO compact files are copied into `./FCCee/ALLEGRO/compact/ALLEGRO_o1_v01`
  	- The detector constructors are placed in the corresponding directory inside `./detector`, and therein in the appropriate sub-directory;
  		- `calorimeter`: CaloEndcapDiscs_o1_v01_geo.cpp, HCalThreePartsEndcap_o1_v01_geo.cpp, HCalTileBarrel_o1_v01_geo.cpp
  		- `other`: SimpleCylinder_geo_o1_v01.cpp
  		- `tracker`: parametrised_DriftChamber_o1_v01.cpp
  	- DD4hep segmentations from FCCDetectors were also migrated into the directory `detectorSegmentations`, and renamed by appending the suffix `_k4geo`

* 2023-10-03 Andre Sailer ([PR#291](https://github.com/key4hep/k4geo/pull/291))
  * CLD_o2_v05: Beampipe_o4_v05: correct the material for the ConicalChamber from Beryllium to Albemet162

* 2023-10-03 Armin Fehr ([PR#273](https://github.com/key4hep/k4geo/pull/273))
  - Adding the IDEA_o1_v02 detector in key4hep, for the start just the vertex detector (vertex barrel + outer barrel) and the machine components, some place holder detectors otherwise.

* 2023-08-25 Daniel Jeans ([PR#285](https://github.com/key4hep/k4geo/pull/285))
  - Add a simple calorimeter stack, based on CaloPrototype_v02, which uses MultiSegmentation
     - e.g. can be used for a stack with alternating sensitive layer types

* 2023-08-23 Andre Sailer ([PR#287](https://github.com/key4hep/k4geo/pull/287))
  - Examples: make arcsim example always create DD4hep ROOT output

* 2023-08-23 Andre Sailer ([PR#284](https://github.com/key4hep/k4geo/pull/284))
  - CMake: Install detector XML files in share/k4geo instead of additional compact

* 2023-08-23 Andre Sailer ([PR#280](https://github.com/key4hep/k4geo/pull/280))
  - FCCee_o2_v03: Copy of FCCDetectors FCC_o2_v03, and renaming of XML files to make diff with existing models easier, this model is obsolete
  - FCCee_o2_v04: Copy of FCCDetectors FCC_o2_v04, and renaming of XML files to make diff easier
  
  - CLD_o2_v05: New implementation of the CLD detector with small beampipe design according to latest standard design; modified vertex detector to fit into the beampipe constraints (reduced length of all barrel layers, based on FCCee_o2_v02 from k4Geo); fixed overlaps in Inner and Outer Tracker

* 2023-07-14 jmcarcell ([PR#281](https://github.com/key4hep/k4geo/pull/281))
  - Rename `CMAKE_{SOURCE,BIN}_DIR` to `PROJECT_{SOURCE,BIN}_DIR`

* 2023-07-05 Andre Sailer ([PR#279](https://github.com/key4hep/k4geo/pull/279))
  - CLD Detector models moved from FCCee to FCCee/CLD, for #267

* 2023-06-27 jmcarcell ([PR#278](https://github.com/key4hep/k4geo/pull/278))
  - Remove a couple of unneeded lines in CMakeLists.txt

* 2023-06-27 jmcarcell ([PR#277](https://github.com/key4hep/k4geo/pull/277))
  - Fix setting RPATH for the installed binaries

* 2023-06-20 BrieucF ([PR#276](https://github.com/key4hep/k4geo/pull/276))
  - Migrate the noble liquid ECAL barrel C++ builder from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/main/Detector/DetFCChhECalInclined/src/ECalBarrelInclined_geo.cpp) together with its documentation

* 2023-06-13 tmadlener ([PR#275](https://github.com/key4hep/k4geo/pull/275))
  - Move to `dd4hep::CellID` for cellIDs, after [AIDASoft/DD4hep#1125](https://github.com/AIDASoft/DD4hep/pull/1125)

# v00-18-01

* 2023-05-26 jmcarcell ([PR#274](https://github.com/key4hep/k4geo/pull/274))
  - Change eps_min to .001 in the example steering file

* 2023-04-10 Frank Gaede ([PR#266](https://github.com/key4hep/k4geo/pull/266))
  -  add definitions of regions to the SEcal06_hybrid barrel and endcap models of ILD:
       - `EcalBarrelRegion` and `EcalEndcapRegion`
       - these can be used for fast simulation (ML generation)

* 2023-01-13 Daniel Jeans ([PR#263](https://github.com/key4hep/k4geo/pull/263))
  - deal with case when previousStep is not yet defined (could cause crash)
  - add timing information to the lowPt hits

# v00-18

* 2022-12-02 Daniel Jeans ([PR#261](https://github.com/ilcsoft/lcgeo/pull/261))
  - bug fix: no longer carry over energy to next event (sometimes affected lowpt hits)
  - address one FIXME: create highpt hits even if steps have not passed the padrow centre
  - some slight code cleanup and reorganisation to hopefully make it clearer to read

* 2022-11-07 Thomas Madlener ([PR#260](https://github.com/ilcsoft/lcgeo/pull/260))
  - Remove gcc8 based workflow because the underlying nightly builds are no longer available
  - Update github actions to latest available versions.

* 2022-11-07 Daniel Jeans ([PR#250](https://github.com/ilcsoft/lcgeo/pull/250))
  - ILD_l5[_o1,2,3,4]_v09 models with CLIC-inspired all silicon outer tracker in place of TPC+SET
  - otherwise identical to ILD_l5_v02 model

# v00-17

* 2022-10-13 Andre Sailer ([PR#259](https://github.com/ilcsoft/lcgeo/pull/259))
  - CLIC_o3_v15: New detector model with corrected TrackerEndcapSupport with the same thicknesses in forward and backward directions
  - FCCee_o1_v05: New detector model with corrected TrackerEndcapSupport, with the same thicknesses in forward and backward directions
  - FCCee_o2_v02: New detector model with corrected TrackerEndcapSupport, with the same thicknesses in forward and backward directions

* 2022-08-19 Andre Sailer ([PR#258](https://github.com/ilcsoft/lcgeo/pull/258))
  - TrackerEndcapSupport_o1_v02: Fix the thickness of the whole layer volume, was twice as thick as needed to be

* 2022-08-18 Nazar Bartosik ([PR#257](https://github.com/ilcsoft/lcgeo/pull/257))
  - Add TrackerEndcapSupport_o1_v02, contributed by @bartosik-hep. Fixes an issue that the reflected sides in TrackerEndcapSupport_o1_v01 were twice as thick as the non-reflected sides.

* 2022-07-28 Andre Sailer ([PR#256](https://github.com/ilcsoft/lcgeo/pull/256))
  - LinearSortingPolicy: change return value to 1 for success, adapt to recent changes in DD4hep that took the actual return value from the plugin, (https://github.com/AIDASoft/DD4hep/pull/936)

* 2022-07-27 Andre Sailer ([PR#255](https://github.com/ilcsoft/lcgeo/pull/255))
  - TPCSDAction: adapt to new const G4step interface in dd4hep 1.21
  - CaloPreShowerSDAction: adapt to new const G4step interface in dd4hep 1.21

# v00-16-08

* 2022-06-14 Valentin Volkl ([PR#253](https://github.com/iLCSoft/lcgeo/pull/253))
  - Move SiD_o2_v04 beampipe constants to global list to fix an error in key4hep builds

* 2022-06-09 Dan Protopopescu ([PR#252](https://github.com/iLCSoft/lcgeo/pull/252))
  - Added SiD_o2_v04, which is an updated version of o2_v03, containing a few fixes, among which
      - correct size and position of ECal layer 0 via new driver ECalBarrel_o2_v04_geo.cpp
      - new beam pipe by Chris Potter
      - removed brass HCal option

* 2022-03-19 Valentin Volkl ([PR#251](https://github.com/iLCSoft/lcgeo/pull/251))
  - Fix more hyphens in xml comments

* 2022-03-09 Andre Sailer ([PR#249](https://github.com/iLCSoft/lcgeo/pull/249))
  - Rebrand LCGEO as Lepton Collider GEOmetry, Fix #248

# v00-16-07

* 2021-11-05 scott snyder ([PR#246](https://github.com/iLCSoft/lcgeo/pull/246))
  - Fixed XML comment syntax in SiD_o2_v03 XML files.

* 2021-09-01 Andre Sailer ([PR#247](https://github.com/iLCSoft/lcgeo/pull/247))
  - CI: Run against LCG_100, clang10, gcc10 and LCG_99python2 gcc8
  - CMake: restructure main CMake file, more use of targets
  - Remove requirement for Boost, wasn't actually used for some time (still needed by DD4hep)
  - GenericEndcapCalo: move setting of sensitive type because of error in SID simulation with newer ddsim

* 2020-09-21 Valentin Volkl ([PR#242](https://github.com/iLCSoft/lcgeo/pull/242))
  - Add standard cpp/dd4hep .gitignore

* 2020-09-18 vvolkl ([PR#244](https://github.com/iLCSoft/lcgeo/pull/244))
  - Fix print statement for python3 compatibiltiy

* 2020-05-25 Valentin Volkl ([PR#243](https://github.com/iLCSoft/lcgeo/pull/243))
  - Add INSTALL_COMPACT_FILES option (default: OFF) to copy compact files to install area

# v00-16-06

* 2020-03-02 Remi Ete ([PR#241](https://github.com/iLCSoft/lcgeo/pull/241))
  - `Yoke05_Barrel` driver:
       - Removed slices DetElement construction
       - Moved layer DetElement construction to external loop: fixes the DetElement hierarchy.
  This fixes the access to the DetElement given a CellID in our reconstruction. 
  Related issue: https://github.com/iLCSoft/ILDConfig/issues/88

# v00-16-05

* 2019-12-11 Ete Remi ([PR#240](https://github.com/iLCSoft/lcgeo/pull/240))
  - Fixed steering file for unit test. Get units from g4units.

* 2019-12-11 Daniel Jeans ([PR#239](https://github.com/iLCSoft/lcgeo/pull/239))
  - SEcal06 driver: add absorber layer thickness to layer information (this had been missing)

# v00-16-04

* 2019-09-11 Andre Sailer ([PR#237](https://github.com/iLCSoft/lcgeo/pull/237))
  - CLIC_o4_v14: New CLIC detector model with DECal, based on CLIC_o3_v14. Ported from https://github.com/robbie-bosley/CLIC_DECAL

* 2019-05-23 Emilia Leogrande ([PR#236](https://github.com/iLCSoft/lcgeo/pull/236))
  - FCCee_o2_v01: new detector model with smaller beampipe radius (10 mm instead of 15 mm) and closer vertex subdetector

* 2019-02-20 Emilia Leogrande ([PR#235](https://github.com/iLCSoft/lcgeo/pull/235))
  - FCCee_o1_v04: Implemented changes the new sorting policy
     - Analogous of #234

* 2019-02-19 Andre Sailer ([PR#234](https://github.com/iLCSoft/lcgeo/pull/234))
  - Plugin: lcgeo_LinearSortingPolicy : new plugin to set the sorting policy for surfaces identified by placement path based on their z position and a linear function, requires AIDASoft/DD4hep#486
  - CLIC_o3_v14: add use of lcgeo_LinearSortingPolicy for TrackerEndcaps
  - TrackerEndcapSupport_o1_v01: attach support volumes to different DetElements so we can attach different extensions to them

* 2018-11-07 Marko Petric ([PR#232](https://github.com/iLCSoft/lcgeo/pull/232))
  - Add vis attribute for better visualisation of VertexEndcap_o1_v06

* 2018-11-01 Marko Petric ([PR#231](https://github.com/iLCSoft/lcgeo/pull/231))
  - Add vis attribute for better visualisation of `TrackerEndcap_o1_v02`

* 2018-08-23 Oleksandr Viazlo ([PR#229](https://github.com/iLCSoft/lcgeo/pull/229))
  - new FCCee_o1_v04 detector model: 
    - update on MDI and shielding from Anna
    - merge Materials_v01.xml and materials.xml into one file and remove unused items

* 2018-08-09 Andre Sailer ([PR#228](https://github.com/iLCSoft/lcgeo/pull/228))
  - ddsim program moved to DD4hep, AIDAsoft/DD4hep#420

* 2018-07-05 Andre Sailer ([PR#225](https://github.com/iLCSoft/lcgeo/pull/225))
  - Tests: Using other cmake variables to be able to pick up the lcgeo tests from outside of lcgeo

# v00-16-02

* 2018-06-25 Frank Gaede ([PR#223](https://github.com/ilcsoft/lcgeo/pull/223))
  - fix travis CI: use wget --no-check-certificate

* 2018-06-25 Daniel Jeans ([PR#222](https://github.com/ilcsoft/lcgeo/pull/222))
  - apply anti-DID field map to small ILD models. (assume same anti-DID field as for large models)

* 2018-05-23 Dan Protopopescu ([PR#216](https://github.com/ilcsoft/lcgeo/pull/216))
  - Added separate GlobalForwardCaloReadoutID for LumiCal and BeamCal and reverted the XMLs to use the DD4hep drivers from o2_v02
  - Deleted now unused driver CylCalEndcap_o1_v01_geo.cpp
  - Reverted to 'systemID=20' in plugin definition of the Barrel ECal

* 2018-06-21 Daniel Jeans ([PR#220](https://github.com/ilcsoft/lcgeo/pull/220))
  - QD0 and QDEX1A magnet strengths for 1 TeV
  - new models
    ILD_(sl)5_v07 : 1 TeV fwd magnets, solenoid field map, no anti-DID
    ILD_(sl)5_v08 : 1 TeV fwd magnets, solenoid field map + with anti-DID

* 2018-06-01 Marko Petric ([PR#219](https://github.com/ilcsoft/lcgeo/pull/219))
  - Accommodate AIDASoft/DD4hep#397 - `DD4hep.py` was dropped in favor of `dd4hep.py`

# v00-16-01

* 2018-05-22 Frank Gaede ([PR#218](https://github.com/ilcsoft/lcgeo/pull/218))
  - ILD:  add 4T solenoid field map for ILD_s5_o?_v03 models

* 2018-05-18 Emilia Leogrande ([PR#217](https://github.com/ilcsoft/lcgeo/pull/217))
  - Vertex_o4_v05.xml: doubled support material in the vertex (barrel and endcap)
  - This model is made for testing purposes only

* 2018-04-26 Oleksandr Viazlo ([PR#215](https://github.com/ilcsoft/lcgeo/pull/215))
  - new FCCee_o1_v03 detector model:
    - extend ECAL endcap
    - shrink HCAL ring
    - reduce magnetic field in barrel yoke from 1.5 T to 1.0 T
    - fix in HCAL layer layout (swap steel absorber with steel wall of cassette)

# v00-16

* 2018-02-22 Dan Protopopescu ([PR#201](https://github.com/iLCSoft/lcgeo/pull/201))
  - Removed some unused parameters

* 2018-02-22 Daniel Jeans ([PR#200](https://github.com/iLCSoft/lcgeo/pull/200))
  - to ensure that consistent maximum step lengths are used in all models, moved their definition to common directory

* 2018-02-27 Daniel Jeans ([PR#199](https://github.com/iLCSoft/lcgeo/pull/199))
  - implement downstream magnets in ILD models:
      - thickened beampipe in regions of final focus magnets
      - should reproduce the Mokka models
      - currently verifying positions with accel. experts, so details may still change

* 2018-02-26 Andre Sailer ([PR#203](https://github.com/iLCSoft/lcgeo/pull/203))
  - DDSim: allow passing arbitrary arguments to add_argument, needed for #196 , purely behind the scenes

* 2018-02-28 Daniel Jeans ([PR#204](https://github.com/iLCSoft/lcgeo/pull/204))
  - implement step limits in SEcal06 driver (SEcal06_Helpers.cpp)
  - apply cal_limits to ECAL descriptions (barrel, endcap, ECring) in ILD_common_v02

* 2018-02-28 Ete Remi ([PR#196](https://github.com/iLCSoft/lcgeo/pull/196))
  - Added new command line parameters to ddsim
      - --eventParameters (list) to write event parameters to every output lcio events
      - --runNumberOffset (int) to offset the run number counter in lcio output file
      - --eventNumberOffset (int) to offset the event number counter in lcio output file

* 2018-02-13 Andre Sailer ([PR#197](https://github.com/iLCSoft/lcgeo/pull/197))
  - DDSim: Print the startUp and event processing time separately in addition to total runtime

* 2017-12-01 TiborILD ([PR#186](https://github.com/iLCSoft/lcgeo/pull/186))
  - update ILD models w/ SDHcal:
         - change the names of original SimCalorimeterHit collections of SDHCAL to the same ones as introduced in hybrid models, i.e. HCalBarrelRPCHits, HCalEndcapRPCHits,HCalECRingRPCHits

* 2017-11-23 Dan Protopopescu ([PR#185](https://github.com/iLCSoft/lcgeo/pull/185))
  * SiD_o2_v02: Fix BeamCal segmentation (thanks to @shaojunlu), Fixes #184

* 2017-11-23 Shaojun Lu ([PR#183](https://github.com/iLCSoft/lcgeo/pull/183))
  - Update lcgeoTests for ILD to ILD_l5_v02 and ILD_s5_v02.
      - the current models for optimisation studies, will be tested during the build.

* 2018-03-16 Frank Gaede ([PR#210](https://github.com/iLCSoft/lcgeo/pull/210))
  - improve TPCSDAction
     - enable writing of TPCSpacePoint and TPCLowPtCollections
     - TPCLowPtCollections is only written if `TPCLowPtStepLimit==true`
     - add steering properties to TPCSDAction 
     - use  in ddsim:
             `SIM.action.mapActions['tpc'] = ('TPCSDAction', {'TPCLowPtCut': 10*MeV , 'TPCLowPtStepLimit' : True })`
     - use Mokka default values:
  ```
  	Control.TPCLowPtCut = CLHEP::MeV ;
  	Control.TPCLowPtStepLimit = false ;
  	Control.TPCLowPtMaxHitSeparation = 5. * CLHEP::mm ;
  ```
    -  properly access layerID through IDDescriptor

* 2018-03-14 Daniel Jeans ([PR#209](https://github.com/iLCSoft/lcgeo/pull/209))
  - new driver for tube of square cross-section (BoxSupport_o1_v01_geo.cpp)
  - new driver for yoke endcaps with square central hole (Yoke06_Endcaps.cpp)
  - adjust description of QD0 and first extraction quad, and add its cryostat
  - remove more upstream magnets
  - include forward region support tubes: around QD0, outer tube froum outside -> LHCAL, and between ecal endcap and ring
  - update yoke design to fit with this support structure
  - implement the above in the ILD_?5_v02 models

* 2017-12-12 Dan Protopopescu ([PR#188](https://github.com/iLCSoft/lcgeo/pull/188))
  Version 3 of the SiD option 2 model (not final though). Main changes and updates:
  - cleaned up all XMLs
  - using DetIDs
  - unified segmentation and cell ID encoding
  - added type_flags where missing
  - updated plugins
  - added envelopes where missing
  - stepped design Muon calorimeters use now scintillator instead of RPCs
  - merged barrel+endcap envelopes for the Muon calorimeter
  - added or updated certain drivers
  - included alternate brass HCal 
  - added common tracker barrel+endcap envelope
  - implemented some ECal barrel driver fixes
  - checked: no overlaps
  - added test_SiD_o2_v03
  
  To check:
  - Muon caloData definitions
  - ECal towers visibility attributes
  - LumiCal and BeamCal envelopes
  
  A detailed table of changes is available in
  https://www.evernote.com/l/AJ2Q2yTuLwRGVqYpSNB8HRKx-iBxYQ6Ry0s

* 2018-03-05 Frank Gaede ([PR#207](https://github.com/iLCSoft/lcgeo/pull/207))
  -  activate G4StepLimiterPhysics in DDSim by default
      - needed to limit the steps in detector volumes and regions

* 2018-03-05 Dan Protopopescu ([PR#202](https://github.com/iLCSoft/lcgeo/pull/202))
  Minor corrections to SiD_o2_v03:
  - Corrected LumiCal grid size parameters
  - Added hole in the BeamCal (beam pipe parameters not used with the current driver)
  - Removed DetType_AUXILIARY flag from BeamCal
  - Removed some unused parameters

* 2018-03-28 Frank Gaede ([PR#214](https://github.com/iLCSoft/lcgeo/pull/214))
  - Fix for the removal of DDSurfaces which have been merged into DDRec 
    -  includes from `DDSurfaces` -> `DDRec`
    - namespace `DDSurfaces` -> `dd4hep::rec`

* 2017-12-11 Oleksandr Viazlo ([PR#193](https://github.com/iLCSoft/lcgeo/pull/193))
  - new FCCee_o1_v02 detector model
  - Material budget of VTX was increased for 50%
  - change of LumiCal_cell_size from 1.0 to 1.81 mm, change of LumiCall offset
  - add test for FCCee_o1_v02 model, remove test for FCCee_o1_v01

* 2017-12-04 Andre Sailer ([PR#189](https://github.com/iLCSoft/lcgeo/pull/189))
  - Add CLIC_o3_v14 with fixed LumiCal segmentation

* 2017-12-04 Frank Gaede ([PR#187](https://github.com/iLCSoft/lcgeo/pull/187))
  - Set DDSim defaults to Geant4 defaults and retain `largest_step = 10 * m` as recommended in the Geant4 manual (~ detector size). This improves simulation time for single muons approximately a factor 2 and for single pi+ for 25%

* 2018-03-27 Ete Remi ([PR#213](https://github.com/iLCSoft/lcgeo/pull/213))
  - DD4hepSimulation:
     - Fixed division by zero when no event has been processed, Fixes #208

* 2018-03-22 Oleksandr Viazlo ([PR#211](https://github.com/iLCSoft/lcgeo/pull/211))
  - create FCCee_dev detector model (extended ECAL endcap, shrinked HCAL ring, reduced magnetic field in yoke from 1.5 T to 1.0 T)
  - add test for the model

* 2018-03-23 Daniel Jeans ([PR#212](https://github.com/iLCSoft/lcgeo/pull/212))
  - writeAllILDCompactDescriptions.py: python script to create compact detector descriptions with different variations: large/small, ideal/realistic fields at 250/500 GeV, with or without anti-DID fields, different reconstruction options
  - compact descriptions all kept in one directory ILD_sl5_v02
  - each model is linked individually
  - some models are not yet functional: in particular small models with realistic fields, for which no field maps are yet available. [ ie ILD_s5_v03, v04, v05, v06 ]
  - ILD_?5_v02 models and their options should be completely unchanged

# v00-15-03

* 2017-11-22 Frank Gaede ([PR#182](https://github.com/ilcsoft/lcgeo/pull/182))
  - update the documentation for the ILD detector models in `./ILD/compact/Readme.md`

* 2017-11-22 TiborILD ([PR#180](https://github.com/ilcsoft/lcgeo/pull/180))
  - Bug fixing and modifications in ILD_l2_v02 model 
      - Removed 1st comment line from Hcal_EndcapRing_SD_v01.xml  preventing to appear this detector in the geometry
      -  alligned Hcal_EndcapRing_SD to the back envelope and replaced the material of the remaining free space at the start of this detector by air
      - reduced the thickness of HcalSD_back_plate  from 15->10 mm; gain of additional layer in Endcaps; now at (47/48)
      - ILD_l2_v02  stands for SDHCAL technology in VIDEAU geometry
      - some modifications in Hcal_Endcaps_SD_v01, Hcal_Endcaps_SD_v02
  - Added  ILD_l6_v02 for  SDHCAL technology in TESLA geometry
  - Added  ILD_s6_v02 a small version w/SDHCAL in TESLA geomtry

# v00-15-02

* 2017-11-20 Daniel Jeans ([PR#181](https://github.com/ilcsoft/lcgeo/pull/181))
  - bug fix for ILD_l/s5 models
     - fix calculation of layer positions in dd4hep::rec::LayeredCalorimeterData in SEcal06_Helpers

# v00-15-01

* 2017-11-10 Shaojun Lu ([PR#179](https://github.com/iLCSoft/lcgeo/pull/179))
  - update HcalEndcapRing in ILD models
       - Fix HcalEndcapRing first layer placement Z offset.
       - Update HcalEndcapRing alignment to the back of the envelope.

# v00-15

* 2017-11-02 Shaojun Lu ([PR#177](https://github.com/ilcsoft/lcgeo/pull/177))
  - Added one independent parameter "Hcal_endcap_lateral_structure_thickness"
      - to make Hcal endcap more flexible.
      - tower: |5mm|2.5mm|360mm|2.5mm|5mm|
      - tower: |LateralStructure|AirGap|Active(HBU)|AirGap|LateralStructure|
      - user may update them from compact file.

* 2017-11-02 Shaojun Lu ([PR#176](https://github.com/ilcsoft/lcgeo/pull/176))
  - Fix ILD SDHCal endcaps segmentation offset.
      - use half of SDHCal_cell_size offset to make sure that cells equally distributed in the sensitive area from the center to the left, to the right, to the top and to the bottom.
      - to avoid the cell to be placed at (0,0) in xy-plane, then generate the fragmentation cell at the boundary.

* 2017-10-17 Oleksandr Viazlo ([PR#174](https://github.com/ilcsoft/lcgeo/pull/174))
  - add new FCCee_o1_v01 detector model (based only on DD4hep and lcgeo drivers)
  - add test for this model

* 2017-10-27 Shaojun Lu ([PR#175](https://github.com/ilcsoft/lcgeo/pull/175))
  - Fix ILD HCAL endcaps drivers to use half of Hcal_lateral_structure_thickness for each side.
      - and half of Hcal_endcap_layer_air_gap size for each side too.
      - these two parameters could be modified inside the compact files "hcal_defs.xml"

# v00-14

* 2017-08-09 Marko Petric ([PR#152](https://github.com/iLCSoft/lcgeo/pull/152))
  - New detector model CLIC_o3_v13
  - CLIC_o3_v13: Set the cell size in HCalEndcap, HCalRing and YokeEndcap to be read from main compact file
  - CLIC_o3_v13: Extend the support in the vertex of layer 4 to the edge of the sensitive

* 2017-10-10 Ete Remi ([PR#172](https://github.com/iLCSoft/lcgeo/pull/172))
  - Split ILD ecal segmentation selection into two values: value0 and value1 while parsing geometry for ILD_l/s5_* models. 
  - Changed both drivers and compact xml files.

* 2017-10-12 Frank Gaede ([PR#173](https://github.com/iLCSoft/lcgeo/pull/173))
  - use unbounded surfaces in the ILD TPC
      - assign surfaces only to the forward half of the TPC
      - as they are unbounded they will also be 'visible' in the backward half

* 2017-08-25 Shaojun Lu ([PR#156](https://github.com/iLCSoft/lcgeo/pull/156))
  - Fix the ILD_l4/s4_v02 Hybird Hcal RPC Anode and Cathode FloatGlass place order.
    - For the Hybird Hcal Endcaps and EndcapRing, they have been fixed inside compact file.
      - "SHcalSc04_EndcapRing_v01.xml" for both ILD_l4_v02 and ILD_s4_v02
      - "SHcalSc04_Endcaps_v01_LARGE.xml" for ILD_l4_v02
      - "SHcalSc04_Endcaps_v01_SMALL.xml" for ILD_s4_v02
    - For the Hybird Hcal Barrel, they have been fixed in both the driver and compact file
      - "SHcalSc04_Barrel_v04.xml" for both ILD_l4_v02 and ILD_s4_v02.
      - "SHcalSc04_Barrel_v04.cpp" the slices place order has been fixed to start from IP.
   - All the compact files are inside the "ILD_common_v02" folder, which are used by ILD_l5/s5_v02, too.
      - So both ILD_l5_v02 and ILD_s5_v02 have the updeted Hybird Hcal, too.

* 2017-07-28 TiborILD ([PR#150](https://github.com/iLCSoft/lcgeo/pull/150))
  - added inner box cut off for Hcal_Endcaps_SD_v01
  - added SDHCAL Barrel module with TESLA geometry

* 2017-08-30 TiborILD ([PR#157](https://github.com/iLCSoft/lcgeo/pull/157))
  - update of ILD Hcal endcap readout segmentation
          - in order to unify the readouts for all Hcal detectors  the CartesianGridXZ readout of SHCalsc04_Endcaps_v01 was changed to CartesianGridXY readout 
         - affects models ILD_l/s4_v02 and ILD_l/s5_v02

* 2017-08-30 Jan Strube ([PR#153](https://github.com/iLCSoft/lcgeo/pull/153))
  - allow to change the name of the MCParticle input collection in ddsim
          - using parameter 'mcParticleCollectionName'

* 2017-08-30 Jan Strube ([PR#153](https://github.com/iLCSoft/lcgeo/pull/153))
  /

* 2017-08-02 TiborILD ([PR#151](https://github.com/iLCSoft/lcgeo/pull/151))
  - modifications to comply with a new xml scheme in ILD_common_v02

* 2017-09-26 Daniel Jeans ([PR#167](https://github.com/iLCSoft/lcgeo/pull/167))
  - adjust parameters of hybrid ECAL model: reduce scintillator 2 -> 1.5 mm, add 0.5mm air gap
  - take hybrid ECAL segmentation from parameters (not hard-coded numbers)

* 2017-08-14 Andre Sailer ([PR#154](https://github.com/iLCSoft/lcgeo/pull/154))
  - The CommandLine information was no longer filled in the lcio runheader, it is now filled again

* 2017-09-29 Ete Remi ([PR#169](https://github.com/iLCSoft/lcgeo/pull/169))
  - Ecal (Hcal) segmentation selection variable moved from ecal(hcal)_defs.xml files to top level detector xml file
  - Added 12 new detectors for different reconstruction combination : 
    - ILD_l/s4_o1_v02 : SiWEcal + AHcal (as for standard ILD_l/s4_v02)
    - ILD_l/s4_o2_v02 : SiWEcal + SDHCAL
    - ILD_l/s5_o1_v02 : SiWEcal + AHcal (as for standard ILD_l/s5_v02)
    - ILD_l/s5_o2_v02 : SiWEcal + SDHCAL
    - ILD_l/s5_o3_v02 : ScEcal + AHcal
    - ILD_l/s5_o4_v02 : ScEcal + SDHCAL
  - One can use simulated files from i.e ILD_l4_v02 and reconstruct them with ILD_l4_o1_v02 or ILD_l4_o2_v02 models

* 2017-09-29 Frank Gaede ([PR#168](https://github.com/iLCSoft/lcgeo/pull/168))
  - adapt to changes in dd4hep::BitField64 in https://github.com/AIDASoft/DD4hep/pull/238

* 2017-07-27 Shaojun Lu ([PR#148](https://github.com/iLCSoft/lcgeo/pull/148))
  - Update ILD Documentation tools for ILD_l4/s4_v02.
        - Following the DDG4 development, and update ILD Documentation tool "extractParameters.py"
        - Following the ILD_l4/s4_v02 setup, and clean unused parameter in "ILD_envelope_parameters.txt" for the envelope parameters list.

* 2017-08-23 Daniel Jeans ([PR#155](https://github.com/iLCSoft/lcgeo/pull/155))
  - In case of multi-readout setup, give the two (or more) sensitive layers within a single calorimeter layer the same layer index

* 2017-09-14 Shaojun Lu ([PR#166](https://github.com/iLCSoft/lcgeo/pull/166))
  - Fix ILD_l4_v02 Hybird model collections key_value.
      - follow the update of the slice placement and correct the collections key_value.
      - the segmentation and collection should use the same key_value for both SDHCAL and  AHcal.

* 2017-09-06 Shaojun Lu ([PR#160](https://github.com/iLCSoft/lcgeo/pull/160))
  - Fix ILD SEcal05_ECRing for filling dd4hep::rec::LayeredCalorimeterData.
    - For the reconstruction, we fill the LayeredCalorimeterData at runtime.
    - Both side ECRings have identical layer structure, fill the one from module_num==0, and used for both.

* 2017-09-06 Shaojun Lu ([PR#159](https://github.com/iLCSoft/lcgeo/pull/159))
  - Fix ILD_l4_v02 target readout slice number for AHCAL.
    - after re-order the Hybird HCAL layers, we should also update Hcal_readout_segmentation_slice_barrel to "3" as "Hcal_readout_segmentation_slice_endcap" now.

* 2017-09-12 Daniel Jeans ([PR#165](https://github.com/iLCSoft/lcgeo/pull/165))
  Fix the filling of dd4hep::rec::LayeredCalorimeterData::Layer data for SEcal05, SEcal06 models:
  - layer.distance now defined as the distance to the *start* of a calo layer (at the start of the absorber, except in the case of the 1st layer, in which it's the start of the ECAL).
  Some additional fixes to SEcal05_ECRing: 
  - use correct absorber thickness in stack transitions (previously not guaranteed to be correct)
  - use CarbonFibre when calculating material properties for LayeredCalorimeterData (previously using air for structural material)
  - cosmetic changes: indentations, commented-out code, ...

* 2017-09-12 Shaojun Lu ([PR#164](https://github.com/iLCSoft/lcgeo/pull/164))
  - For backward consistency, update all ILD AHCAL Barrel slice placement.
      - update all ILD AHCAL Barrel drivers to place slice in the same consistent way.
      - update the slice order in all the ILD AHCAL Barrel compact files too.

* 2017-09-12 Dan Protopopescu ([PR#163](https://github.com/iLCSoft/lcgeo/pull/163))
  - Option 3 of the SiD model, with Cu absorber plates in the HCal, instead of Steel235

* 2017-09-12 Shaojun Lu ([PR#162](https://github.com/iLCSoft/lcgeo/pull/162))
  -  Fix ILD_l4_v02 Barrel MultiSegmentation slice number for different technology.
     - After change the order of slice, the segmentation slice should be updated too.
     - ILD_s4_v02 use the same compact file, should be fine.

* 2017-10-06 Andre Sailer ([PR#171](https://github.com/iLCSoft/lcgeo/pull/171))
  - Drop unused and no longer existing header includes AidaSoft/DD4hep#241

* 2017-10-05 Frank Gaede ([PR#170](https://github.com/iLCSoft/lcgeo/pull/170))
  - fill `tpcData->zMinReadout` in `TPC10_geo.cpp`
    - describing the cathode thickness

* 2017-09-11 Frank Gaede ([PR#161](https://github.com/iLCSoft/lcgeo/pull/161))
  - fix spelling of ddsim parameter:
       - rename mcParticleCollectionName to MCParticleCollectionName

# v00-13-04

* 2017-07-21 Shaojun Lu ([PR#147](https://github.com/iLCSoft/lcgeo/pull/147))
  - Fine tuning FTD envelope and SIT Cone cables to fix overlap that was not detected on Z-axis for ILD_l4_v02.

* 2017-07-21 Shaojun Lu ([PR#145](https://github.com/iLCSoft/lcgeo/pull/145))
  - Update lcgeoTests for ILD.
    - ILD_l4/s4_v02 are the current optimisation models.
    - Remove test of ILD_l1_v01, and added tests for ILD_l4/s4_v02.

* 2017-07-21 Shaojun Lu ([PR#144](https://github.com/iLCSoft/lcgeo/pull/144))
  - Fix VXD04 overlap by changing offset_phi - used in all ILD models

* 2017-07-20 Marko Petric ([PR#146](https://github.com/iLCSoft/lcgeo/pull/146))
  - Change starting position of inner tracker barrel to remove gaps and overlaps in CLIC_o3_v12

# v00-13-03

* 2017-07-19 Daniel Jeans ([PR#143](https://github.com/iLCSoft/lcgeo/pull/143))
  - introduce SEcal06 drivers, which can deal with multi-readout (requires at least DD4hep v01-01-01 for updated MegatileLayerGridXY class). Should be completely back-compatible with SEcal05 drivers.
  - defined hybrid ECAL, with both silicon and scintillator readout in each layer
  - defined large and small ILD models using such a hybrid ECAL: ILD_l5_v02 and ILD_s5_v02

* 2017-07-19 Daniel Jeans ([PR#142](https://github.com/iLCSoft/lcgeo/pull/142))
  - Hcal_readout_segmentation_slice_barrel/endcap paramters were specified in both the ILD_?v_v02.xml and hcal_defs.xml files, to inconsistent values. the values in ILD_?v_v02.xml were correct, and were used. The values in hcal_defs.xml were wrong, but were not used. 
  I have corrected values in hcal_defs.xml, removed duplicate definitions in ILD_?v_v02.xml.

* 2017-07-19 Marko Petric ([PR#141](https://github.com/iLCSoft/lcgeo/pull/141))
  - New model CLIC_o3_v12 that fixes an overlap that was considered a false positive in outer tracker barrel due to shift in z0
  - In file `OuterTracker_o2_v06_01.xml` L6 value of `OuterTracker_Barrel_half_length` change from 1264mm -> 1264.2mm

* 2017-07-19 Remi Ete ([PR#140](https://github.com/iLCSoft/lcgeo/pull/140))
  * Modified DD4hepSimulation.parseOptions() to receive argv argument, allowing user scripts to pass custom command line arguments

* 2017-07-14 Remi Ete ([PR#139](https://github.com/iLCSoft/lcgeo/pull/139))
  - CLIC_o3_v11: move HCal_cell_size constant to main file to avoid duplicate definition
  - CLIC_o3_v11: BeamCal: set correct size of incoming beampipe, only changes of absorber shape, reduce Beampipe wall thickness in beamcal by 0.02 mm to avoid overlaps

# v00-13-02

* 2017-07-11 Frank Gaede ([PR#132](https://github.com/iLCSoft/lcgeo/pull/132))
  - add material surfaces for ILD _*_v02 models to SIT and VXD cables in (SService00)
        - use a cylinder for VXD cables ( rather than outwards pointing cone)

* 2017-07-12 Frank Gaede ([PR#136](https://github.com/iLCSoft/lcgeo/pull/136))
  - fix beamcal overlaps for ILD_*_v02 models
      - add new version BeamCal_o1_v02 that uses assemblies in first two levels 
         of gemoetry hierarchy
      - use this for all ILD_*_v02 models

* 2017-07-12 Frank Gaede ([PR#135](https://github.com/iLCSoft/lcgeo/pull/135))
  - fix overlap in ILD TPC w/ outer cathode grip ring
      - reduce cathode grip ring height fromm 20mm to 18mm
      - introduce symmetrical inner and outer service areas
        of 18.1 mm
      - results in 220(163) layers for ILD_l* (ILDs*) model

* 2017-07-12 Daniel Jeans ([PR#134](https://github.com/iLCSoft/lcgeo/pull/134))
  changes to SEcal05_ECRing.cpp:
  - include first layer in LayeredCalorimeterData (if not preshower)
  - take layering within slab from compact description to ensure correct treatment (removing some assumptions from driver): should now be identical layering to ecal endcaps from SEcal05_Endcaps

* 2017-07-12 Daniel Jeans ([PR#134](https://github.com/iLCSoft/lcgeo/pull/134))
  changes to SEcal05_ECRing.cpp: 
  
      * include first layer in LayeredCalorimeterData (if not preshower) 
      * take layering within slab from compact description to ensure correct treatment (removing some assumptions from driver): should now be identical layering to ecal endcaps from SEcal05_Endcaps

* 2017-07-12 Shaojun Lu ([PR#133](https://github.com/iLCSoft/lcgeo/pull/133))
  - Fix ILD_*_v02 OverLap
    - Add one more env_safety to avoid YokeEndcaps and HcalEndcaps overlap.
    - To adjust the FTD envelope to fix the SIT cable cone overlap.
    - Added one more gap_thickness to fix the YokeEndcap extruded by the last layer.

# v00-13-01

* 2017-07-06 Frank Gaede ([PR#131](https://github.com/iLCSoft/lcgeo/pull/131))
  - fix the pixel SIT for ILD_*_v02 models
        - add new driver SIT_Simple_Pixel_geo.cpp
        - use sit_simple_pixel_sensors_01.xml

* 2017-07-06 Shaojun Lu ([PR#130](https://github.com/iLCSoft/lcgeo/pull/130))
  - Update ILD services driver to access the ILD compact file "env_safety", and apply it in the services driver to avoid services overlap with sub-detector envelope .

# v00-13

* 2017-07-06 Frank Gaede ([PR#129](https://github.com/iLCSoft/lcgeo/pull/129))
  - ILD_*_v02 models: remove stereo angle in SIT ( move to pixel readout)

* 2017-07-06 Daniel Jeans ([PR#128](https://github.com/iLCSoft/lcgeo/pull/128))
  - many tweaks to ECAL dimensions in ILD_*_V02 models, referring to technical design document from H.Videau et al.
  - a couple of small bug fixes to ECAL drivers (previously ignored request for extra CF thickness in front face, if no preshower was set)
  - adjusted TPC dimensions and barrel-endcap gap to accommodate new ECAL

* 2017-07-05 Shaojun Lu ([PR#127](https://github.com/iLCSoft/lcgeo/pull/127))
  - Implement ILD HcalEndcap FrontEnd readout electronics.
  - Implement it into both ILD_common_v01 and ILD_common_v02.
  - ILD models l4/s4 will use the same services as l1/s1.

* 2017-07-05 Shaojun Lu ([PR#126](https://github.com/iLCSoft/lcgeo/pull/126))
  - Update for ILD services cables
    - Implement VXD cable cone services.
    - Implement the VXD cable services into ILD_common_v01 and ILD_common_v02.
    - Comment out the services in ILD_o1_v05.
    - Remove the ILD services reference code (not necessary).
    - The SServices00 has been updated, and used by ILD.
    - Update the ILD VXD cable configuration by expert.

* 2017-07-05 Daniel Jeans ([PR#125](https://github.com/iLCSoft/lcgeo/pull/125))
  - previously several ECAL sub-layers were combined into an averaged material to simplify simulation.
  This requires to redefine this material every time thicknesses change (quite often at the moment), which was not done consistently. 
  For the time being, separate the components to ensure correct material description. 
  (Consider recombining when sub-layer thicknesses are fixed and stable.)

* 2017-07-04 Daniel Jeans ([PR#124](https://github.com/iLCSoft/lcgeo/pull/124))
  - now use DD4hep_Beampipe_o1_v01 driver (gets rid of "deprecated" warning when using Beampipe_o1_v01)

* 2017-07-04 Shaojun Lu ([PR#123](https://github.com/iLCSoft/lcgeo/pull/123))
  - Update ILD_o1_v05 compact files for the services.
    - Change the FTD envelope shape to place SIT cable.
    - Add missing service parameters.

* 2017-07-04 Daniel Jeans ([PR#122](https://github.com/iLCSoft/lcgeo/pull/122))
  - updates to TPC in ILD_*_v02 models from Dimitra Tsionou
  - increased material in outer field cage
  - increased material in cathode plane
  - reduced thickness of outer field cage by 5mm (to accommodate thicker ECAL: final numbers to be confirmed)

* 2017-07-03 Daniel Jeans ([PR#121](https://github.com/iLCSoft/lcgeo/pull/121))
  - in ILD models, material of TPC cooling pipes is accounted for in TPC driver. They were also added in sservices00 driver, causing double counting. Solution: turn them off in SServices00.

# v00-12

* 2017-06-12 Daniel Jeans ([PR#110](https://github.com/iLCSoft/lcgeo/pull/110))
  - First "development" set of ILD "v02" models, with thicker ECAL and reduced TPC outer radius: Large and Small models have been implemented in "1", "2", and "4" technology options.
  - new dimensions not yet frozen, so will probably need some adjustment
  - re-organisation of code for v2 models: readout definition moved together with detector definitions (previously defined in top-level ILD_*.xml).
  - TPC services (all cooling pipes except one) removed for v02 models. Awaiting input from TPC experts on correct position/size of this pipe.
  - some change to definition of ECAL rails: move to currently favoured 2-rail design. Some update/cleanup of SServices00 driver. This also affects v01 models.

* 2017-06-28 Shaojun Lu ([PR#117](https://github.com/iLCSoft/lcgeo/pull/117))
  - Implement SITCables services for ILD detector models
  - Update the FTD envelope shape to get free space to place SITCables.
  - Add SITCables services parameters into XML files for ILD_common_v01.
  - Add SITCables services parameters into XML files for ILD_common_v02.

* 2017-06-29 Andre Sailer ([PR#119](https://github.com/iLCSoft/lcgeo/pull/119))
  - Split libraries into one without linking against Geant4 and one with Geant4 dependency. Avoid loading of Geant4 libraries in reconstruction

* 2017-06-29 Shaojun Lu ([PR#118](https://github.com/iLCSoft/lcgeo/pull/118))
  - Follow the TUBE update, and use the TUBE envelope radius to fix the services SitDisk overlaps. 
  - Define the SitDisk Rmin with TUBE_IPOuterBulge_end_envradius instead of TUBE_IPOuterBulge_end_radius.
  - Implement HcalBarrel_EndcapServices which also include the services coming from  Ecal and TPC.
  - Implement HcalBarrel_EndcapServices parameters into both ILD_comon_v01 and ILD_common_v02.
  - Add parameter 'TPC_cooling_nRings' in SServices00.xml to allow change wrt large and small ILD model for this moment.

* 2017-06-08 Daniel Jeans ([PR#107](https://github.com/iLCSoft/lcgeo/pull/107))
  - For ILD models: increase inner radius of FTD envelope, to avoid overlap with beampipe of finite thickness.

* 2017-06-16 Frank Gaede ([PR#112](https://github.com/iLCSoft/lcgeo/pull/112))
  - remove unused and deprecated (https://github.com/AIDASoft/DD4hep/pull/165 ) DDRec extensions LayeringExtension/SubdetectorExtensions from Ecal drivers

* 2017-06-15 Shaojun Lu ([PR#111](https://github.com/iLCSoft/lcgeo/pull/111))
  - Ported ILD Mokka class SServices00 line by line into a new class SServices00_v01 into lcgeo.
  - Used in ILDSServices.cpp to create services geometry for ILD_o1_v05.

* 2017-05-15 Marko Petric ([PR#102](https://github.com/iLCSoft/lcgeo/pull/102))
  - Introduce new model CLIC_o3_v11 and shift segmentation in HCal for half cell size

* 2017-05-15 Frank Gaede ([PR#101](https://github.com/iLCSoft/lcgeo/pull/101))
  - fix the logic for assigning sensitive slices in ILD Hcal drivers

* 2017-06-30 Frank Gaede ([PR#120](https://github.com/iLCSoft/lcgeo/pull/120))
  - update the ILD_(ls)(124)_v02 models
          - move the un-instrumented gaseaous volume in the TPC to the inner field cage 
          - cleaned up the use of plugins:
                 - CaloFace plugins are now in the SEcal04 drivers (in ILD_common_v02)
                 - add DD4hepVolumeManager (for consistency) and InstallSurfaceManager
                     to main compact files for all models

* 2017-06-20 Frank Gaede ([PR#114](https://github.com/iLCSoft/lcgeo/pull/114))
  -  replace DDSurfaces w/ dd4hep::rec

* 2017-06-20 Andre Sailer ([PR#113](https://github.com/iLCSoft/lcgeo/pull/113))
  - Adapt to changes in namespaces in DD4hep

* 2017-06-21 Marko Petric ([PR#115](https://github.com/iLCSoft/lcgeo/pull/115))
  - Remove using namespace from header

* 2017-06-26 Frank Gaede ([PR#116](https://github.com/iLCSoft/lcgeo/pull/116))
  - fix ILD_o4 models:
  - ILD_l4_v01:  
     - reverse order of RPC and scintillator  ( RPC first)
  - ILD_l/s4_v02:
     - reverse order of RPC and scintillator  ( RPC first)
     - add parameters for tracking volume and reconstruction geometry
     - add surface plugins for tracking

* 2017-04-28 Daniel Jeans ([PR#100](https://github.com/iLCSoft/lcgeo/pull/100))
  - Changes to HCAL endcap envelopes in ILD_l/s* models
          - in all ILD_l/s* models, reduce radial size of hcal endcap envelopes to leave space for cables. cryo-hcal endcap gap is controlled by parameter Hcal_endcap_cryostat_gap, set to 170 mm.
         - in ILD_?2 models, change SDHCAL endcap envelope shape to tube (same as other models)

* 2017-05-29 TiborILD ([PR#106](https://github.com/iLCSoft/lcgeo/pull/106))
  - updated ILD SDHcal drivers and geometry:
      - reverse ordering of slices in Barrel, Endcaps_v01
      - correction for z-cracks in Barrel
      - towers numbering from inner to outer radius (Endcaps_v02)
      - staves numbering 0-3 (Endcaps_v02)
      - removed not used variable and commented line (Endcaps_v02)

* 2017-05-30 luisaleperez ([PR#103](https://github.com/iLCSoft/lcgeo/pull/103))
  - Updated Field maps readers, FieldMapBrBz and FieldMapXYZ, to automatically get from the root file the maps parameters
    - Coordinates ranges, step size and ordering

* 2017-04-22 StrahinjaLukic ([PR#96](https://github.com/iLCSoft/lcgeo/pull/96))
  - Outgoing beam tube radius was reduced by 1 mm from "BeamCal_min_z-5*mm" to 5999*mm because the BeamCal inner radius was reduced as part of the adjustments to L*=4.1m
  - The protruding outer radius of QDEX1AFront link between tube sections was reduced to the outer radius of the larger tube.
  - The changes affect all ILD models that use ILD_common.

* 2017-04-22 Frank Gaede ([PR#95](https://github.com/iLCSoft/lcgeo/pull/95))
  - fix the creation of the DDRec::LayeredCalorimeterData data structures for the the 
    multisegmentation readout model LD_l4_v01
  - reverse the order of the slices in the HcalEndcaps and HCalEndcapRing for this model

* 2017-05-24 Frank Gaede ([PR#104](https://github.com/iLCSoft/lcgeo/pull/104))
  - set the correct DetType_ENDCAP for all ILD HcalEndcapRing models
        - was DetType_BARREL

* 2017-04-25 Marko Petric ([PR#97](https://github.com/iLCSoft/lcgeo/pull/97))
  - Introduce new model CLIC_o3_v10
  - Unify all encoding strings for ECal, HCal, Yoke with the introduction of GlobalCalorimeterReadoutID

* 2017-04-26 Daniel Jeans ([PR#98](https://github.com/iLCSoft/lcgeo/pull/98))
  - use correct LHCal01 detector setup in ILD_l1_v01

* 2017-06-09 Marko Petric ([PR#109](https://github.com/iLCSoft/lcgeo/pull/109))
  - Update the location of DD4hep Handle.inl include file to adapt to DD4hep v00.24

* 2017-06-09 Daniel Jeans ([PR#108](https://github.com/iLCSoft/lcgeo/pull/108))
  - improve SEcal05 drivers and compact description for ILD models
       - remove hard-coded numbers of carbon fiber layers (in alveolii, around absorber plates) in SEcal05* drivers: more flexible.
      - now parameterize in terms of total CF thickness, rather than thickness of one sheet multiplied by number of sheets: easier to understand, more flexible.
      - associated changes required in compact descriptions

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
