# ALFA detector 

## ALFA\_o1\_v00 
First version of integrated ALFA geometry, including: 
- Barrel Outer Tracker: OuterTrackerBarrel\_o1\_v00, comprised of the middle (2) and outer (2) layers of the ALFA tracker, corresponding to the BAR354\_EC353\_VTX108 version.
Modules consist of 8 sensors, 21 modules along z, 2 modules along rphi.
Engineering design discussed at [FCC Week 2026](https://indico.cern.ch/event/1552126/contributions/7128854/).
- Barrel Grainita ECAL: Grainita\_ECAL\_Barrel\_v02, composed with 72(phi) * 15(theta) quasi-pointing modules. 
Each module is fullfilled with averaged material: 79% ZnWO4 + 21% heavy liquid + fibers. 
No real fiber volume in the geometry. 
Inactive material include: carbon-fiber front and backword supporting, carbon-fiber module frame. 
Readout is achieved with DD4hep::Segmentation FCCSWGridRhoPhiTheta\_k4geo, granularity 3.2 mrad (theta) * 3.2 mrad (phi) in transverse and 4 in longitudinal. 
Note: segmentation is not matched with the quasi-pointing module geometry. 






