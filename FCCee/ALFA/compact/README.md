# ALFA detector 

## ALFA_o1_v00 
First version of integrated ALFA geometry, including: 
- Barrel Outer Tracker: OuterTrackerBarrel_o1_v00, 
- Barrel Grainita ECAL: Grainita_ECAL_Barrel_v02, composed with 72(phi) * 15(theta) quasi-pointing modules. 
Each module is fullfilled with averaged material: 79% ZnWO4 + 21% heavy liquid + fibers. 
No real fiber volume in the geometry. 
Inactive material include: carbon-fiber front and backword supporting, carbon-fiber module frame. 
Readout is achieved with DD4hep::Segmentation FCCSWGridRhoPhiTheta_k4geo, granularity 3.2 mrad (theta) * 3.2 mrad (phi) in transverse and 4 in longitudinal. 
Note: segmentation is not matched with the quasi-pointing module geometry. 






