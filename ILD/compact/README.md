# ILD detector models

The following ILD detector models are available in lcgeo:



| Model         |  Description                              | Status            |
| ------------- | ----------------------------------------- |-------------------|
| ILD_l1_v01    | large DBD like model, AHcal and SiW-Ecal  |  under validation |
| ILD_s1_v01    | small DBD like model, AHcal and SiW-Ecal  |  under validation |
| ILD_l2_v01    | large DBD like model, SDHcal and SiW-Ecal |  under validation |
| ILD_s2_v01    | small DBD like model, SDHcal and SiW-Ecal |  under validation |
| ILD_l4_v01    | large DBD like model, A/SDHcal, SiW-Ecal  |  experimental     |
| ILD_o3_v05    | large DBD like model, AHcal and SciW-Ecal |  experimental     |
| ILD_o1_v05    | DBD model (Mokka) , AHcal and SiW-Ecal    |  historic         |
   



## ILD_l1_v01
 - optimisation large ILD model
 - with SiWEcal and Tesla AHCAL
 - TPC_outer_radius = 1808*mm
 - compact file: [./ILD_l1_v01/ILD_l1_v01.xml](./ILD_l1_v01/ILD_l1_v01.xml)

## ILD_s1_v01
 - optimisation small ILD model
 - with SiWEcal and Tesla AHCAL
 - TPC_outer_radius = 1460*mm
 - compact file: [./ILD_s1_v01/ILD_s1_v01.xml](./ILD_s1_v01/ILD_s1_v01.xml)

## ILD_l2_v01
 - optimisation large ILD model
 - with SiWEcal and Videau SDHcal
 - TPC_outer_radius = 1808*mm
 - compact file: [./ILD_l2_v01/ILD_l2_v01.xml](./ILD_l2_v01/ILD_l2_v01.xml)

## ILD_s2_v01
 - optimisation small ILD model
 - with SiWEcal and Videau SDHcal
 - TPC_outer_radius = 1460*mm
 - compact file: [./ILD_s2_v01/ILD_s2_v01.xml](./ILD_s2_v01/ILD_s2_v01.xml)


## ILD_l4_v01
 - large ILD model
 - with SiWEcal and generic polyhedral Hcal
	 - with scintillator and RPC readout
	 - creates two sets of hit collections
 - TPC_outer_radius = 1808*mm
 - **experimental**
 - compact file: [./ILD_l4_v01/ILD_l4_v01.xml](./ILD_l4_v01/ILD_l4_v01.xml)


## ILD_o1_v05
- DBD model
- ported from Mokka
- kept for historic reference 
- compact file: [./ILD_o1_v05/ILD_o1_v05.xml](./ILD_o1_v05/ILD_o1_v05.xml)


## ILD_o3_v05
 - large ILD model
 - with **Sci-W**-Ecal and Tesla AHCAL
 - TPC_outer_radius = 1808*mm
 - experimental, Ecal only
 - compact file: [./ILD_o3_v05/ILD_o3_v05.xml](./ILD_o3_v05/ILD_o3_v05.xml)

 
