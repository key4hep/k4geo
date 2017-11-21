
# ILD detector models - Overview

The following ILD detector models are available in lcgeo ( current production models **ILD_ls5_v02** )



| Model         |  Description               | Hcal   |  Ecal   | geometry | Status            |
| ------------- | ---------------------------|--------|---------|----------|-------------------|
| ILD_l5_v02    | large simulation model     | hybrid | hybrid  | Tesla    |  validated        |
| ILD_l5_o1_v02 | large reconstruction model | sci    | Si      | Tesla    |  validated        |
| ILD_l5_o2_v02 | large reconstruction model | RPC    | Si      | Tesla    |  validated        |
| ILD_l5_o3_v02 | large reconstruction model | sci    | sci     | Tesla    |  validated        |
| ILD_l5_o4_v02 | large reconstruction model | RPC    | sci     | Tesla    |  validated        |
| ------------- | ---------------------------|--------|---------|----------|-------------------|
| ILD_s5_v02    | small simulation model     | hybrid | hybrid  | Tesla    |  validated        |
| ILD_s5_o1_v02 | small reconstruction model | sci    | Si      | Tesla    |  validated        |
| ILD_s5_o2_v02 | small reconstruction model | RPC    | Si      | Tesla    |  validated        |
| ILD_s5_o3_v02 | small reconstruction model | sci    | sci     | Tesla    |  validated        |
| ILD_s5_o4_v02 | small reconstruction model | RPC    | sci     | Tesla    |  validated        |
| ------------- | ---------------------------|--------|---------|----------|-------------------|
| ILD_l4_v02    | large simulation model     | hybrid | Si      | Tesla    |  validated        |
| ILD_l4_o1_v02 | large reconstruction model | sci    | Si      | Tesla    |  validated        |
| ILD_l4_o2_v02 | large reconstruction model | RPC    | Si      | Tesla    |  validated        |
| ------------- | ---------------------------|--------|---------|----------|-------------------|
| ILD_s4_v02    | small simulation model     | hybrid | Si      | Tesla    |  validated        |
| ILD_s4_o1_v02 | small reconstruction model | sci    | Si      | Tesla    |  validated        |
| ILD_s4_o2_v02 | small reconstruction model | RPC    | Si      | Tesla    |  validated        |
| ------------- | ---------------------------|--------|---------|----------|-------------------|
| ILD_l1_v02    | large model                | sci    | Si      | Tesla    |  under validation |
| ILD_s1_v02    | small model                | sci    | Si      | Tesla    |  under validation |
| ILD_l2_v02    | large model                | RPC    | Si      | Videau   |  under validation |
| ILD_s2_v02    | small model                | RPC    | Si      | Videau   |  under validation |
| ILD_l6_v02    | large model                | RPC    | Si      | Tesla    |  under validation |
| ILD_s6_v02    | small model                | RPC    | Si      | Tesla    |  under validation |
| ILD_o1_v05    | DBD model (Mokka)          | sci    | Si      | Tesla    |  historic         |
| ------------- | ---------------------------|--------|---------|----------|-------------------|



## Details

### ILD_l5_v02
 - large ILD model
 - Ecal 
 	 - with scintillator **and** Si readout  
 - Tesla Hcal 
	 - with scintillator **and** RPC readout
	 - creates two sets of hit collections
 - TPC_outer_radius = 1808*mm
 - compact files:
	 - [./ILD_l5_v02/ILD_l5_v02.xml](./ILD_l5_v02/ILD_l5_v02.xml)
	 - [./ILD_l5_o1_v02/ILD_l5_o1_v02.xml](./ILD_l5_o1_v02/ILD_l5_o1_v02.xml)
	 - [./ILD_l5_o2_v02/ILD_l5_o2_v02.xml](./ILD_l5_o2_v02/ILD_l5_o2_v02.xml)
	 - [./ILD_l5_o3_v02/ILD_l5_o3_v02.xml](./ILD_l5_o3_v02/ILD_l5_o3_v02.xml)
	 - [./ILD_l5_o4_v02/ILD_l5_o4_v02.xml](./ILD_l5_o4_v02/ILD_l5_o4_v02.xml)

### ILD_s5_v02
 - small ILD model
 - Ecal 
 	 - with scintillator **and** Si readout  
 - Tesla Hcal 
	 - with scintillator **and** RPC readout
	 - creates two sets of hit collections
 - TPC_outer_radius = 1460*mm
 - compact files:
	 - [./ILD_s5_v02/ILD_s5_v02.xml](./ILD_s5_v02/ILD_s5_v02.xml)
	 - [./ILD_s5_o1_v02/ILD_s5_o1_v02.xml](./ILD_s5_o1_v02/ILD_s5_o1_v02.xml)
	 - [./ILD_s5_o2_v02/ILD_s5_o2_v02.xml](./ILD_s5_o2_v02/ILD_s5_o2_v02.xml)
	 - [./ILD_s5_o3_v02/ILD_s5_o3_v02.xml](./ILD_s5_o3_v02/ILD_s5_o3_v02.xml)
	 - [./ILD_s5_o4_v02/ILD_s5_o4_v02.xml](./ILD_s5_o4_v02/ILD_s5_o4_v02.xml)


### ILD_l4_v02
 - large ILD model
 - with SiWEcal and Tesla Hcal 
	 - with scintillator **and** RPC readout
	 - creates two sets of hit collections
 - TPC_outer_radius = 1808*mm
 - compact files:
	 - [./ILD_l4_v02/ILD_l4_v02.xml](./ILD_l4_v02/ILD_l4_v02.xml)
	 - [./ILD_l4_o1_v02/ILD_l4_o1_v02.xml](./ILD_l4_o1_v02/ILD_l4_o1_v02.xml)
	 - [./ILD_l4_o2_v02/ILD_l4_o2_v02.xml](./ILD_l4_o2_v02/ILD_l4_o2_v02.xml)

### ILD_s4_v02
 - small ILD model
 - with SiWEcal and Tesla Hcal 
	 - with scintillator **and** RPC readout
	 - creates two sets of hit collections
 - TPC_outer_radius = 1460*mm
 - compact files:
	 - [./ILD_s4_v02/ILD_s4_v02.xml](./ILD_s4_v02/ILD_s4_v02.xml)
	 - [./ILD_s4_o1_v02/ILD_s4_o1_v02.xml](./ILD_s4_o1_v02/ILD_s4_o1_v02.xml)
	 - [./ILD_s4_o2_v02/ILD_s4_o2_v02.xml](./ILD_s4_o2_v02/ILD_s4_o2_v02.xml)



### ILD_l1_v02
 - optimisation large ILD model
 - with SiWEcal and Tesla AHCAL
 - TPC_outer_radius = 1808*mm
 - compact file: [./ILD_l1_v02/ILD_l1_v02.xml](./ILD_l1_v02/ILD_l1_v02.xml)

### ILD_s1_v02
 - optimisation small ILD model
 - with SiWEcal and Tesla AHCAL
 - TPC_outer_radius = 1460*mm
 - compact file: [./ILD_s1_v02/ILD_s1_v02.xml](./ILD_s1_v02/ILD_s1_v02.xml)

### ILD_l2_v02
 - optimisation large ILD model
 - with SiWEcal and Videau SDHcal
 - TPC_outer_radius = 1808*mm
 - compact file: [./ILD_l2_v02/ILD_l2_v02.xml](./ILD_l2_v02/ILD_l2_v02.xml)

### ILD_s2_v02
 - optimisation small ILD model
 - with SiWEcal and Videau SDHcal
 - TPC_outer_radius = 1460*mm
 - compact file: [./ILD_s2_v02/ILD_s2_v02.xml](./ILD_s2_v02/ILD_s2_v02.xml)

### ILD_l6_v02
 - optimisation large ILD model
 - with SiWEcal and Tesla SDHcal
 - TPC_outer_radius = 1808*mm
 - compact file: [./ILD_l6_v02/ILD_l6_v02.xml](./ILD_l6_v02/ILD_l6_v02.xml)

### ILD_s6_v02
 - optimisation small ILD model
 - with SiWEcal and Tesla SDHcal
 - TPC_outer_radius = 1460*mm
 - compact file: [./ILD_s6_v02/ILD_s6_v02.xml](./ILD_s6_v02/ILD_s6_v02.xml)



### ILD_o1_v05
- DBD model
- ported from Mokka
- kept for historic reference 
- compact file: [./ILD_o1_v05/ILD_o1_v05.xml](./ILD_o1_v05/ILD_o1_v05.xml)



 
