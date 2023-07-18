!-- aciarma 17.07.2023

Engineered beam pipe model + copper cooling
CAD design courtesy of F. Fransesini

To display the beampipe:
geoDisplay standalone_CADbp.xml

To use in simulation, include in the main xml file:
<include ref="PATH.TO.K4GEO/FCCee/CAD_beampipe/compact/Beampipe_CADimport.xml"/>
<include ref="PATH.TO.K4GEO/FCCee/CAD_beampipe/compact/Beampipe_onlylegs.xml"/>



STL files stored in:
PATH.TO.K4GEO/FCCee/CAD_beampipe/CAD_beampipe

STL files present (upd. 18.07.2023):
Beampipe:		  Pipe_17032023/AlBeMet162_17032023.stl
Cooling - central:        Pipe_18072023/Copper_cent_18072023.stl
Cooling - trapezoidal:    Pipe_18072023/Copper_trap_18072023_v2_brep_bin.stl
Paraffin:       	  Pipe_18072023/Paraffin_18072023.stl
Gold layer:     	  Pipe_17032023/Gold_17032023.stl


CAD model only up to -2.5m/+2.5m from the IP, after beampipe separation
