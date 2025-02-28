The Pandora Layered Calorimeter Plugin is a generic DD4hep-based tool designed to facilitate integration and processing of Layered Calorimeter Data for different experiments irrespective of their geometries. Layered Calorimeter Data stores and provides detector material properties e.g. Radiation length, interaction length to PandoraPFA using which it calculates shower depths. In this Plugin Material Manager is used to average materials between given two points and calculate the material properties of the averaged material. This approach helps the Plugin to stay geometry-agnostic.

STEP 1:
IMPORTANT CHECK: The extension for Layered Calorimeter Data can either be inside the geometry-driver code or in the Plugin. It cannot be in the both the places. If using the plugin check step 2. 

STEP 2:
Update the input values for the Plugin in the geometry compact XML file

e.g.
<plugin name="Pandora_LayeredCalorimeterPlugin">
        <argument value="ECalBarrel"/>
        <argument value="DetType_CALORIMETER + DetType_ELECTROMAGNETIC + DetType_BARREL" />
        <argument value="EMBarrel_rmin"/> <!-- Rmin ECalBarrel_inner_radius -->
        <argument value="EMBarrel_rmax"/> <!-- Rmax ECalBarrel_outer_radius -->
        <argument value="EMBarrel_dz"/>  <!-- half_length -->
       	<argument value="LAr"/>  <!-- Material -->
	<!--Input layer heights RADIALLY. Number of layers can be manipulated-->
        <argument value="1.5*cm"/>   <!-- Layer 1 height -->
        <argument value="3.5*cm"/>   <!-- Layer 2 height -->
</plugin>

VERY IMPORTANT: The layer heights should be provided along the RADIUS. For calorimeters like the Nobel Liquid Argon, the layer heights maybe defined along the electrode. It is important that these values are converted before providing them as an input.

STEP 3: Using a simulated file one can run reconstruction using the respective steering file (Config file) e.g https://github.com/SwathiSasikumar/LAr_scripts/blob/main/CLD_LAr.py. Fine tuning for PandoraPFA parameters can be made inside the steering file. 

STEP 4: CHECK if Pandora collections exist. 

