// Test the setting of the sensitiveThickness for the DataStructs in some Tracker detectorDrivers

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DDTest.h>
#include <DD4hep/Detector.h>
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>

#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static dd4hep::DDTest test( "SensThickness" ) ;

template<typename T> void checkExtensions( dd4hep::Detector& lcdd, std::vector<std::string> detectorNames, std::vector<double> sensThickness) {

  int counter = 0;
  for (auto const& detName : detectorNames ) {
    const double thickness = sensThickness[counter]; ++counter;
    auto detElement = lcdd.detector(detName);
    auto* theExtension = detElement.extension<T>();
    const std::vector<typename T::LayerLayout>& layers= theExtension->layers;
    int layerCount = 0;
    for (auto const& layer : layers) {
      std::stringstream msg;
      msg << detName << " layer " << layerCount << ": "
	  << "layer sensitive thickness:  " << layer.thicknessSensitive/dd4hep::micrometer
	  << "  expected " << thickness/dd4hep::micrometer;
      test( fabs(layer.thicknessSensitive - thickness)  < 1e-11, msg.str() );
      ++layerCount;
    }
  }
}


int main (int argc, char **args) {

  if ( argc != 4 ){
    throw std::runtime_error( "need to provide compact file and thickness for tracker and vertex silicon thicknesses in micrometer");
  }
  std::string compactFile = std::string(args[1]);

  const double trackerThick = atof(args[2]);
  const double vertexThick = atof(args[3]);

  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
  lcdd.fromCompact( compactFile );

  std::vector<std::string> endcaps = { "InnerTrackerEndcap",
				       "OuterTrackerEndcap",
				       "VertexEndcap"};

  std::vector<std::string> barrels = {"InnerTrackerBarrel",
				      "OuterTrackerBarrel",
				      "VertexBarrel" };

  std::vector<double> sensThickness = { trackerThick*dd4hep::micrometer,
					trackerThick*dd4hep::micrometer,
					vertexThick*dd4hep::micrometer };

  checkExtensions<dd4hep::rec::ZDiskPetalsData>( lcdd, endcaps, sensThickness );
  checkExtensions<dd4hep::rec::ZPlanarData>( lcdd, barrels, sensThickness );


}
