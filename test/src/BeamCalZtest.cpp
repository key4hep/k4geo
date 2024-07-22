#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>
#include "DD4hep/DD4hepUnits.h"

#include <string>
#include <vector>

int main (int argc, char **args) {

  if ( argc < 2 ){
    std::cout << "Usage: BeamCalZTest <compact file name>.xml\n";
    exit(0);
  }
  std::string compactFile = std::string(args[1]);

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  theDetector.fromXML( compactFile );

  const std::string detectorName(theDetector.header().name());

  const std::vector< dd4hep::DetElement > &detElements = theDetector.detectors("calorimeter", true);

  const std::string beamCalName("BeamCal");

  bool foundBeamCal = false;
  double bcalOuterZ = 0;
  for(unsigned i=0; i<detElements.size(); i++) {

    std::string detName(detElements.at(i).name());
    if ( detName.compare(beamCalName) == 0 ) {
      foundBeamCal = true;
      dd4hep::rec::LayeredCalorimeterData* bcData = detElements.at(i).extension<dd4hep::rec::LayeredCalorimeterData>() ;
      bcalOuterZ = bcData->extent[3];
    }
  }


  if (foundBeamCal) {
    std::cout << "\n" << beamCalName << " Z_outer = " << bcalOuterZ/dd4hep::mm << " mm\n\n";
    if ( detectorName.find("ILD") != std::string::npos ) {

      if (bcalOuterZ/dd4hep::mm > 3420) {
        std::cout << "ERROR: " << beamCalName << " is too far out in the z-direction.\n";
        std::cout << "ILD detector models should have " << beamCalName << " with Z_outer <= 3420 mm.\n\n";
      }
      else {
        std::cout << "OK.\n\n";
      }
    }
    else {
      std::cout << "No specific requirements are set on " << beamCalName << " z-position for " << detectorName << "\n\n";
    }
  }
  else {
    std::cout << "\nCould not find detector with name " << beamCalName << ".\n\n";
  }


  return 0;

}
