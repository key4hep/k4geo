//==========================================================================
// iLCSoft - linear collider geometry
//--------------------------------------------------------------------------
//
// For the licensing terms see lcgeo/LICENSE.
//
// Author     : A. Sailer
//
//==========================================================================
//
// Linear Sorting Policy Extension
//
// Adds the sorting policy variable to surface DetElements following a
// linear function, surfaces are selected depending on the placement path
// of their DetElement
//
//==========================================================================

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DD4hep/Factories.h>
#include <DD4hep/Printout.h>

#include <DDRec/DetectorData.h>
#include <DDRec/SurfaceHelper.h>

#include <map>
#include <string>
#include <tuple>
#include <vector>

using dd4hep::DetElement;
using dd4hep::PrintLevel;



namespace {

  /** Plugin for adding the SortingPolicy parameter to surface
   * 
   * The sorting policy is calculated from the z Position of the surface and a first order polynomial
   *  sp = p[2] * (zPosition - p[0]) + p[1]
   * Arguments are:
   *  - Path indentifying the DetElement
   *  - zOffset for the function
   *  - rOffset for the function
   *  - Slope for the function
   *
   * any number of these four arguments can be given. The evaluation order is in the order the
   * four-tuple is provided, the first match will be used to calculate the sorting policy
   *
   * @author A.Sailer, CERN
   */
  static long addSortingPolicy(dd4hep::Detector& description, int argc, char** argv) {
    const std::string LOG_SOURCE("SortingPolicyPlugin");

    // use a vector to keep the same order as in the arguments given
    std::vector<std::pair<std::string, std::vector<double>>> pathToLinear;

    {
      std::string path = "";
      for(int i=0; i<argc; ++i)  {
        if(i % 4 == 0){
          path = std::string(argv[i]);
          dd4hep::printout(PrintLevel::DEBUG, LOG_SOURCE, "argument[%d]: %s", i, argv[i]);
          pathToLinear.emplace_back(path, std::vector<double>(3, 0.0));
        } else {
          const int variableID = (i % 4) - 1;
          pathToLinear.back().second[variableID] = dd4hep::_toDouble(argv[i]);
          dd4hep::printout(PrintLevel::DEBUG, LOG_SOURCE, "argument[%d]: %s --> %3.5f", i, argv[i],
                           pathToLinear.back().second[variableID]);
        }
      }
    }

    auto world = description.world();

    dd4hep::rec::SurfaceHelper ds( world );
    dd4hep::rec::SurfaceList const& detSL = ds.surfaceList();
    for(dd4hep::rec::ISurface* surf: detSL){
      dd4hep::Volume volume = ((dd4hep::rec::Surface*)surf)->volume();
      auto* ddsurf  = ((dd4hep::rec::Surface*)surf);
      auto const volumeName = std::string(volume->GetName());
      if(not ddsurf->detElement().isValid()) continue;
      std::string path = ddsurf->detElement().path();
      for (auto const& pathAndValues : pathToLinear) {
        if(path.find(pathAndValues.first) != std::string::npos) {
          std::vector<double> const& parameters = pathAndValues.second;
          double zPosition = std::fabs(surf->origin()[2]);
          double rValue = parameters.at(2) * (zPosition-parameters.at(0)) + parameters.at(1);

          dd4hep::rec::DoubleParameters* para = nullptr;
          try { // use existing map, or create a new one
            para = ddsurf->detElement().extension<dd4hep::rec::DoubleParameters>(false);
            if(not para) {
              para = new dd4hep::rec::DoubleParameters;
              ddsurf->detElement().addExtension<dd4hep::rec::DoubleParameters>(para);
            }
            para->doubleParameters["SortingPolicy"] = rValue;
          } catch(...) {
            para = new dd4hep::rec::DoubleParameters;
            para->doubleParameters["SortingPolicy"] = rValue;
            ddsurf->detElement().addExtension<dd4hep::rec::DoubleParameters>(para);
          }
          printout(PrintLevel::DEBUG, LOG_SOURCE, "Added extension to %s, matching %s "
                   " path %s, type %s, zPos %3.5f, value %3.5f",
                   volumeName.c_str(), pathAndValues.first.c_str(),
                   ddsurf->detElement().path().c_str(),
                   ddsurf->detElement().type().c_str(),
                   zPosition, para->doubleParameters["SortingPolicy"]);
          break;
        }
      }// for all paths
    }

    return 1;
  }


} // namespace


DECLARE_APPLY(lcgeo_LinearSortingPolicy, ::addSortingPolicy)
