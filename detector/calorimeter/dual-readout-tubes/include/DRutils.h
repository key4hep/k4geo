#ifndef DRutils_h
#define DRutils_h 1

#include <vector>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Objects.h"

// Utility functions for the construction of the barrel dual-readout calorimeter (based on tubes)
namespace DRBarrelTubes {
// Quick rounding functions
int fast_floor(double x);
int fast_ceil(double x);
// Make sure a given double is an integer (within a tolerance)
bool check_for_integer(double x);

// Functions which are used to calculate the tube lengths using the crossing point of a line and a plane
std::vector<double> get_plane_equation(const dd4hep::Position& point1, const dd4hep::Position& point2,
                                       const dd4hep::Position& point3);
dd4hep::Position get_intersection(const std::vector<double>& plane_coefficients, const dd4hep::Position& line_point,
                                  const dd4hep::Direction& line_direction);
dd4hep::Position get_intersection(const dd4hep::Direction& plane_normal, const dd4hep::Position& plane_point,
                                  const dd4hep::Position& line_point, const dd4hep::Direction& line_direction);
} // namespace DRBarrelTubes

#endif // DRutils_h
