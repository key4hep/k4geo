#ifndef DRutils_h
#define DRutils_h 1

#include <vector>

#include "DD4hep/Objects.h"
#include "DD4hep/DD4hepUnits.h"

using namespace dd4hep;


// Utility functions for the construction of the barrel dual-readout calorimeter (based on tubes)
namespace DRBarrelTubes
{
    // Quick rounding functions
    int fast_floor(double x);
    int fast_ceil(double x);
    // Make sure a given double is an integer (within a tolerance)
    bool check_for_integer(double x);

    // Functions which are used to calculate the tube lengths using the crossing point of a line and a plane
    std::vector<double> get_plane_equation(const Position& point1, const Position& point2, const Position& point3);
    Position get_intersection(const std::vector<double>& plane_coefficients, const Position& line_point, const Direction& line_direction);
    Position get_intersection(const Direction& plane_normal, const Position& plane_point, const Position& line_point, const Direction& line_direction);
} // namespace DRBarrelTubes

#endif // DRutils_h
