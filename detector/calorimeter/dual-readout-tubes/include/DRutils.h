#ifndef DRutils_h
#define DRutils_h 1

#include <vector>

#include "DD4hep/Objects.h"
#include "DD4hep/DD4hepUnits.h"

using namespace dd4hep;

namespace DDDRCaloTubes
{
    int fast_floor(double x);
    int fast_ceil(double x);
    bool check_for_integer(double x);

    std::vector<double> get_plane_equation(const Position& point1, const Position& point2, const Position& point3);
    Position get_intersection(const std::vector<double>& plane_coefficients, const Position& line_point, const Direction& line_direction);
    Position get_intersection(const Direction& plane_normal, const Position& plane_point, const Position& line_point, const Direction& line_direction);
    double distance_from_plane(const std::vector<double>& plane_coefficients, const Position& point);
} // namespace DDDRCaloTubes

#endif // DRutils_h
