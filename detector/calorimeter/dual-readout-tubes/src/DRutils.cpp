#include <cmath>
#include "DRutils.h"

using namespace dd4hep;

namespace DDDRCaloTubes
{
    int fast_floor(double x)
    {
        return (int) x - (x < (int) x);
    }

    int fast_ceil(double x)
    {
        return (int) x + (x > (int) x);
    }

    bool check_for_integer(double x)
    {
        return (std::abs(x - std::round(x)) < 1e-6);
    }

    std::vector<double> get_plane_equation(const Position& point1, const Position& point2, const Position& point3)
    {
        Direction normal = (point2 - point1).Cross(point3 - point1).Unit();
        double A = normal.x();
        double B = normal.y();
        double C = normal.z();
        double D = -1.0 * (A * point1.x() + B * point1.y() + C * point1.z());
        std::vector<double> coefficients = {A, B, C, D};
        return coefficients;
    }

    Position get_intersection(const std::vector<double>& plane_coefficients, const Position& line_point, const Direction& line_direction)
    {
        double A = plane_coefficients[0];
        double B = plane_coefficients[1];
        double C = plane_coefficients[2];
        double D = plane_coefficients[3];
        double t = (-1.0*(A * line_point.x() + B * line_point.y() + C * line_point.z() + D)) / (A * line_direction.x() + B * line_direction.y() + C * line_direction.z());
        Position intersection = line_point + t * line_direction;
        return intersection;
    }

    Position get_intersection(const Direction& plane_normal, const Position& plane_point, const Position& line_point, const Direction& line_direction)
    {
        double t = (plane_normal.Dot(plane_point - line_point)) / plane_normal.Dot(line_direction);
        Position intersection = line_point + t * line_direction;
        return intersection;
    }

    double distance_from_plane(const std::vector<double>& plane_coefficients, const Position& point)
    {
        double A = plane_coefficients[0];
        double B = plane_coefficients[1];
        double C = plane_coefficients[2];
        double D = plane_coefficients[3];
        return std::abs(A * point.x() + B * point.y() + C * point.z() + D);
    }

} // namespace DDDRCaloTubes
