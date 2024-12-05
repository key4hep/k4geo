//**************************************************************************
// \file DREndcapTube.hh
// \brief: Auxiliary class to handle complex computations for the IDEA
//         dual-readout calorimeter with the capillary tubes technology
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 12 July 2024
//**************************************************************************

#ifndef DREndcapTube_H
#  define DREndcapTube_H

// This struct represents a plane corresponding to a tower's face.
// I construct it with the 4 edges (points) of the tower's face,
// however only 3 points will be used to define the plane.
struct Plane
{
    Vector3D P1, P2, P3, P4;
    Plane(Vector3D p1, Vector3D p2, Vector3D p3, Vector3D p4) : P1(p1), P2(p2), P3(p3), P4(p4) {};
};

// This struct represents a line towards negatize Z
// given its starting point (origin).
// The starting points will be distributed over the back face
// of the a tower. The corresponding lines will be used to find
// the intersection with the tower's faces.
struct ZLine
{
    Vector3D origin;
    Vector3D fuZ = Vector3D(0, 0, -1);
    ZLine(Vector3D P) : origin(P) {};
};

// Custom exception class for intersecting ZLines with Planes
class IntersectException : public std::exception
{
  public:
    const char* what() const noexcept override
    {
      return "IntersectLinePlane: Invalid intersection of ZLine with a Plane";
    }
};

class DREndcapTubeHelper
{
  public:
    // Constructor and de-constructor
    DREndcapTubeHelper() = default;
    ~DREndcapTubeHelper() = default;

  private:
    // Class fields
    //
    // Fields to be set
    double fTubeRadius;  // Radius of tubes (same for Scin and Cher)
    double fPhiZRot;  // tower phi dimension
    bool fRbool;  // set to true to build right endcap
    bool fCalBasicBool;  // has CalBasic already being called?
    double fInnerR;  // inner radius at theta = 0. rad
    double fTowerHeight;  // tower height/length
    double fNumZRot;  // number of rotations around Z axis
    double fDeltaTheta;  // theta dimension of current tower
    double fDeltaTheta2;  // theta dimension of next tower
    double fThetaOfCenter;  // theta angle pointing to center of current tower
    double fThetaOfCenter2;  // theta angle pointing to center of next tower
    // Fields to be computed
    double finnerR_new;
    double finnerR_new2;
    double fTrns_Length;
    Vector3D fTrns_Vector;
    Vector3D fV1;
    Vector3D fV2;
    Vector3D fV3;
    Vector3D fV4;
    double Ratio;
    double Ratio2;

  public:
    // Methods to set class fields
    void SetRbool(bool rBool) { fRbool = rBool; }  // set to True if building right endcap
    void SetInnerR(double innerR) { fInnerR = innerR; }
    void SetTowerHeight(double towerHeight) { fTowerHeight = towerHeight; }
    void SetNumZRot(int num)
    {
      fNumZRot = num;
      fPhiZRot = 2 * M_PI / (double)num;
    }
    void SetDeltaTheta(double theta) { fDeltaTheta = theta; }
    void SetThetaOfCenter(double theta) { fThetaOfCenter = theta; }
    void SetDeltaTheta2(double theta) { fDeltaTheta2 = theta; }
    void SetThetaOfCenter2(double theta) { fThetaOfCenter2 = theta; }
    void SetTubeRadius(double tubeRadius) { fTubeRadius = tubeRadius; }

    // Methods to calculate towers parameters
    //
    void CalBasic();
    // Get origin on tower within its phi slice/stave
    Vector3D GetOrigin(int i);
    // Get 8 vector/points for the tower Trap edges
    void Getpt(Vector3D* pt);

    // Methods to calculate tubes lengths
    //
    // This method finds the "intersection" between a ZLine (tube) and a Plane (tower face)
    void IntersectLinePlane(const ZLine& line, const Plane& plane, Vector3D& intersection);
    // This method calculates a tube length given the cylindrical structure of the tube
    // and each tower face
    double GetTubeLength(const Vector3D (&pt)[8], const Vector3D& point);
};

inline void DREndcapTubeHelper::CalBasic()
{
  fCalBasicBool = 1;

  // distance from center to front face of first tower
  // radius is 2500 mm which is same as distance from center to front face of very
  // last tower.
  // To get distance to first tower I divide this radius by cosine of first tower.
  // To understand impact of sin*tan part
  finnerR_new = fInnerR / (cos(fThetaOfCenter) - sin(fThetaOfCenter) * tan(fDeltaTheta / 2.));
  // distance from center to front face of second tower (same as above)
  finnerR_new2 = fInnerR / (cos(fThetaOfCenter2) - sin(fThetaOfCenter2) * tan(fDeltaTheta2 / 2.));

  // Half size at front face of first tower given by effective radius times tan of half delta
  double innerSide_half = finnerR_new * tan(fDeltaTheta / 2.);
  // Half size at back face
  // double outerSide_half = (finnerR_new+fTowerHeight)*tan(fDeltaTheta/2.);

  // Half size at front face of second tower (same as above)
  double innerSide_half2 = finnerR_new2 * tan(fDeltaTheta2 / 2.);

  // Distance from origin to center of first tower
  fTrns_Length = fTowerHeight / 2. + finnerR_new;

  // Vector from origin to center of first tower
  // Remember towers are placed starting from the one laying on the x-axis
  fTrns_Vector = dd4hep::rec::Vector3D(cos(fThetaOfCenter) * fTrns_Length, 0,
                                       sin(fThetaOfCenter) * fTrns_Length);

  double dx1 = fInnerR;  // inner radius hardcoded in detector construction
  double dxi = sin(fThetaOfCenter) * finnerR_new + innerSide_half * cos(fThetaOfCenter);
  double dxi2 = sin(fThetaOfCenter2) * finnerR_new2 + innerSide_half2 * cos(fThetaOfCenter2);
  Ratio = dxi / dx1;
  Ratio2 = dxi2 / dx1;

  // These vectors are distributed over the plane of the tower closest to the barrel
  // The difference between fV1 and fV2 (fV3 and fV4) is the correction of finnerR_new+fTowerHeight
  // that meake the vector start at the begin and end of the tower surface
  fV1 = dd4hep::rec::Vector3D((Ratio2)
                                * (cos(fThetaOfCenter) * finnerR_new
                                   + sin(fThetaOfCenter) * finnerR_new * tan(fDeltaTheta / 2.)),
                              0,
                              sin(fThetaOfCenter) * finnerR_new
                                - sin(fThetaOfCenter) * finnerR_new * tan(fDeltaTheta / 2.));

  fV2 = dd4hep::rec::Vector3D(
    (Ratio2)
      * (cos(fThetaOfCenter) * (finnerR_new + fTowerHeight)
         + sin(fThetaOfCenter) * (finnerR_new + fTowerHeight) * tan(fDeltaTheta / 2.)),
    0,
    sin(fThetaOfCenter) * (finnerR_new + fTowerHeight)
      - sin(fThetaOfCenter) * (finnerR_new + fTowerHeight) * tan(fDeltaTheta / 2.));

  fV3 = dd4hep::rec::Vector3D((Ratio)
                                * (cos(fThetaOfCenter) * finnerR_new
                                   - sin(fThetaOfCenter) * finnerR_new * tan(fDeltaTheta / 2.)),
                              0,
                              sin(fThetaOfCenter) * finnerR_new
                                + sin(fThetaOfCenter) * finnerR_new * tan(fDeltaTheta / 2.));

  fV4 = dd4hep::rec::Vector3D(
    (Ratio)
      * (cos(fThetaOfCenter) * (finnerR_new + fTowerHeight)
         - sin(fThetaOfCenter) * (finnerR_new + fTowerHeight) * tan(fDeltaTheta / 2.)),
    0,
    sin(fThetaOfCenter) * (finnerR_new + fTowerHeight)
      + sin(fThetaOfCenter) * (finnerR_new + fTowerHeight) * tan(fDeltaTheta / 2.));
}

// Get 8 vector/points for the tower Trap edges
inline void DREndcapTubeHelper::Getpt(dd4hep::rec::Vector3D* pt)
{
  double innerSide_half = finnerR_new * tan(fDeltaTheta / 2.);
  double outerSide_half = (finnerR_new + fTowerHeight) * tan(fDeltaTheta / 2.);

  if (fRbool == 1) {
    pt[0] =
      dd4hep::rec::Vector3D(-(fV1.x() * tan(fPhiZRot / 2.)), -innerSide_half, -fTowerHeight / 2.);
    pt[1] =
      dd4hep::rec::Vector3D((fV1.x() * tan(fPhiZRot / 2.)), -innerSide_half, -fTowerHeight / 2.);
    pt[2] =
      dd4hep::rec::Vector3D(-(fV3.x() * tan(fPhiZRot / 2.)), innerSide_half, -fTowerHeight / 2.);
    pt[3] =
      dd4hep::rec::Vector3D((fV3.x() * tan(fPhiZRot / 2.)), innerSide_half, -fTowerHeight / 2.);
    pt[4] =
      dd4hep::rec::Vector3D(-(fV2.x() * tan(fPhiZRot / 2.)), -outerSide_half, fTowerHeight / 2.);
    pt[5] =
      dd4hep::rec::Vector3D((fV2.x() * tan(fPhiZRot / 2.)), -outerSide_half, fTowerHeight / 2.);
    pt[6] =
      dd4hep::rec::Vector3D(-(fV4.x() * tan(fPhiZRot / 2.)), outerSide_half, fTowerHeight / 2.);
    pt[7] =
      dd4hep::rec::Vector3D((fV4.x() * tan(fPhiZRot / 2.)), outerSide_half, fTowerHeight / 2.);
  }
  else {
    pt[0] =
      dd4hep::rec::Vector3D(-(fV1.x() * tan(fPhiZRot / 2.)), -innerSide_half, -fTowerHeight / 2.);
    pt[1] =
      dd4hep::rec::Vector3D((fV1.x() * tan(fPhiZRot / 2.)), -innerSide_half, -fTowerHeight / 2.);
    pt[2] =
      dd4hep::rec::Vector3D(-(fV3.x() * tan(fPhiZRot / 2.)), innerSide_half, -fTowerHeight / 2.);
    pt[3] =
      dd4hep::rec::Vector3D((fV3.x() * tan(fPhiZRot / 2.)), innerSide_half, -fTowerHeight / 2.);
    pt[4] =
      dd4hep::rec::Vector3D(-(fV2.x() * tan(fPhiZRot / 2.)), -outerSide_half, fTowerHeight / 2.);
    pt[5] =
      dd4hep::rec::Vector3D((fV2.x() * tan(fPhiZRot / 2.)), -outerSide_half, fTowerHeight / 2.);
    pt[6] =
      dd4hep::rec::Vector3D(-(fV4.x() * tan(fPhiZRot / 2.)), outerSide_half, fTowerHeight / 2.);
    pt[7] =
      dd4hep::rec::Vector3D((fV4.x() * tan(fPhiZRot / 2.)), outerSide_half, fTowerHeight / 2.);
  }
}

// Get origin on tower within its phi slice/stave
inline dd4hep::rec::Vector3D DREndcapTubeHelper::GetOrigin(int i)
{
  if (fCalBasicBool == 0) {
    std::cout << "fCalBasicBool = 0" << std::endl;
    return dd4hep::rec::Vector3D();
  }
  else {
    if (fRbool == 1) {
      double x, y, z;
      x = cos(i * fPhiZRot) * fTrns_Vector.x();
      y = sin(i * fPhiZRot) * fTrns_Vector.x();
      z = fTrns_Vector.z();

      return dd4hep::rec::Vector3D(x, y, z);
    }
    else {
      double x, y, z;
      x = cos(i * fPhiZRot) * fTrns_Vector.x();
      y = sin(i * fPhiZRot) * fTrns_Vector.x();
      z = -fTrns_Vector.z();

      return dd4hep::rec::Vector3D(x, y, z);
    }
  }
}

// This method finds the "intersection" between a ZLine (tube) and a Plane (tower face).
// Only the first 3 points of the plane and needed to define the plane.
// The "intersection" is updated with the Vector3D at the intersection found.
// Exceptions are handled by this function if the line is parallel to the plane or the
// intersection is behind the line origin, i.e. below the backface of the tower.
// Both cases are to be considered errors.
inline void DREndcapTubeHelper::IntersectLinePlane(const ZLine& line, const Plane& plane,
                                                   Vector3D& intersection)
{
  Vector3D p0 = plane.P1;
  Vector3D p1 = plane.P2;
  Vector3D p2 = plane.P3;

  // Compute the plane normal
  Vector3D u = Vector3D(p1.x() - p0.x(), p1.y() - p0.y(), p1.z() - p0.z());
  Vector3D v = Vector3D(p2.x() - p0.x(), p2.y() - p0.y(), p2.z() - p0.z());
  Vector3D normal = Vector3D(u.y() * v.z() - u.z() * v.y(), u.z() * v.x() - u.x() * v.z(),
                             u.x() * v.y() - u.y() * v.x());

  double denominator =
    normal.x() * line.fuZ.x() + normal.y() * line.fuZ.y() + normal.z() * line.fuZ.z();

  if (std::fabs(denominator) < 1e-6) {
    throw IntersectException();  // The line is parallel to the plane
  }

  double d = normal.x() * p0.x() + normal.y() * p0.y() + normal.z() * p0.z();
  double t =
    (d - normal.x() * line.origin.x() - normal.y() * line.origin.y() - normal.z() * line.origin.z())
    / denominator;

  if (t < 0) {
    throw IntersectException();  // The intersection is behind the line origin
  }

  intersection = Vector3D(line.origin.x() + t * line.fuZ.x(), line.origin.y() + t * line.fuZ.y(),
                          line.origin.z() + t * line.fuZ.z());
};

// This method calculates a tube length given the cylindrical structure of the tube
// and each tower face.
// The intersection is calculated over 5 planes (tower front, up, down, left and right faces).
// For each plane the intersection is calculated over four points on the tube surface (up, down,
// left, right). The shortest intersection is the one to be used.
inline double DREndcapTubeHelper::GetTubeLength(const Vector3D (&pt)[8], const Vector3D& point)
{
  Plane FirstPlane(pt[0], pt[1], pt[2], pt[3]);  // front plane (smallest one)
  Plane SecondPlane(pt[1], pt[5], pt[3], pt[7]);  // right plane (looking at the front face)
  Plane ThirdPlane(pt[2], pt[3], pt[6], pt[7]);  // top plane
  Plane FourthPlane(pt[0], pt[4], pt[2], pt[6]);  // left plane
  Plane FifthPlane(pt[0], pt[1], pt[5], pt[4]);  // bottom plane
  std::array<Plane, 5> Planes{FirstPlane, SecondPlane, ThirdPlane, FourthPlane, FifthPlane};

  ZLine PointZLine(point);  // Axis of the tube with origin on tower back plane
  ZLine UpPointZLine(Vector3D(point.x(), point.y() + fTubeRadius, point.z()));
  ZLine DownPointZLine(Vector3D(point.x(), point.y() - fTubeRadius, point.z()));
  ZLine LeftPointZLine(Vector3D(point.x() - fTubeRadius, point.y(), point.z()));
  ZLine RightPointZLine(Vector3D(point.x() + fTubeRadius, point.y(), point.z()));
  std::array<ZLine, 4> Lines{UpPointZLine, DownPointZLine, LeftPointZLine, RightPointZLine};

  Vector3D intersection;  // intersection to be found
  double tubeLength{0.};  // tube length to be returned

  for (std::size_t k = 0; k < Lines.size(); k++) {  // loop over cylinder surface points
    for (std::size_t i = 0; i < Planes.size(); i++) {  // loop over tower's planes
      IntersectLinePlane(Lines.at(k), Planes.at(i), intersection);
      double fiberLengthPlane = (Lines.at(k).origin - intersection).r();
      if (k == 0 && i == 0)
        tubeLength = fiberLengthPlane;
      else if (tubeLength > fiberLengthPlane)
        tubeLength = fiberLengthPlane;
      else {
        ;
      }
    }  // end of loop over tower's planes
  }  // end of loop over cylinder surface points

  return tubeLength;
};

#endif  // DREndcapTube_H

//**************************************************************************
