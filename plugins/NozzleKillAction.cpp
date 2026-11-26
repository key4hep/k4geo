//==========================================================================
//  NozzleKillAction
//--------------------------------------------------------------------------
//  Geant4 stepping action that stops tracking particles deep inside the
//  forward nozzles of Muon Collider detectors.
//
//  Motivation
//  ----------
//  Muon Collider detector simulations of forward going particles hitting the
//  nozzles are expensive, and the part of the shower developing deep inside
//  the nozzle produces no useful information.
//
//  Kill region
//  -----------
//  Two cones with their tips on the beam axis at z = +KillZ and z = -KillZ,
//  each opening away from the interaction point with half-angle OpeningAngle
//  (measured from the beam axis), and ending at the end of the nozzle on
//  that side, zEnd.  A track whose post-step point satisfies
//
//      KillZ < |z| < zEnd   and   r < (|z| - KillZ) * tan(OpeningAngle)
//
//  is stopped, whatever it is and whatever material it is in (bore, tungsten,
//  borated polyethylene).  Nothing is killed beyond the end of the nozzle.
//
//  Nozzle geometry
//  ---------------
//  At the beginning of each run the action collects the polycone volumes of
//  the Geant4 region NozzleRegion (the MuColl compact files put every nozzle
//  piece in "NozzleRegion"), and builds the outer radius of the nozzle as a
//  function of |z| on each side.  zEnd is the largest |z| of those volumes.
//
//  Clearance check: the action computes the widest opening angle for which
//  the cones stay at least NozzleMargin inside the nozzle outer radius
//  everywhere between KillZ and zEnd.  If OpeningAngle is wider, it is
//  reduced to that limit and a warning is printed.  Set NozzleMargin < 0 to
//  skip the check (the cones are still capped at zEnd).
//
//  Configuration (from the ddsim steering file, using g4units):
//      SIM.action.step = ("NozzleKillAction", {
//          "KillZ":        200 * cm,        # cone tips at +-KillZ
//          "OpeningAngle": 5 * deg,         # cone half-angle from the beam axis
//          "NozzleMargin": 5 * cm,          # min. distance to the nozzle outer radius
//          "NozzleRegion": "NozzleRegion",  # region holding the nozzle volumes
//      })
//
//  Values are in Geant4 internal units (mm, rad).  The properties are read
//  and validated at the beginning of each run, after the steering file has
//  been applied.  The action only ever kills tracks; it never creates hits.
//==========================================================================

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "DD4hep/InstanceCount.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Volumes.h"
#include "DDG4/Geant4Context.h"
#include "DDG4/Geant4RunAction.h"
#include "DDG4/Geant4SteppingAction.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "G4ParticleDefinition.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"

#include "TGeoMatrix.h"
#include "TGeoPcon.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

class G4SteppingManager;

namespace dd4hep {
namespace sim {

  /// Stepping action stopping tracks deep inside the forward nozzles
  //
  // Units: dd4hep::cm (TGeo, = 1) and CLHEP::cm (Geant4, = 10) are both visible
  // in this namespace, so every unit below is written with its namespace.
  class NozzleKillAction : public Geant4SteppingAction {
  public:
    /// Standard constructor
    NozzleKillAction(Geant4Context* ctxt, const std::string& nam) : Geant4SteppingAction(ctxt, nam) {
      declareProperty("KillZ", m_killZ);
      declareProperty("OpeningAngle", m_openingAngle);
      declareProperty("NozzleMargin", m_nozzleMargin);
      declareProperty("NozzleRegion", m_nozzleRegion);
      // Property values are only final once the steering file has been
      // applied, so derived quantities are computed at begin of run.
      context()->runAction().callAtBegin(this, &NozzleKillAction::beginRun);
      InstanceCount::increment(this);
    }

    /// Destructor: report what the action did
    virtual ~NozzleKillAction() {
      if (m_nCalls)
        info("NozzleKill: %llu steps, killed gamma=%llu e+-=%llu mu+-=%llu neutron=%llu other=%llu", m_nCalls, m_nGamma,
             m_nElectron, m_nMuon, m_nNeutron, m_nOther);
      InstanceCount::decrement(this);
    }

    /// Begin-of-run callback: validate the properties and derive the cuts.
    /// Runs after the steering-file properties have been applied.
    void beginRun(const G4Run* /* run */) {
      if (m_killZ <= 0.0)
        except("NozzleKill: KillZ must be positive, got %.1f cm.", m_killZ / CLHEP::cm);
      if (m_openingAngle < 0.0 || m_openingAngle >= 90.0 * CLHEP::deg)
        except("NozzleKill: OpeningAngle must be in [0, 90) deg, got %.3f deg.", m_openingAngle / CLHEP::deg);

      const Sides sides = collectNozzleSegments();
      for (int side : {0, 1})
        if (m_zEnd[side] <= m_killZ)
          except("NozzleKill: KillZ = %.1f cm is beyond the %s end of the nozzle (|z| = %.1f cm).", m_killZ / CLHEP::cm,
                 kSideName[side], m_zEnd[side] / CLHEP::cm);

      double angle = m_openingAngle;
      if (m_nozzleMargin >= 0.0) {
        const double maxAngle = maxOpeningAngle(sides);
        if (angle > maxAngle) {
          warning("NozzleKill: OpeningAngle %.2f deg brings the cone within %.1f cm of the outer "
                  "radius of %s; limited to %.2f deg",
                  angle / CLHEP::deg, m_nozzleMargin / CLHEP::cm, m_nozzleRegion.c_str(), maxAngle / CLHEP::deg);
          angle = maxAngle;
        }
      }
      m_tanAngle = std::tan(angle);

      info("NozzleKill: active - killing tracks in cones from |z| = %.1f cm "
           "with opening angle %.2f deg",
           m_killZ / CLHEP::cm, angle / CLHEP::deg);
      info("NozzleKill: cones end with the nozzle at z = %.1f cm and z = -%.1f cm", m_zEnd[0] / CLHEP::cm,
           m_zEnd[1] / CLHEP::cm);
    }

    /// Stepping callback
    void operator()(const G4Step* step, G4SteppingManager* /* mgr */) override {
      ++m_nCalls;
      const G4ThreeVector& pos = step->GetPostStepPoint()->GetPosition();

      // Cheap test first: most steps are nowhere near the nozzle.
      const double absZ = std::fabs(pos.z());
      const double dz = absZ - m_killZ;
      if (dz <= 0.0 || absZ >= m_zEnd[sideOf(pos.z())])
        return;

      // Inside a cone: kill whatever it is.
      const double rCone = dz * m_tanAngle;
      if (pos.perp2() >= rCone * rCone)
        return;

      G4Track* track = step->GetTrack();
      track->SetTrackStatus(fStopAndKill);
      const int pdg = std::abs(track->GetParticleDefinition()->GetPDGEncoding());
      ++(pdg == 22 ? m_nGamma : pdg == 11 ? m_nElectron : pdg == 13 ? m_nMuon : pdg == 2112 ? m_nNeutron : m_nOther);
    }

  private:
    /// Straight piece of a nozzle outer surface, in |z| (Geant4 units, z1 < z2)
    struct Segment {
      double z1, z2, r1, r2;
    };
    /// Outer-surface segments of the nozzle volumes on each side (0: z > 0, 1: z < 0)
    using Sides = std::array<std::vector<Segment>, 2>;

    static constexpr const char* kSideName[2] = {"+z", "-z"};
    static int sideOf(double z) { return z > 0.0 ? 0 : 1; }

    /// Outer radius of the nozzle at |z| on one side: the largest rmax of the pieces covering z, -1 in a gap
    static double outerRadius(const std::vector<Segment>& segs, double z) {
      double rOut = -1.0;
      for (const Segment& s : segs)
        if (z >= s.z1 && z <= s.z2)
          rOut = std::max(rOut, s.r1 + (s.r2 - s.r1) * (z - s.z1) / (s.z2 - s.z1));
      return rOut;
    }

    /// Collect the outer-surface segments of the nozzle volumes and set the nozzle end on each side
    Sides collectNozzleSegments() {
      const double toG4 = CLHEP::cm / dd4hep::cm; // TGeo length -> Geant4 length
      Detector& description = context()->detectorDescription();
      Sides sides;
      m_zEnd = {0.0, 0.0};
      int nVolumes = 0;
      for (const auto& entry : description.detectors()) {
        DetElement de(entry.second);
        if (!de.isValid() || !de.placement().isValid())
          continue;
        Volume vol = de.volume();
        Region reg = vol.region();
        if (!reg.isValid() || reg.name() != m_nozzleRegion)
          continue;
        const TGeoPcon* pcon = dynamic_cast<const TGeoPcon*>(vol.solid().ptr());
        if (!pcon) {
          warning("NozzleKill: %s is in %s but is not a polycone; ignored.", de.name(), m_nozzleRegion.c_str());
          continue;
        }

        // World position of the local beam axis: the piece must be coaxial with z.
        const Alignment nominal = de.nominal(); // named, so the matrix reference below outlives no temporary
        const TGeoHMatrix& m = nominal.worldTransformation();
        const double* rot = m.GetRotationMatrix();
        const double* tr = m.GetTranslation();
        if (std::hypot(rot[2], rot[5]) > 1e-9 || std::hypot(tr[0], tr[1]) > 1e-6 * dd4hep::mm)
          except("NozzleKill: %s is not coaxial with the beam axis; NozzleKillAction needs it to be.", de.name());

        ++nVolumes;
        for (int i = 0; i + 1 < pcon->GetNz(); ++i) {
          double za = (tr[2] + rot[8] * pcon->GetZ(i)) * toG4, ra = pcon->GetRmax(i) * toG4;
          double zb = (tr[2] + rot[8] * pcon->GetZ(i + 1)) * toG4, rb = pcon->GetRmax(i + 1) * toG4;
          if (za * zb < 0.0)
            continue; // crosses the interaction point, not a forward nozzle piece
          if (std::fabs(za) > std::fabs(zb)) {
            std::swap(za, zb);
            std::swap(ra, rb);
          }
          if (std::fabs(zb) <= std::fabs(za))
            continue; // zero length
          const int side = sideOf(za + zb);
          sides[side].push_back({std::fabs(za), std::fabs(zb), ra, rb});
          m_zEnd[side] = std::max(m_zEnd[side], std::fabs(zb));
        }
      }
      if (!nVolumes)
        except("NozzleKill: no polycone volumes found in region %s. Set NozzleRegion to the region "
               "holding the nozzle volumes.",
               m_nozzleRegion.c_str());
      return sides;
    }

    /// Widest opening angle keeping both cones NozzleMargin inside the nozzle outer radius
    double maxOpeningAngle(const Sides& sides) const {
      const double step = 1.0 * CLHEP::mm;
      double maxTan = std::numeric_limits<double>::max();
      double limitZ = 0.0;
      for (int side : {0, 1}) {
        const std::vector<Segment>& segs = sides[side];
        const double zEnd = m_zEnd[side];

        // Sample the outer radius on a fine grid plus every segment end point.
        std::vector<double> zs;
        for (double z = m_killZ + step; z < zEnd; z += step)
          zs.push_back(z);
        for (const Segment& s : segs)
          for (double z : {s.z1, s.z2})
            if (z > m_killZ)
              zs.push_back(z);

        for (double z : zs) {
          const double rOut = outerRadius(segs, z);
          if (rOut < 0.0)
            continue; // gap between nozzle pieces
          const double room = rOut - m_nozzleMargin;
          if (room <= 0.0)
            except("NozzleKill: at |z| = %.1f cm the nozzle is only %.1f cm in outer radius, less than "
                   "NozzleMargin = %.1f cm. This usually means the nozzle has a gap there (e.g. a cut for "
                   "a forward detector, bridged only by a thin beam tube). The cones end at the end of the "
                   "nozzle (|z| = %.1f cm), so they would run through the gap and kill tracks in whatever "
                   "sits in it; NozzleKillAction cannot be used with this nozzle geometry.",
                   z / CLHEP::cm, rOut / CLHEP::cm, m_nozzleMargin / CLHEP::cm, zEnd / CLHEP::cm);
          const double t = room / (z - m_killZ);
          if (t < maxTan) {
            maxTan = t;
            limitZ = side == 0 ? z : -z;
          }
        }
      }
      const double maxAngle = std::atan(maxTan);
      info("NozzleKill: widest opening angle keeping %.1f cm clearance from %s is %.2f deg (limited at z = %.1f cm)",
           m_nozzleMargin / CLHEP::cm, m_nozzleRegion.c_str(), maxAngle / CLHEP::deg, limitZ / CLHEP::cm);
      return maxAngle;
    }

    /// Position of the cone tips on the beam axis (+-KillZ)
    double m_killZ = 200.0 * CLHEP::cm;
    /// Cone half-angle, measured from the beam axis
    double m_openingAngle = 5.0 * CLHEP::deg;
    /// Minimum distance between the cone and the nozzle outer radius; < 0 skips the check
    double m_nozzleMargin = 5.0 * CLHEP::cm;
    /// Region holding the nozzle volumes
    std::string m_nozzleRegion = "NozzleRegion";

    /// Derived at begin of run; with zEnd = 0 the action kills nothing until then
    double m_tanAngle = 0.0;
    /// End of the nozzle in |z| on each side (0: z > 0, 1: z < 0)
    std::array<double, 2> m_zEnd = {0.0, 0.0};

    unsigned long long m_nCalls = 0;
    unsigned long long m_nGamma = 0;
    unsigned long long m_nElectron = 0;
    unsigned long long m_nMuon = 0;
    unsigned long long m_nNeutron = 0;
    unsigned long long m_nOther = 0;
  };

} // namespace sim
} // namespace dd4hep

#include "DDG4/Factories.h"
DECLARE_GEANT4ACTION(NozzleKillAction)
