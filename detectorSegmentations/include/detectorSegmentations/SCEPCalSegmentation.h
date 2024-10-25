//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#ifndef SCEPCalSegmentation_h
#define SCEPCalSegmentation_h 1
#include "DDSegmentation/Segmentation.h"
#include "TVector3.h"
#include "DD4hep/DetFactoryHelper.h"
#include <vector>
#include <cmath>

namespace dd4hep {
namespace DDSegmentation {

class SCEPCalSegmentation : public Segmentation {
    public:
        SCEPCalSegmentation(const std::string& aCellEncoding);
        SCEPCalSegmentation(const BitFieldCoder* decoder);
        virtual ~SCEPCalSegmentation() override;

        virtual Vector3D position(const CellID& aCellID) const;

        virtual Vector3D myPosition(const CellID& aCellID) ;

        virtual CellID cellID(const Vector3D& aLocalPosition,
                              const Vector3D& aGlobalPosition,
                              const VolumeID& aVolumeID) const;

        VolumeID setVolumeID(int System, int Eta, int Phi, int Depth) const;
        CellID setCellID(int System, int Eta, int Phi, int Depth) const;

        int System(const CellID& aCellID) const;
        int Eta(const CellID& aCellID) const;
        int Phi(const CellID& aCellID) const;
        int Depth(const CellID& aCellID) const;

        int getFirst32bits(const CellID& aCellID) const { return (int)aCellID; }
        int getLast32bits(const CellID& aCellID) const;

        CellID convertFirst32to64(const int aId32) const { return (CellID)aId32; }
        CellID convertLast32to64(const int aId32) const;

        int System(const int& aId32) const { return System( convertFirst32to64(aId32) ); }
        int Eta(const int& aId32) const { return Eta( convertFirst32to64(aId32) ); }
        int Phi(const int& aId32) const { return Phi( convertFirst32to64(aId32) ); }
        int Depth(const int& aId32) const { return Depth( convertFirst32to64(aId32) ); }

        inline void setGeomParams(double Fdz, double Rdz,
                                  double nomfw, double nomth, 
                                  double EBz, double Rin,
                                  double sipmth,
                                  int phiSegments, int NProjectiveFill) {
            f_Fdz = Fdz;
            f_Rdz = Rdz;
            f_nomfw = nomfw;
            f_nomth = nomth;
            f_EBz = EBz;
            f_Rin = Rin;
            f_sipmth = sipmth;
            f_phiSegments = phiSegments;
            f_NProjectiveFill = NProjectiveFill;
        }

        inline double getFdz() const{ return f_Fdz; }
        inline double getRdz() const{ return f_Rdz; }
        inline double getnomfw() const{ return f_nomfw; }
        inline double getnomth() const{ return f_nomth; }
        inline double getEBz() const{ return f_EBz; }
        inline double getRin() const{ return f_Rin; }
        inline double getSipmth() const{ return f_sipmth; }
        inline int    getphiSegments() const{ return f_phiSegments; }
        inline int    getNProjectiveFill() const{ return f_NProjectiveFill; }

    protected:
        std::string fSystemId;
        std::string fEtaId;
        std::string fPhiId;
        std::string fDepthId;

        double f_Fdz;
        double f_Rdz;
        double f_nomfw;
        double f_nomth;
        double f_EBz;
        double f_Rin;
        double f_sipmth;
        int    f_phiSegments;
        int    f_NProjectiveFill;

    private:
        std::unordered_map<int, Vector3D> fPositionOf;

};
}
}

#endif
