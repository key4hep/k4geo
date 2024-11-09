#ifndef DETECTORSEGMENTATIONS_HCALPHITHETA_K4GEO_H
#define DETECTORSEGMENTATIONS_HCALPHITHETA_K4GEO_H

// FCCSW
#include "detectorSegmentations/GridTheta_k4geo.h"

/** FCCSWHCalPhiTheta_k4geo Detector/detectorSegmentations/detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h FCCSWHCalPhiTheta_k4geo.h
 *
 *  Segmentation in theta and phi.
 *  Based on GridTheta_k4geo, addition of azimuthal angle coordinate.
 *
 *  Rectangular shape cells are defined in r-z plan and all the hits within a defined cell boundaries are assigned the same cellID.
 *
 */

namespace dd4hep {
  namespace DDSegmentation {
    class FCCSWHCalPhiTheta_k4geo : public GridTheta_k4geo {
      public:
        /// default constructor using an arbitrary type
        FCCSWHCalPhiTheta_k4geo(const std::string& aCellEncoding);
        /// Default constructor used by derived classes passing an existing decoder
        FCCSWHCalPhiTheta_k4geo(const BitFieldCoder* decoder);

        /// destructor
        virtual ~FCCSWHCalPhiTheta_k4geo() = default;

        /**  Get the postion of the geometric center of the cell based on the cellID
        *   @param[in] aCellID
        *   return the global coordinates of cell center
        */
        virtual Vector3D position(const CellID& aCellID) const;

        /**  Assign a cellID based on the global position.
        *  @param[in] aLocalPosition (not used).
        *  @param[in] aGlobalPosition position in the global coordinates.
        *  @param[in] aVolumeId ID of a volume.
        *  return Cell ID.
        */
        virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                        const VolumeID& aVolumeID) const;

        /**  Find neighbours of the cell
        *   @param[in] aCellId ID of a cell.
        *   return vector of neighbour cellIDs.
        */
	std::vector<uint64_t> neighbours(const CellID& cID) const;

        /**  Calculate layer radii and edges in z-axis, then define cell edges in each layer using defineCellEdges().
        */
        void defineCellsInRZplan() const;

        /**  Define cell edges in z-axis for the given layer.
        *   Logic:
        *      1) Find theta bin centers that fit within the given layer; 
        *      2) Define a cell edge in z-axis as the middle of each pair of theta bin centers
        *   @param[in] layer index 
        */
	void defineCellEdges(const unsigned int layer) const;

        /**  Determine the azimuthal angle of HCal cell based on the cellID.
        *   @param[in] aCellId ID of a cell.
        *   return Phi.
        */
        double phi(const CellID& aCellID) const;

        /**  Get the grid size in phi.
        *   return Grid size in phi.
        */
        inline double gridSizePhi() const { return 2 * M_PI / static_cast<double>(m_phiBins); }

        /**  Get the number of bins in azimuthal angle.
        *   return Number of bins in phi.
        */
        inline int phiBins() const { return m_phiBins; }

        /**  Get the coordinate offset in azimuthal angle.
        *   return The offset in phi.
        */
        inline double offsetPhi() const { return m_offsetPhi; }

        /**  Get the vector of theta bins (cells) in a give layer.
        */
	inline std::vector<int> thetaBins(const uint layer) const {
          if(m_radii.empty()) defineCellsInRZplan();
          if(!m_thetaBins.empty()) return m_thetaBins[layer];
          else return std::vector<int>();
        }

        /**  Get the coordinate offset in z-axis.
        *   return The offset in z.
        */
        inline std::vector<double> offsetZ() const { return m_offsetZ; }

        /**  Get the z width of the layer.
        *   return the z width.
        */
	inline std::vector<double> widthZ() const { return m_widthZ; }

        /**  Get the coordinate offset in radius.
        *   return the offset in radius.
        */
	inline std::vector<double> offsetR() const { return m_offsetR; }

        /**  Get the number of layers.
        *   return the number of layers.
        */
	inline std::vector<int> numLayers() const { return m_numLayers; }

        /**  Get the dR of layers.
        *   return the dR.
        */
	inline std::vector<double> dRlayer() const { return m_dRlayer; }

        /**  Get the field name for azimuthal angle.
        *   return The field name for phi.
        */
        inline const std::string& fieldNamePhi() const { return m_phiID; }

        /**  Determine the minimum and maximum polar angle of HCal cell based on the cellID.
        *   @param[in] aCellId ID of a cell.
        *   return Theta.
        */
        std::array<double, 2> cellTheta(const CellID& cID) const;

        /**  Get the min and max layer indexes of each HCal part.
        * For Endcap, returns the three elements vector, while for Barrel - single element vector.
        */
	std::vector<std::pair<uint,uint> > getMinMaxLayerId() const ;

        /**  Set the number of bins in azimuthal angle.
        *   @param[in] aNumberBins Number of bins in phi.
        */
        inline void setPhiBins(int bins) { m_phiBins = bins; }

        /**  Set the coordinate offset in azimuthal angle.
        *   @param[in] aOffset Offset in phi.
        */
        inline void setOffsetPhi(double offset) { m_offsetPhi = offset; }

        /**  Set the coordinate offset in z-axis.
        *   @param[in] offset in z (centre of the layer).
        */
        inline void setOffsetZ(std::vector<double> const&offset) { m_offsetZ = offset; }

        /**  Set the z width.
        *   @param[in] width in z.
        */
        inline void setWidthZ(std::vector<double> const&width) { m_widthZ = width; }

        /**  Set the coordinate offset in radius.
        *   @param[in] offset in radius.
        */
        inline void setOffsetR(std::vector<double> const&offset) { m_offsetR = offset; }

        /**  Set the number of layers.
        *   @param[in] number of layers
        */
        inline void setNumLayers(std::vector<int> const&num) { m_numLayers = num; }

        /**  Set the dR of layers.
        *   @param[in] dR of layers.
        */
        inline void setdRlayer(std::vector<double> const&dRlayer) { m_dRlayer = dRlayer; }


        /**  Set the field name used for azimuthal angle.
        *   @param[in] aFieldName Field name for phi.
        */
        inline void setFieldNamePhi(const std::string& fieldName) { m_phiID = fieldName; }

      protected:
        /// determine the azimuthal angle phi based on the current cell ID
        double phi() const;
        /// the number of bins in phi
        int m_phiBins;
        /// the coordinate offset in phi
        double m_offsetPhi;
        /// the field name used for phi
        std::string m_phiID;
        /// the z offset of middle of the layer
        std::vector<double> m_offsetZ;
        /// the z width of the layer
        std::vector<double> m_widthZ;
        /// the radius of the layer
        std::vector<double> m_offsetR;
        /// number of layers
        std::vector<int> m_numLayers;
        /// dR of the layer
        std::vector<double> m_dRlayer;
        /// radius of each layer
        mutable std::vector<double> m_radii;
        /// z-min and z-max of each layer
        mutable std::vector<std::pair<double, double> >  m_layerEdges;
        /// dR of each layer
        mutable std::vector<double>  m_layerDepth;
        /// theta bins (cells) in each layer
        mutable std::vector<std::vector<int> > m_thetaBins;
        /// z-min and z-max of each cell (theta bin) in each layer
        mutable std::vector<std::unordered_map<int,std::pair<double, double> > > m_cellEdges;
    };
  }
}
#endif /* DETSEGMENTATION_HCALPHITHETA_H */
