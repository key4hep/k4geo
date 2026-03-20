#ifndef WIRE_TRACKER_INFO_H
#define WIRE_TRACKER_INFO_H

#include "DDRec/DetectorData.h"
#include <Evaluator/DD4hepUnits.h>
#include <map>
#include <Math/Vector3D.h>
#include <Math/GenVector/RotationZ.h>
#include <vector>

namespace dd4hep {  namespace rec {


/* Data structure to store Wire Trackers geometry parameters
 *
 * Global parameters are members of this class
 *
 * Parameters of each layer are stored in the helper
 * class LayerInfo.
 *
 * The member database is a container which stores one
 * LayerInfo object per layer.
 *
 * To use this class, instantiate an object and define the global parameters.
 * Then call the method BuildLayerDatabase, which calculates
 * the parameters of each layer based on the global parameters.
 */
struct WireTracker_info_struct
{
public:

    //--------------------------------------------------------------
    // Use alias of types to show more clearly what the variable is
    // if everything is double, the code is less readable

    /// alias for root 3D vectors
    using Vector3D = ROOT::Math::XYZVector;
    /// type for layer number
    using layer_t = int;
    /// type for lengths
    using length_t = double;
    /// type for angles
    using angle_t = double;

    //--------------------------------------------------------------
    // Per-layer geometry record

    // Internal helper struct for defining the layer layout
    struct LayerInfo
    {
        /// layer number
        layer_t layer = {0};
        /// 2x number of cells in that layer
        int nwires = {0};
        /// cell parameter
        double height_z0 = {0};
        /// cell parameter
        double width_z0 = {0};
        /// stereoangle at z=0
        angle_t stereo_sw_z0 = {0};
        /// radius (cylindrical coord) of sensitive wire
        length_t radius_sw_z0 = {0};
        /// radius (cylindrical coord) of 'down' field wires
        length_t radius_fdw_z0 = {0};
        /// radius (cylindrical coord) of 'up' field wires
        length_t radius_fuw_z0 = {0};
    };

    //--------------------------------------------------------------
    // Global geometry parameters

    /// Half length of the active volume
    length_t Lhalf = {0};
    /// Inner radius of the active volume
    length_t rin = {0};
    /// Outer radius of the active volume
    length_t rout = {0};

    /// input number of layers in each superlayer
    layer_t nlayersPerSuperlayer = {0};
    /// input number of superlayers
    /// superlayer is an abstract level of grouping layers used to
    /// parametrize the increment of cells in each layer
    layer_t nsuperlayers = {0};
    /// Calculated as dch_nlayersPerSuperlayer * nsuperlayers
    layer_t nlayers = {0};

    /// global twist angle
    /// alternating layers will change its sign
    angle_t twist_angle = {0};

    // ------------------------------------------------------------------
    //  Convenience setters

    void Set_lhalf(length_t _Lhalf){Lhalf= _Lhalf;}
    void Set_rin  (length_t _rin  ){rin  = _rin;  }
    void Set_rout (length_t _rout ){rout = _rout; }

    void Set_nlayersPerSuperlayer(int _nlayersPerSuperlayer){nlayersPerSuperlayer = _nlayersPerSuperlayer;}
    void Set_nsuperlayers        (int _nsuperlayers        ){nsuperlayers         = _nsuperlayers;        }

    void Set_twist_angle (length_t _dch_twist_angle ){twist_angle = _dch_twist_angle;}

    /// ------------------------------------------------------------------
    // Wire-distance calculations  (used by digitisation / reconstruction)
    // Notation: `ilayer` is the correlative number of the layer, `layer` is reserved to number within a superlayer
    layer_t  CalculateILayerFromCellIDFields(int layer, int superlayer) const { layer_t ilayer = layer + (this->nlayersPerSuperlayer)*superlayer + 1; return ilayer;}
    Vector3D Calculate_hitpos_to_wire_vector(int superlayer, int ilayer, int sector, int nphi, const Vector3D& hit_position /*in cm*/) const;
    Vector3D Calculate_wire_vector_ez       (int superlayer, int ilayer, int sector, int nphi) const;
    Vector3D Calculate_wire_z0_point        (int superlayer, int ilayer, int sector, int nphi) const;

    //--- virtual interface methods
    virtual angle_t Get_cell_phi_angle      (int superlayer, int ilayer, int sector, int nphi) const = 0;

    // ------------------------------------------------------------------
    //  Per-layer database

    /// map to store parameters for each layer
    std::map<layer_t, LayerInfo> database;

    bool IsDatabaseEmpty() const { return database.empty(); }
    /// Check if outer volume is not zero (0 < Lhalf*rout), and if the database was filled
    bool IsValid() const {return ((0 < Lhalf*rout) && (not IsDatabaseEmpty() ) );  }

    //--- virtual interface methods
    virtual void BuildLayerDatabase() = 0;
    virtual void ShowDatabase(std::ostream& io) const = 0;

    virtual ~WireTracker_info_struct(){};
};

typedef StructExtension<WireTracker_info_struct> WireTracker_info ;
inline std::ostream& operator<<( std::ostream& io , const WireTracker_info& d ){d.ShowDatabase(io); return io;}

//--------------------------------------------------------------
// inline implementations of wire-distance calculations 
//

///
/// Calculate the postion vector of the wire centre (i.e. at z=0)
///
inline WireTracker_info_struct::Vector3D
WireTracker_info_struct::Calculate_wire_z0_point(int superlayer, int ilayer, int sector, int nphi) const {
  auto&    l   = this->database.at(ilayer);
  double   rz0 = l.radius_sw_z0;
  Vector3D p1(rz0, 0, 0);
  double   phi_z0 = Get_cell_phi_angle(superlayer, ilayer, sector, nphi);
  ROOT::Math::RotationZ z_rotation = ROOT::Math::RotationZ(phi_z0);
  z_rotation(p1);
  return p1;
}

///
/// Calculate wire direction vector
///
inline WireTracker_info_struct::Vector3D
WireTracker_info_struct::Calculate_wire_vector_ez(int superlayer, int ilayer, int sector, int nphi) const {
  auto& l = this->database.at(ilayer);

  // See original paper Hoshina et al, Computer Physics Communications 153 (2003) 3
  // eq. 2.9, for the definition of ez, vector along the wire

  // initialize some variables
  double rz0        = l.radius_sw_z0;
  double stereo     = l.stereo_sw_z0;
  // kappa is the same as in eq. 2.9 
  // using the relation: tan(stereoangle) = R(z=0)   / (L/2) * tan( twist_angle/2))
  double kappa = tan(stereo) / rz0;

  //--- calculating wire position
  // the points p1 and p2 correspond to the ends of the wire

  // point 1
  double dy = rz0 * kappa * Lhalf;  // m

  // (x, y, z)  [m]
  Vector3D p1(rz0, -dy, -Lhalf); // point 1
  Vector3D p2(rz0, +dy, +Lhalf); // point 2

  // calculate phi rotation of whole twisted tube, ie, rotation at z=0
  double phi_z0 = Get_cell_phi_angle(superlayer, ilayer, sector, nphi);
  auto z_rotation = ROOT::Math::RotationZ(phi_z0);
  z_rotation(p1);
  z_rotation(p2);

  //--- end calculating wire position

  return (p2 -p1).Unit();
}

///
/// Calculate vector from hit position to wire
///
inline WireTracker_info_struct::Vector3D
WireTracker_info_struct::Calculate_hitpos_to_wire_vector(int isuperlayer, int ilayer, int isector, int nphi, const Vector3D& hit_position /*in cm*/) const {
  // Solution distance from a point to a line given here:
  // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
  Vector3D n = this->Calculate_wire_vector_ez(isuperlayer, ilayer, isector, nphi);
  Vector3D a = this->Calculate_wire_z0_point(isuperlayer, ilayer, isector, nphi);
  // Remember using cm as natural units of DD4hep consistently!
  // Vector3D p {hit_position.x()*MM_TO_CM,hit_position.y()*MM_TO_CM,hit_position.z()*MM_TO_CM};

  Vector3D a_minus_p       = a - hit_position;
  double   a_minus_p_dot_n = a_minus_p.Dot(n);
  Vector3D scaled_n        = a_minus_p_dot_n * n;
  return (a_minus_p - scaled_n);
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//                                                              //
//                     Derived Structures                       //
//                                                              //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//===========================================================
// Drift chamber info
//===========================================================

/* Data structure to store Drift Chamber geometry parameters,
 * implementing virtual methods of WireTracker_info_struct.
 */
struct DCH_info_struct : WireTracker_info_struct {

    //--------------------------------------------------------------
    //  DCH-specific global parameters

    /// Inner and outer guard wires radius
    length_t guard_inner_r_at_z0 = {0};
    length_t guard_outer_r_at_zL2 = {0};

    /// number of cells of first layer
    int ncell0 = {0};
    /// increment the number of cells for each superlayer as:
    ///   ncells(ilayer) = dch_ncell0 + increment*superlayer(ilayer)
    ///   See WireTracker_info::Get_nsuperlayer_minus_1(ilayer)
    int ncell_increment = {0};
    /// cells within the same layer may be grouped into sectors, not in use atm
    int ncell_per_sector = {0};
    /// Cell width for the first layer
    double first_width = {0};
    /// Cell radius for the first layer
    length_t first_sense_r = {0};

    //--------------------------------------------------------------
    // Setters

    void Set_guard_rin_at_z0  (length_t _rin_z0_guard  ){guard_inner_r_at_z0  = _rin_z0_guard;   }
    void Set_guard_rout_at_zL2(length_t _rout_zL2_guard){guard_outer_r_at_zL2 = _rout_zL2_guard; }

    void Set_ncell0          (int _ncell0          ){ncell0           = _ncell0;          }
    void Set_ncell_increment (int _ncell_increment ){ncell_increment  = _ncell_increment; }
    void Set_ncell_per_sector(int _ncell_per_sector){ncell_per_sector = _ncell_per_sector;}

    void Set_first_width  (double _first_width  ){first_width   = _first_width;   }
    void Set_first_sense_r(double _first_sense_r){first_sense_r = _first_sense_r; }


    //--------------------------------------------------------------
    // DCH geo helpers

    /// calculate superlayer for a given ilayer.
    /// WARNING: division of integers on purpose!
    int Get_nsuperlayer_minus_1(int ilayer){ return int((ilayer-1)/nlayersPerSuperlayer);}

    /// Calculate radius at z=L/2 given at z=0
    length_t Radius_zLhalf(length_t r_z0) const {
        return r_z0/cos(twist_angle/2/dd4hep::rad);
    }

    /// tan(stereoangle) = R(z=0)   / (L/2) * tan( twist_angle/2)
    angle_t stereoangle_z0(length_t r_z0) const {
        return atan( r_z0/Lhalf*tan(twist_angle/2/dd4hep::rad));
    }

    /// tan(stereoangle) = R(z=L/2) / (L/2) * sin( twist_angle/2)
    angle_t stereoangle_zLhalf(length_t r_zLhalf) const {
        return atan( r_zLhalf/Lhalf*sin(twist_angle/2/dd4hep::rad));
    }

    /// separation between wires (along the circle)
    length_t Pitch_z0(LayerInfo layer, length_t r_z0) const {
        return dd4hep::twopi*r_z0/layer.nwires;
    };

    /// WireLength = 2*dch_Lhalf/cos(atan(Pitch_z0(r_z0)/(2*dch_Lhalf)))/cos(stereoangle_z0(r_z0))
    length_t WireLength(int nlayer, length_t r_z0) const {
        auto Pitch_z0 = this->Pitch_z0(database.at(nlayer), r_z0);
        return  2*this->Lhalf/cos(atan(Pitch_z0/(2*this->Lhalf)))/cos(stereoangle_z0(r_z0)/dd4hep::rad) ;
    };

    ///  stereo angle is positive for odd layer number
    bool IsStereoPositive(LayerInfo layer_info) const {
        return (1 == layer_info.layer%2);
    };

    /// calculate sign based on IsStereoPositive
    int StereoSign(LayerInfo layer) const {
        return (IsStereoPositive(layer)*2 - 1);
    };

    /// number of cells per layer, calculated as the total number of wires / 2
    int Get_ncells(int ilayer) const {
        return database.at(ilayer).nwires/2;
    };

    /// dphi between two sensing wires
    angle_t Get_phi_width(int ilayer) const {
        return (dd4hep::twopi/Get_ncells(ilayer))*dd4hep::rad;
    };

    //--------------------------------------------------------------
    // DCH geo interface implementation

    /// phi positioning, adding offset for odd ilayers
    /// there is a staggering in phi for al {return database.at(ilayer).nwires/2;}ternating layers, 0.25*cell_phi_width*(ilayer%2);
    angle_t Get_cell_phi_angle(int, int ilayer, int, int iphi) const override final { 
        return (Get_phi_width(ilayer) * (iphi + 0.25*(ilayer%2)));
    }

    //--------------------------------------------------------------
    // Database construction

    // DCH database builder
    void BuildLayerDatabase() override final {
        // do not fill twice the database
        if( not this->IsDatabaseEmpty() ) return;

        auto ff_check_positive_parameter = [](double p, std::string pname) -> void {
            if(p<=0)throw std::runtime_error(Form("DCH: %s must be positive",pname.c_str()));
            return;
        };
        ff_check_positive_parameter(this->rin  ,"inner radius");
        ff_check_positive_parameter(this->rout ,"outer radius");
        ff_check_positive_parameter(this->Lhalf,"half length" );

        ff_check_positive_parameter(this->guard_inner_r_at_z0 ,"inner radius of guard wires" );
        ff_check_positive_parameter(this->guard_outer_r_at_zL2,"outer radius of guard wires" );

        ff_check_positive_parameter(this->ncell0,"ncells in the first layer" );
        ff_check_positive_parameter(this->ncell_increment,"ncells increment per superlayer" );
        // ff_check_positive_parameter(this->ncell_per_sector,"ncells per sector" );

        // if dch_ncell_per_sector is not divisor of dch_ncell0 and dch_ncell_increment
        // throw an error
        if( 0 != (ncell0 % ncell_per_sector) || 0 != (ncell_increment % ncell_per_sector) )
            throw std::runtime_error("dch_ncell_per_sector is not divisor of dch_ncell0 or dch_ncell_increment");

        ff_check_positive_parameter(this->nsuperlayers,"number of superlayers" );
        ff_check_positive_parameter(this->nlayersPerSuperlayer,"number of layers per superlayer" );

        /// nlayers = nsuperlayers * nlayersPerSuperlayer
        /// default: 112 = 14 * 8
        this->nlayers = this->nsuperlayers * this->nlayersPerSuperlayer;

        ff_check_positive_parameter(this->first_width,"width of first layer cells" );
        ff_check_positive_parameter(this->first_sense_r,"radius of first layer cells" );

        // initialize layer 1 from input parameters
        {
            LayerInfo layer1_info;
            layer1_info.layer         = 1;
            layer1_info.nwires        = 2*this->ncell0;
            layer1_info.height_z0     = first_width;    // cell height = cell width
            layer1_info.radius_sw_z0  = first_sense_r;
            layer1_info.radius_fdw_z0 = first_sense_r - 0.5*first_width;
            layer1_info.radius_fuw_z0 = first_sense_r + 0.5*first_width;
            layer1_info.width_z0      = dd4hep::twopi*first_sense_r/this->ncell0;
            layer1_info.stereo_sw_z0  = StereoSign(layer1_info) * this->stereoangle_z0(first_sense_r);

            this->database.emplace(layer1_info.layer, layer1_info);
        }

        // some parameters of the following layer are calculated based on the previous ones
        // the rest are left as methods of WireTracker_info or WireTracker_info_layer class
        // loop over all layers, skipping the first one
        for(int ilayer = 2; ilayer<= this->nlayers; ++ilayer) {
            // initialize empty object, parameters are set later
            LayerInfo layer_info;

            // the loop counter actually corresponds to the layer number
            layer_info.layer = ilayer;
            // nwires is twice the number of cells in this particular layer (ilayer)
            layer_info.nwires = 2*(this->ncell0 + this->ncell_increment*Get_nsuperlayer_minus_1(ilayer) );

            // the previous layer info is needed to calculate parameters of current layer
            const auto& previousLayer = this->database.at(ilayer-1);

            //calculate height_z0, radius_sw_z0, stereo_sw_z0
            {
                double h  = previousLayer.height_z0;
                double ru = previousLayer.radius_fuw_z0;
                double rd = previousLayer.radius_fdw_z0;

                if(0 == Get_nsuperlayer_minus_1(ilayer))
                    layer_info.height_z0 = h*ru/rd;
                else
                    // we calculate the height (=width) as  h =  (ru + 1/2 * h) * (2pi / n_cells)
                    layer_info.height_z0 = ru*dd4hep::twopi/(0.5*layer_info.nwires - dd4hep::pi);

                layer_info.radius_sw_z0 = 0.5*layer_info.height_z0 + ru;
                layer_info.stereo_sw_z0  = StereoSign(layer_info) * this->stereoangle_z0(layer_info.radius_sw_z0);
            }

            //calculate radius_fdw_z0, radius_fuw_z0, width_z0
            layer_info.radius_fdw_z0 = previousLayer.radius_fuw_z0;
            layer_info.radius_fuw_z0 = previousLayer.radius_fuw_z0 + layer_info.height_z0;
            layer_info.width_z0 = dd4hep::twopi*layer_info.radius_sw_z0/(0.5*layer_info.nwires);

            // according to expert prescription, width_z0 == height_z0
            if(fabs(layer_info.width_z0 - layer_info.height_z0)>1e-4)
                throw std::runtime_error("fabs(l.width_z0 - l.height_z0)>1e-4");

            this->database.emplace(ilayer, layer_info);
        }

        std::cout << "\t+ Total size of DCH database = " << database.size() << std::endl;
        return;
    }

    // DCH database dumper
    void ShowDatabase(std::ostream & oss) const override final {
        oss << "\n";
        oss << "Global parameters of DCH:\n";
        oss << "\tGas, half length/mm = " << Lhalf/dd4hep::mm << '\n';
        oss << "\tGas, radius in/mm  = " << rin/dd4hep::mm << '\n';
        oss << "\tGas, radius out/mm = " << rout/dd4hep::mm<< '\n';
        oss << "\tGuard, radius in(z=0)/mm  = " << guard_inner_r_at_z0/dd4hep::mm << '\n';
        oss << "\tGuard, radius out(z=L/2)/mm = " << guard_outer_r_at_zL2/dd4hep::mm << '\n';
        oss << "\n";
        oss << "\tTwist angle (2*alpha) / deg = " << twist_angle/dd4hep::deg << '\n';
        oss << "\n";
        oss << "\tN superlayers = " << nsuperlayers << '\n';
        oss << "\tN layers per superlayer = " << nlayersPerSuperlayer << '\n';
        oss << "\tN layers = " << nlayers << '\n';
        oss << "\n";
        oss << "\tN cells layer1 = " << ncell0 << '\n';
        oss << "\tN cells increment per superlayer = " << ncell_increment << '\n';
        oss << "\tN cells per sector = " << ncell_per_sector << '\n';
        oss << "\n";
        oss << "Layer parameters of DCH:\n";
        oss
                << "\t" << "layer"
                << "\t" << "nwires"
                << "\t" << "height_z0/mm"
                << "\t" << "width_z0/mm"
                << "\t" << "radius_fdw_z0/mm"
                << "\t" << "radius_sw_z0/mm"
                << "\t" << "radius_fuw_z0/mm"
                << "\t" << "stereo_sw_z0/deg"
                << "\t" << "stereoangle_z0/deg"
                << "\t" << "Pitch_z0/mm"
                << "\t" << "radius_sw_zLhalf/mm"
                << "\t" << "WireLength/mm"
                << "\n" << std::endl;

        if( this->IsDatabaseEmpty() )
        {
            oss << "\nDatabase empty\n";
            return;
        }

        for(const auto& [nlayer, l]  : database )
        {
            oss
                    << "\t" << l.layer
                    << "\t" << l.nwires
                    << "\t" << l.height_z0/dd4hep::mm
                    << "\t" << l.width_z0/dd4hep::mm
                    << "\t" << l.radius_fdw_z0/dd4hep::mm
                    << "\t" << l.radius_sw_z0/dd4hep::mm
                    << "\t" << l.radius_fuw_z0/dd4hep::mm
                    << "\t" << l.stereo_sw_z0/dd4hep::deg
                    << "\t" << this->StereoSign(l)*this->stereoangle_z0(l.radius_sw_z0)/dd4hep::deg
                    << "\t" << this->Pitch_z0(l, l.radius_sw_z0)/dd4hep::mm
                    << "\t" << this->Radius_zLhalf(l.radius_sw_z0)/dd4hep::mm
                    << "\t" << this->WireLength(l.layer,l.radius_sw_z0)/dd4hep::mm
                    << "\n" << std::endl;
        }
        return;
    }
};  // end of DCH_info_struct
typedef StructExtension<DCH_info_struct> DCH_info;



//===========================================================
// StrawTubeTracker info
//===========================================================

/* Data structure to store StrawTubeTracker geometry parameters,
 * implementing virtual methods of WireTracker_info_struct.
 */
struct STT_info_struct : WireTracker_info_struct {

    //--------------------------------------------------------------
    //  STT-specific global parameters

    /// number of sectors
    std::vector<int> nsectors {};
    /// Half angular distance between two consecutive tubes within the same layer, for each superlayer.
    /// This also corresponds to the delta_phi staggering of tubes across layers.
    std::vector<angle_t> delta_phi {};
    /// stereo angle of superlayers
    std::vector<angle_t> stereo {};
    // radius of sensitive in wire in innermost layers of each superlayer
    std::vector<length_t> innermost_radius {};

    //--------------------------------------------------------------
    // STT geo helpers

    //--------------------------------------------------------------
    // STT geo interface implementation

    /// Global phi positioning in STT,
    /// sum of sector phi + phi offset due to layer staggering + tube phi within the layer
    angle_t Get_cell_phi_angle(int superlayer, int ilayer, int sector, int tube) const override final {
        //const auto & l = this->database.at(ilayer);
        angle_t phi_start = sector * (2 * dd4hep::pi  / this->nsectors.at(superlayer));
        angle_t phi_rel = (ilayer + 2 * tube) * this->delta_phi.at(superlayer) * pow(-1, superlayer);
        return phi_start + phi_rel;
    }

    //--------------------------------------------------------------
    // Database construction

    // STT database builder
    void BuildLayerDatabase() override final {
        // do not fill twice the database
        if( not this->IsDatabaseEmpty() ) return;

        /// nlayers = nsuperlayers * nlayersPerSuperlayer
        /// default: 80 = 8 * 10
        this->nlayers = this->nsuperlayers * this->nlayersPerSuperlayer;

        int layer = 1; //
        for(int superlayer = 0; superlayer < this->nsuperlayers; ++superlayer) {

            // get some superlayer parameters
            length_t R = this->innermost_radius.at(superlayer);
            length_t cell_diameter = 2 * R * sin(this->delta_phi.at(superlayer));
            angle_t sl_stereo = this->stereo.at(superlayer);

            for(int ilayer = 0; ilayer < this->nlayersPerSuperlayer; ++ilayer) {
                // initialize empty object
                LayerInfo layer_info;

                // set parameters
                layer_info.layer = layer;
                // in consecutive layers, minimal radial distance between staggered tubes centers
                // is cell_diameter * sqrt(3)/2 (=0.866)
                layer_info.radius_sw_z0 = R + (ilayer * cell_diameter * 0.866);
                layer_info.stereo_sw_z0 = sl_stereo;

                this->database.emplace(layer, layer_info);
                layer++;
            }
        }

        std::cout << "\t+ Total size of DCH database = " << database.size() << std::endl;
        return;
    }

    // STT database dumper
    void ShowDatabase(std::ostream & oss) const override final {
        oss << "\n";
        oss << "Global parameters of DCH:\n";
        oss << "\tGas, half length/mm = " << Lhalf/dd4hep::mm << '\n';
        oss << "\tGas, radius in/mm  = " << rin/dd4hep::mm << '\n';
        oss << "\tGas, radius out/mm = " << rout/dd4hep::mm<< '\n';
        oss << "\n";
        oss << "\tN superlayers = " << nsuperlayers << '\n';
        oss << "\tN layers per superlayer = " << nlayersPerSuperlayer << '\n';
        oss << "\tN layers = " << nlayers << '\n';
        oss << "\n";
        oss << "Layer parameters of DCH:\n";
        oss
                << "\t" << "layer"
                << "\t" << "radius_sw_z0/mm"
                << "\t" << "stereo_sw_z0/mm"
                << "\n" << std::endl;

        if( this->IsDatabaseEmpty() )
        {
            oss << "\nDatabase empty\n";
            return;
        }

        for(const auto& [nlayer, l]  : database )
        {
            oss
                    << "\t" << l.layer
                    << "\t" << l.radius_sw_z0/dd4hep::mm
                    << "\t" << l.stereo_sw_z0/dd4hep::deg
                    << "\n" << std::endl;
        }
        return;
    }

};  // end of STT_info_struct
typedef StructExtension<STT_info_struct> STT_info;


}} // end namespace dd4hep::rec::

#endif // WIRE_TRACKER_INFO_H
