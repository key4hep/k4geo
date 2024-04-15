#ifndef DCH_INFO_H_INCLUDED
#define DCH_INFO_H_INCLUDED

#include "TMath.h"

#include "DDRec/DetectorData.h"
#include <map>

namespace dd4hep {  namespace rec {


/* Data structure to store DCH excel table
 * to be filled before building the geometry!
 *
 * Global variables are members of this class
 * Helper class DCH_info contains one entry in the table
 * which corresponds to one layer
 *
 * @author A. Tolosa Delgado, CERN
 * @date April 2024
 * @version Drift chamber v2
 *
 */
struct DCH_info_struct
{
public:
    // use alias of types to show more clear what the variable is
    // if everything is double, the code is not readable
    /// type for layer number
    using DCH_layer = int;
    /// tpye for lengths
    using MyLength_t = double;
    /// tpye for angles
    using MyAngle_t = double;
    /// Half length of the active volume
    MyLength_t dch_Lhalf = {0};
    /// Inner radius of the active volume
    MyLength_t dch_rin = {0};
    /// Outer radius of the active volume
    MyLength_t dch_rout = {0};

    /// Inner guard wires radius
    MyLength_t dch_rin_z0_guard = {0};
    /// Outer guard wires radius
    MyLength_t dch_rout_z0_guard = {0};

    /// number of cells of first layer
    int dch_ncell0 = {0};
    /// increment the number of cells for each superlayer as:
    ///   ncells(ilayer) = dch_ncell0 + increment*superlayer(ilayer)
    ///   See DCH_info::Get_nsuperlayer_minus_1(ilayer)
    int dch_ncell_increment = {0};

    /// cells within the same layer may be grouped into sectors, not in use atm
    int dch_ncell_per_sector = {0};

    /// input number of layers in each superlayer
    int dch_nlayersPerSuperlayer = {0};
    /// input number of superlayers
    /// superlayer is an abstract level of grouping layers used to
    /// parametrize the increment of cells in each layer
    int dch_nsuperlayers = {0};
    /// Calculated as dch_nlayersPerSuperlayer * dch_nsuperlayers
    int dch_nlayers = {0};

    /// global twist angle
    /// alternating layers will change its sign
    MyAngle_t  dch_twist_angle = {0};

    /// Cell width for the first layer
    double dch_first_width = {0};
    /// Cell radius for the first layer
    MyLength_t dch_first_sense_r = {0};

    void Set_lhalf(MyLength_t _dch_Lhalf){dch_Lhalf=_dch_Lhalf;};
    void Set_rin  (MyLength_t _dch_rin  ){dch_rin  = _dch_rin; };
    void Set_rout (MyLength_t _dch_rout ){dch_rout = _dch_rout;};

    void Set_guard_rin (MyLength_t _dch_rin_z0_guard ){dch_rin_z0_guard  = _dch_rin_z0_guard; };
    void Set_guard_rout(MyLength_t _dch_rout_z0_guard){dch_rout_z0_guard = _dch_rout_z0_guard;};

    void Set_ncell0              (int _ncell0              ){dch_ncell0               = _ncell0;              };
    void Set_ncell_increment     (int _ncell_increment     ){dch_ncell_increment      = _ncell_increment;     };

    void Set_nlayersPerSuperlayer(int _nlayersPerSuperlayer){dch_nlayersPerSuperlayer = _nlayersPerSuperlayer;};
    void Set_nsuperlayers        (int _nsuperlayers        ){dch_nsuperlayers         = _nsuperlayers;        };

    void Set_ncell_per_sector(int _ncell_per_sector){dch_ncell_per_sector = _ncell_per_sector;};

    void Set_twist_angle (MyLength_t _dch_twist_angle ){dch_twist_angle = _dch_twist_angle;};

    void Set_first_width  (double _first_width  ){dch_first_width   = _first_width;   };
    void Set_first_sense_r(double _first_sense_r){dch_first_sense_r = _first_sense_r; };


    /// Get number of cells in a given layer
    ///   ncells = 2x number of wires
    inline int Get_ncells(int ilayer){return database.at(ilayer).nwires/2;}

    /// Get phi width for the twisted tube and the step (phi distance between cells)
    inline MyAngle_t Get_phi_width(int ilayer){return (TMath::TwoPi()/Get_ncells(ilayer))*dd4hep::rad;}

    ///
    // TODO: check if staggering is just + 0.25*cell_phi_width*(ilayer%2);
    // phi positioning, adding offset for odd ilayers
    MyAngle_t Get_cell_phi_angle(int ilayer, int nphi){ return (Get_phi_width(ilayer) * (nphi + 0.25*(ilayer%2)));}

    /// calculate superlayer for a given ilayer.
    /// WARNING: division of integers on purpose!
    int Get_nsuperlayer_minus_1(int ilayer){ return int((ilayer-1)/dch_nlayersPerSuperlayer);}

    /// Calculate radius at z=L/2 given at z=0
    MyLength_t Radius_zLhalf(MyLength_t r_z0) const {
        return r_z0/cos(dch_twist_angle/2/dd4hep::rad);
    }

    /// tan(stereoangle) = R(z=0)   / (L/2) * tan( twist_angle/2)
    inline MyAngle_t stereoangle_z0(MyLength_t r_z0) const {
        return atan( r_z0/dch_Lhalf*tan(dch_twist_angle/2));
    }

    /// tan(stereoangle) = R(z=L/2) / (L/2) * sin( twist_angle/2)
    inline MyAngle_t stereoangle_zLhalf(MyLength_t r_zLhalf) const {
        return atan( r_zLhalf/dch_Lhalf*sin(dch_twist_angle/2));
    }

    /// WireLength = 2*dch_Lhalf/cos(atan(Pitch_z0(r_z0)/(2*dch_Lhalf)))/cos(stereoangle_z0(r_z0))
    inline MyLength_t WireLength(int nlayer, MyLength_t r_z0) const {
        auto Pitch_z0 = database.at(nlayer).Pitch_z0(r_z0);
        return  2*dch_Lhalf/cos(atan(Pitch_z0/(2*dch_Lhalf)))/cos(stereoangle_z0(r_z0)) ;
    };

    /// Internal helper struct for defining the layer layout
    struct DCH_info_layer
    {
        // each member corresponds to a column in the excel table
        /// layer number
        DCH_layer layer = {0};
        /// 2x number of wires in that layer
        int nwires = {0};
        /// cell parameter
        double height_z0 = {0.};
        /// cell parameter
        double width_z0 = {0.};

        /// radius (cylindral coord) of sensitive wire
        MyLength_t radius_sw_z0 = {0.};

        /// radius (cylindral coord) of 'down' field wires
        MyLength_t radius_fdw_z0 = {0.};

        /// radius (cylindral coord) of 'up' field wires
        MyLength_t radius_fuw_z0 = {0.};

        // some quantityies are derived from previous ones
        // see excel and paper for equations
        ///  stereo angle is positive for odd layer number
        bool IsStereoPositive() const {
            return (1 == layer%2);
        }
        /// calculate sign based on IsStereoPositive
        inline int StereoSign() const {
            return (IsStereoPositive()*2 -1 );
        }

        /// separation between wires (along the circle)
        MyLength_t Pitch_z0(MyLength_t r_z0) const {
            return TMath::TwoPi()*r_z0/nwires;
        };

    };
    /// map to store excel table
    std::map<DCH_layer, DCH_info_layer> database;
    bool IsDatabaseEmpty() const { return (0 == database.size() ); }

    void BuildLayerDatabase();
    void Show_DCH_info_database(std::ostream& io) const;


};
typedef StructExtension<DCH_info_struct> DCH_info ;
std::ostream& operator<<( std::ostream& io , const DCH_info& d ){d.Show_DCH_info_database(io); return io;}

void DCH_info_struct::BuildLayerDatabase()
{
    // do not fill twice the database
    if( not this->IsDatabaseEmpty() ) return;

    auto ff_check_positive_parameter = [](double p, std::string pname) -> void {
        if(p<=0)throw std::runtime_error(Form("DCH: %s must be positive",pname.c_str()));
        return;
    };
    ff_check_positive_parameter(this->dch_rin  ,"inner radius");
    ff_check_positive_parameter(this->dch_rout ,"outer radius");
    ff_check_positive_parameter(this->dch_Lhalf,"half length" );

    ff_check_positive_parameter(this->dch_rin_z0_guard ,"inner radius of guard wires" );
    ff_check_positive_parameter(this->dch_rout_z0_guard,"outer radius of guard wires" );


    ff_check_positive_parameter(this->dch_ncell0,"ncells in the first layer" );
    ff_check_positive_parameter(this->dch_ncell_increment,"ncells increment per superlayer" );
    ff_check_positive_parameter(this->dch_ncell_per_sector,"ncells per sector" );

    // if dch_ncell_per_sector is not divisor of dch_ncell0 and dch_ncell_increment
    // trow an error
    if( 0 != (dch_ncell0 % dch_ncell_per_sector) || 0 != (dch_ncell_increment % dch_ncell_per_sector) )
        throw std::runtime_error("dch_ncell_per_sector is not divisor of dch_ncell0 or dch_ncell_increment");

    ff_check_positive_parameter(this->dch_nsuperlayers,"number of superlayers" );
    ff_check_positive_parameter(this->dch_nlayersPerSuperlayer,"number of layers per superlayer" );

    /// nlayers = nsuperlayers * nlayersPerSuperlayer
    /// default: 112 = 14 * 8
    this->dch_nlayers = this->dch_nsuperlayers * this->dch_nlayersPerSuperlayer;

    ff_check_positive_parameter(this->dch_first_width,"width of first layer cells" );
    ff_check_positive_parameter(this->dch_first_sense_r,"radius of first layer cells" );


    // intialize layer 1 from input parameters
    {
        DCH_info_layer layer1_info;
        layer1_info.layer         = 1;
        layer1_info.nwires        = 2*this->dch_ncell0;
        layer1_info.height_z0     = dch_first_width;
        layer1_info.radius_sw_z0  = dch_first_sense_r;
        layer1_info.radius_fdw_z0 = dch_first_sense_r - 0.5*dch_first_width;
        layer1_info.radius_fuw_z0 = dch_first_sense_r + 0.5*dch_first_width;
        layer1_info.width_z0      = TMath::TwoPi()*dch_first_sense_r/this->dch_ncell0;

        this->database.emplace(layer1_info.layer, layer1_info);
    }

    // some parameters of the following layer are calculated based on the previous ones
    // the rest are left as methods of DCH_info or DCH_info_layer class
    // loop over all layers, skipping the first one
    for(int ilayer = 2; ilayer<= this->dch_nlayers; ++ilayer)
    {
        // initialize here an object that will contain one row from the excel table
        DCH_info_layer layer_info;

        // the loop counter actually corresponds to the layer number
        layer_info.layer = ilayer;
        // nwires is twice the number of cells in this particular layer (ilayer)
        layer_info.nwires = 2*(this->dch_ncell0 + this->dch_ncell_increment*Get_nsuperlayer_minus_1(ilayer) );

        // the previous layer info is needed to calculate parameters of current layer
        auto previousLayer = this->database.at(ilayer-1);

        //calculate height_z0, radius_sw_z0
        {
            double h  = previousLayer.height_z0;
            double ru = previousLayer.radius_fuw_z0;
            double rd = previousLayer.radius_fdw_z0;

            if(0 == Get_nsuperlayer_minus_1(ilayer))
                layer_info.height_z0 = h*ru/rd;
            else
                layer_info.height_z0 = TMath::TwoPi()*ru/(0.5*layer_info.nwires - TMath::Pi());

            layer_info.radius_sw_z0 = 0.5*layer_info.height_z0 + ru;
        }

        //calculate radius_fdw_z0, radius_fuw_z0, width_z0
        layer_info.radius_fdw_z0 = previousLayer.radius_fuw_z0;
        layer_info.radius_fuw_z0 = previousLayer.radius_fuw_z0 + layer_info.height_z0;
        layer_info.width_z0 = TMath::TwoPi()*layer_info.radius_sw_z0/(0.5*layer_info.nwires);

        // according to excel, width_z0 == height_z0
        if(fabs(layer_info.width_z0 - layer_info.height_z0)>1e-4)
            throw std::runtime_error("fabs(l.width_z0 - l.height_z0)>1e-4");



        this->database.emplace(ilayer, layer_info);
    }

    std::cout << "\t+ Total size of DCH database = " << database.size() << std::endl;
    return;
}

void DCH_info_struct::Show_DCH_info_database(std::ostream & oss) const
{
    if( this->IsDatabaseEmpty() )
        throw std::runtime_error("Can not show empty database!");

    oss << "\n";
    oss << "Global parameters of DCH:\n";
    oss << "\tHalf length/mm = " << dch_Lhalf/dd4hep::mm << '\n';
    oss << "\tRadius in/mm  = " << dch_rin/dd4hep::mm << '\n';
    oss << "\tRadius out/mm = " << dch_rout/dd4hep::mm<< '\n';
    oss << "\tRadius guard in/mm  = " << dch_rin_z0_guard/dd4hep::mm << '\n';
    oss << "\tRadius guard out/mm = " << dch_rout_z0_guard/dd4hep::mm << '\n';
    oss << "\n";
    oss << "\tTwis angle (2*alpha) / deg = " << dch_twist_angle/dd4hep::deg << '\n';
    oss << "\n";
    oss << "\tN superlayers = " << dch_nsuperlayers << '\n';
    oss << "\tN layers per superlayer = " << dch_nlayersPerSuperlayer << '\n';
    oss << "\tN layers = " << dch_nlayers << '\n';
    oss << "\n";
    oss << "\tN cells layer1 = " << dch_ncell0 << '\n';
    oss << "\tN cells increment per superlayer = " << dch_ncell_increment << '\n';
    oss << "\tN cells per sector = " << dch_ncell_per_sector << '\n';
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
            << "\t" << "stereoangle_z0/deg"
            << "\t" << "Pitch_z0/mm"
            << "\t" << "radius_sw_zLhalf/mm"
            << "\t" << "WireLength/mm"
            << "\n" << std::endl;

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
                << "\t" << l.StereoSign()*this->stereoangle_z0(l.radius_sw_z0)/dd4hep::deg
                << "\t" << l.Pitch_z0(l.radius_sw_z0)/dd4hep::mm
                << "\t" << this->Radius_zLhalf(l.radius_sw_z0)/dd4hep::mm
                << "\t" << this->WireLength(l.layer,l.radius_sw_z0)/dd4hep::mm
                << "\n" << std::endl;
    }
    return;
}
}} // end namespace dd4hep::rec::

#endif // DCH_INFO_H_INCLUDED
