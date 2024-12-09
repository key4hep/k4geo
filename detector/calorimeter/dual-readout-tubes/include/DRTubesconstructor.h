#ifndef DRconstructor_H
#define DRconstructor_H 1

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Objects.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

#include "DRutils.h"

using namespace dd4hep;

namespace DDDRCaloTubes {

class DRTubesconstructor {
public:
    // Constructor
    DRTubesconstructor(Detector* description,
                  xml_h& entities,
                  SensitiveDetector* sens);

    // Destructor
    ~DRTubesconstructor() {}

    // Function to calculate some (but not all) tower size parameters
    void calculate_tower_parameters();

    // Function to calculate all parameters which only depend on the tower phi
    void calculate_phi_parameters();

    // Function to calculate all parameters which depend on theta (both tower theta and the theta that has been covered by placing towers)
    // Since this is different for each tower, this function is called multiple times
    void calculate_theta_parameters();

    // Function to create tube volumes and store the in a map, if they don't already exist
    void assert_tube_existence(int key, bool cher);

    // Function to calculate the width of the support Trap volume at a given y and z
    // with the option to toggle calculation for the backface or frontface width
    double calculate_trap_width(double given_y, double given_z, bool backface = false);
    // Same as trap width, but for the tower volume (so dimension of "air" insdie the support structure which is the _trap_ volume)
    double calculate_tower_width(int given_row, bool backface = true);

    // Function to calculate the tube lengths and place them to create the actual tower (not the air)
    void assemble_tower(Volume& tower_air_volume);

    // Mostly just a wrapper function
    void construct_tower_trapezoid(Volume& trap_volume);

    // Function to calculate the position of the tower inside the stave
    void calculate_tower_position();

    // Function to construct the trapezoidal support structure for the tower in which fibres are placed
    void construct_tower(Volume& trap_volume);

    void increase_covered_theta(const double& delta_theta) {m_covered_theta += delta_theta;}

    // Function to place the tower in the stave volume
    void place_tower(Volume& stave_volume,
                    Volume& tower_volume,
                    unsigned int layer);

    // Overarching function to construct the calorimeter
    void construct_calorimeter(Volume& calorimeter_volume);

private:
    Detector* m_description;
    xml_h m_entities;
    SensitiveDetector* m_sens;

    // Calorimeter parameters
    double m_calo_inner_r;
    double m_calo_outer_r;
    double m_calo_inner_half_z;

    double m_barrel_endcap_angle; // calculated from m_calo_inner_half_z and m_calo_inner_r

    // Tube parameters
    double    m_capillary_outer_r;
    double    m_scin_clad_outer_r;
    double    m_scin_core_outer_r;
    double    m_cher_clad_outer_r;
    double    m_cher_core_outer_r;
    Material m_capillary_material;
    Material m_scin_clad_material;
    Material m_scin_core_material;
    Material m_cher_clad_material;
    Material m_cher_core_material;
    std::string m_capillary_visString;
    std::string m_scin_clad_visString;
    std::string m_scin_core_visString;
    std::string m_cher_clad_visString;
    std::string m_cher_core_visString;
    bool m_capillary_isSensitive;
    bool m_scin_clad_isSensitive;
    bool m_scin_core_isSensitive;
    bool m_cher_clad_isSensitive;
    bool m_cher_core_isSensitive;

    // Maps to store the tube volumes, so that one volume can be used multiple times
    // The key to the map is an indicator of the tube length (multiple of the tolerance)
    std::unordered_map<int, Volume> m_scin_tube_volume_map;
    std::unordered_map<int, Volume> m_cher_tube_volume_map;

    // Tolerance for which new tube volumes are created
    // e.g 1mm, then all tube lengths are rounded (down) to the nearest mm
    double m_tolerance;


    // Constants used through the function (calculated from other parameters)
    double m_capillary_diameter; // calculated from m_capillary_outer_r
    double m_V; // Vertical spacing for pointy top oriented tubes

    // Tower angle parameters (and derived parameters)
    double m_tower_theta;
    double m_tower_phi;
    
    double m_tower_half_phi;
    double m_tower_tan_half_phi; // calculated from m_tower_phi
    double m_tower_half_length; // calculated from m_calo_inner_r and m_calo_outer_r and m_trap_half_length
    double m_tower_tan_theta;

    // Tower size parameters depending on phi
    unsigned int m_num_phi_towers;       // number of towers in phi direction
    // Tower widths at four edges
    double m_tower_frontface_rightangleedge_x;
    double m_tower_frontface_thetaangleedge_x;
    double m_tower_backface_rightangleedge_x;
    double m_tower_backface_thetaangleedge_x;
    // Angle between the parallel edges of the front/back face
    double m_angle_edges_x;

    // Tower size parameters derived from theta
    double m_tower_frontface_y;
    double m_tower_backface_y;
    double m_tower_polar_angle;
    double m_tower_azimuthal_angle;

    // Trapezoid support parameters
    double m_stave_half_length;
    double m_trap_wall_thickness_sides;
    double m_trap_wall_thickness_front;
    double m_trap_wall_thickness_back;
    double m_trap_frontface_rightangleedge_x;       // width for frontface
    double m_trap_frontface_thetaangleedge_x;
    double m_trap_backface_rightangleedge_x;        // width for backface
    double m_trap_backface_thetaangleedge_x;
    double m_trap_frontface_y;      // height for frontface
    double m_trap_backface_y;       // height for backface
    double m_trap_azimuthal_angle;  // azimuthal angle for the trapezoid
    double m_trap_polar_angle;      // polar angle for the trapezoid
    double m_trap_half_length;      // half length for the trapezoid
    Material m_trap_material;
    std::string m_trap_visString;


    // Construction parameters (which change for each tower)
    double m_covered_theta;
    double m_back_shift;
    Position m_tower_position;

    Material m_air;
    std::string m_air_visString;

};

} // namespace DDDRCaloTubes

#endif // DRCONSTRUCTOR_H
