#include "DRTubesconstructor.h"

#include <TMatrixD.h>

using namespace dd4hep;
using namespace DDDRCaloTubes;

DDDRCaloTubes::DRTubesconstructor::DRTubesconstructor(Detector* description,
                                            xml_h& entities,
                                            SensitiveDetector* sens):
                                            m_entities(entities)
{
    // Initialising the all the calorimeter and tower variables

    m_description = description;
    m_sens = sens;

    xml_dim_t x_dim = ((xml_det_t) entities).dimensions();

    // Calorimeter parameters
    m_calo_inner_r      = x_dim.inner_radius();
    m_calo_outer_r      = x_dim.outer_radius();
    m_calo_inner_half_z = x_dim.z_length();

    
    // Trap parameters
    xml_comp_t  x_trap = entities.child(_Unicode(trap));
    xml_comp_t  x_trap_support = x_trap.child(_Unicode(support));
    m_trap_wall_thickness_front = x_trap_support.depth();
    m_trap_wall_thickness_sides = x_trap_support.width();
    m_trap_wall_thickness_back  = x_trap_support.z2();
    m_trap_material             = m_description->material(x_trap_support.materialStr());
    m_trap_visString            = x_trap_support.visStr();


    // Tube parameters
    xml_comp_t  x_tube = entities.child(_Unicode(tube));

    xml_comp_t  x_capillary = x_tube.child(_Unicode(capillary));
    m_capillary_material    = m_description->material(x_capillary.materialStr()); 
    m_capillary_outer_r     = x_capillary.outer_r();
    m_capillary_visString   = x_capillary.visStr();
    m_capillary_isSensitive = x_capillary.isSensitive();
    m_tolerance             = x_capillary.threshold(50*um);

    xml_comp_t  x_scin_clad = x_tube.child(_Unicode(scin_clad));
    m_scin_clad_material    = m_description->material(x_scin_clad.materialStr()); 
    m_scin_clad_outer_r     = x_scin_clad.outer_r();
    m_scin_clad_visString   = x_scin_clad.visStr();
    m_scin_clad_isSensitive = x_scin_clad.isSensitive();

    xml_comp_t  x_scin_core = x_tube.child(_Unicode(scin_core));
    m_scin_core_material    = m_description->material(x_scin_core.materialStr()); 
    m_scin_core_outer_r     = x_scin_core.outer_r();
    m_scin_core_visString   = x_scin_core.visStr();
    m_scin_core_isSensitive = x_scin_core.isSensitive();

    xml_comp_t  x_cher_clad = x_tube.child(_Unicode(cher_clad));
    m_cher_clad_material    = m_description->material(x_cher_clad.materialStr()); 
    m_cher_clad_outer_r     = x_cher_clad.outer_r();
    m_cher_clad_visString   = x_cher_clad.visStr();
    m_cher_clad_isSensitive = x_cher_clad.isSensitive();

    xml_comp_t  x_cher_core = x_tube.child(_Unicode(cher_core));
    m_cher_core_material    = m_description->material(x_cher_core.materialStr()); 
    m_cher_core_outer_r     = x_cher_core.outer_r();
    m_cher_core_visString   = x_cher_core.visStr();
    m_cher_core_isSensitive = x_cher_core.isSensitive();

    // Some safety checks for the entered tube parameters
    if (m_capillary_outer_r <= 0.0*mm) throw std::runtime_error("Capillary radius needs to be larger than 0");
    if (m_capillary_outer_r < m_scin_clad_outer_r || m_capillary_outer_r < m_cher_clad_outer_r) throw std::runtime_error("Capillary radius needs to be larger than scintillation cladding and cherenkov cladding radii");
    if (m_scin_clad_outer_r < m_scin_core_outer_r) throw std::runtime_error("Scintillation cladding radius needs to be larger than scintillation core radius");
    if (m_cher_clad_outer_r < m_cher_core_outer_r) throw std::runtime_error("Cherenkov cladding radius needs to be larger than cherenkov core radius");


    // Tower parameters
    m_tower_theta = x_dim.deltatheta();
    m_tower_phi   = x_dim.deltaphi();


    // Construction parameters which change for each placed tower in a stave
    m_covered_theta = 0.0*deg;
    m_back_shift = 0.0*mm;

    // Calculate all the derived parameters which are not taken from the xml file
    this->calculate_tower_parameters();
    this->calculate_phi_parameters();

    m_tower_tan_theta = 0.0;

    m_air = m_description->material("Air");

    xml_comp_t  x_air = x_trap.child(_Unicode(air));
    m_air    = m_description->material(x_air.materialStr()); 
    m_air_visString     = x_air.visStr();

}

// Function to calculate all tower parameters which are derived from user given values
void DDDRCaloTubes::DRTubesconstructor::calculate_tower_parameters()
{
    // Angle where endcap and barrel meet
    m_barrel_endcap_angle = std::atan2(m_calo_inner_half_z, m_calo_inner_r);

    m_capillary_diameter = 2*m_capillary_outer_r;

    // Constants used through the function
    double D = 4.0*m_capillary_outer_r/sqrt(3.0); // Long diagonal of hexagaon with capillary_outer_r as inradius
    m_V = 3.0*D/4.0;                        // Vertical spacing for pointy top oriented tubes

    // Some values which are used multiple times are calculated here once
    m_tower_half_phi = m_tower_phi/2.0; // Half of the tower phi angle
    m_tower_tan_half_phi = std::tan(m_tower_half_phi); // Needed several times, calculate once here


    m_stave_half_length  = (m_calo_outer_r*std::cos(m_tower_half_phi) - m_calo_inner_r)/2; // Trapezoid half length

    // Protection against tilted towers in theta going past the outer radius of the calorimeter
    double protect_covered_z = std::tan(m_tower_theta)*m_calo_inner_r;
    double protect_tower_z = std::tan(2*m_tower_theta)*m_calo_inner_r - protect_covered_z;
    double protect_back_shift = std::sin(m_tower_theta)*protect_tower_z;

    // "_trap_" always* refers to the Trap volume including the support structure
    m_trap_half_length = m_stave_half_length/std::cos(m_tower_theta) - protect_back_shift/2;
    
    // "_tower_" always* refers to the collection of tubes forming the tower
    // "Tower" volume encapsulates the air, which is "hollowing out" the trapezoid volume
    m_tower_half_length = m_trap_half_length - m_trap_wall_thickness_front/2.0 - m_trap_wall_thickness_back/2.0; // Tower half length

    // * "always" means that I tried to use it consistently, but there might be exceptions (like m_tower_phi&theta refer to the full strucutre, including support)
    // due to the history of the code, where I started out without the support, only the collection of tubes 
}

// Function to calculate tower parameters specifically for phi direction
void DDDRCaloTubes::DRTubesconstructor::calculate_phi_parameters()
{
    double num_phi_towers_d = 360.0*deg/m_tower_phi;
    if (num_phi_towers_d < 0.0) throw std::runtime_error("Negative tower phi coverage not allowed");
    // Check if num_phi_towers is a whole number
    if (check_for_integer(num_phi_towers_d)) m_num_phi_towers = static_cast<unsigned int>(std::round(num_phi_towers_d));
    else throw std::runtime_error("Not an integer number of towers in phi direction");

    
}

// Function to calculate the size parameters for the trap volume (including support)
// Depends a lot on m_covered_theta, which is increased with each tower placement in theta direction
/*  Sketch of a tower in the z-x (or z-y) plane
     ____________
    |           /       Each tower is a trapezoid with a front face (bottom in the sketch) and a back face (top in the sketch)
    |          /        The front face is narrower than the back face, to ensure that the tower is projective
    |         /         Also each tower has a right angled face with the front and back face (left in the sketch) 
    |        /          and a theta angle side (right in the sketch)
    |       /           Anytime a variable contains "edge" it refers to the edge that would be running in y direction in this sketch with global coordinates
    |      /            (but would be x direction in the local coordinate system of the Trap volume)
 x  |     /
 ^  |    /
 |  |___/
 |
 |____> z

*/
void DDDRCaloTubes::DRTubesconstructor::calculate_theta_parameters()
{
    // How much the front faces of the placed towers cover in z direction in the global coordinate system
    double covered_z = std::tan(m_covered_theta)*m_calo_inner_r;

    // The cummulative theta value which this trapezoid will cover once placed 
    double trap_theta = m_covered_theta+m_tower_theta;
    // Distance the front face of this trapezoid will cover in z 
    double trap_z = std::tan(trap_theta)*m_calo_inner_r - covered_z; 
     // Tower height 
    double trap_frontface_y = std::cos(m_covered_theta)*trap_z;

    if (trap_frontface_y < m_capillary_diameter)
    {
        throw std::runtime_error("Can't construct tower with given tower_theta and calo_inner_radius");
    }

    // Used throughout construction, so calculate once here
    m_tower_tan_theta = std::tan(m_tower_theta);

    // Front and back face of the trapezoid
    m_trap_frontface_y = trap_frontface_y;
    m_trap_backface_y  = trap_frontface_y + 2*m_trap_half_length*m_tower_tan_theta; // calculated based on how much the tower widens in the back

    // Tower sizes based on the wall thickness of the "support" envelope (needs adjustment for angles of tower)
    m_tower_frontface_y = m_trap_frontface_y + m_trap_wall_thickness_front*m_tower_tan_theta - m_trap_wall_thickness_sides*(1+1/std::cos(m_tower_theta));
    m_tower_backface_y  = m_trap_backface_y  - m_trap_wall_thickness_back*m_tower_tan_theta  - m_trap_wall_thickness_sides*(1+1/std::cos(m_tower_theta));
    
    // Distance by which right angle edge of this tower is shifted backwards to ensure inner radius of calorimeter
    m_back_shift = std::tan(m_covered_theta)*m_trap_frontface_y;                      

    // Width in y direction in sketch above (x direction in local coordinates)
    m_trap_frontface_rightangleedge_x = (m_calo_inner_r+std::cos(m_covered_theta)*m_back_shift)*2*m_tower_tan_half_phi;
    m_trap_frontface_thetaangleedge_x = m_calo_inner_r*2*m_tower_tan_half_phi;

    // by how much the backface widens for both tower and trap for the right angled side of the tower
    double tower_backface_phi_increase_rightangleedge = 2*m_trap_half_length*std::cos(m_covered_theta) * 2*m_tower_tan_half_phi;
     // the theta angles side widens differently, because they are not the same length
    double tower_backface_phi_increase_thetaangleedge = 2*m_trap_half_length/std::cos(m_tower_theta)*std::cos(m_covered_theta+m_tower_theta) * 2*m_tower_tan_half_phi;

    m_trap_backface_rightangleedge_x = m_trap_frontface_rightangleedge_x + tower_backface_phi_increase_rightangleedge;
    m_trap_backface_thetaangleedge_x = m_trap_frontface_thetaangleedge_x + tower_backface_phi_increase_thetaangleedge;

    // The effective thickness is used to calculate the size of the "air" tower volume
    // It takes into account, that the wall is tilted in phi direction
    double effective_side_wall_thickness_x = m_trap_wall_thickness_sides/std::cos(m_tower_phi);

    // When looking at the tower from the front (or back), the angle created by a vertical line and the sides of the tower
    m_angle_edges_x = std::atan2((m_trap_backface_rightangleedge_x-m_trap_backface_thetaangleedge_x)/2.0, m_trap_backface_y);
     /*
        ________________
        \              /   So the angle in the bottom left (or right) corner here
         \ _|         /   
       y^ \ |        /  
        |  \|_______/
        |   
        |______> x
    
    */

    // Calculating all tower widths (horizontal lines in sketch above) for both front and back face (in z direction)
    m_tower_frontface_rightangleedge_x = calculate_trap_width(m_trap_wall_thickness_sides, m_trap_wall_thickness_front, false) - 2.0*effective_side_wall_thickness_x;

    m_tower_backface_rightangleedge_x = calculate_trap_width(m_trap_wall_thickness_sides, m_trap_wall_thickness_back, true) - 2.0*effective_side_wall_thickness_x;

    m_tower_frontface_thetaangleedge_x = calculate_trap_width(m_tower_frontface_y+m_trap_wall_thickness_sides, m_trap_wall_thickness_front, false) - 2*effective_side_wall_thickness_x;

    m_tower_backface_thetaangleedge_x = calculate_trap_width(m_tower_backface_y+m_trap_wall_thickness_sides, m_trap_wall_thickness_back, true) - 2*effective_side_wall_thickness_x;

}


// Check if tube of this (half-)length already exists, if not, create it
void DDDRCaloTubes::DRTubesconstructor::assert_tube_existence(int key, bool cher)
{
    std::unordered_map<int, Volume>* tube_volume_map;

    // Select the right volume maps depending on whether we are creating cherenkov or scintillation tubes
    if (cher) {
        tube_volume_map = &m_cher_tube_volume_map;
    } else {
        tube_volume_map = &m_scin_tube_volume_map;
    }
    
    // If this length already exists, return
    if (tube_volume_map->find(key) != tube_volume_map->end()) return;

    // The map key is the multiple of the tolerance forming the tube length
    double length_rounded_down = key*m_tolerance;
    // std::cout << "Creating tube with length " << length_rounded_down/mm << " mm" << std::endl;

    // Creating the volumes
    // Capillary tube
    Tube        capillary_solid(0.0*mm, m_capillary_outer_r, length_rounded_down);
    Volume capillary_volume("capillary", capillary_solid, m_capillary_material);
    if (m_capillary_isSensitive) capillary_volume.setSensitiveDetector(*m_sens);
    capillary_volume.setVisAttributes(*m_description, m_capillary_visString); 

    // Create the right fibres
    if (cher)
    {
        // Cherenkov cladding
        Tube        cher_clad_solid(0.0*mm, m_cher_clad_outer_r, length_rounded_down);
        Volume      cher_clad_volume("cher_clad", cher_clad_solid, m_cher_clad_material);
        if (m_cher_clad_isSensitive) cher_clad_volume.setSensitiveDetector(*m_sens);
        PlacedVolume cher_clad_placed = capillary_volume.placeVolume(cher_clad_volume);
        cher_clad_volume.setVisAttributes(*m_description, m_cher_clad_visString);
        cher_clad_placed.addPhysVolID("clad", 1).addPhysVolID("cherenkov", 1);

        // Chrerenkov core
        Tube        cher_core_solid(0.0*mm, m_cher_core_outer_r, length_rounded_down);
        Volume      cher_core_volume("DRBT_cher_core", cher_core_solid, m_cher_core_material);
        if (m_cher_core_isSensitive) cher_core_volume.setSensitiveDetector(*m_sens);
        PlacedVolume    cher_core_placed = cher_clad_volume.placeVolume(cher_core_volume);
        cher_core_volume.setVisAttributes(*m_description, m_cher_core_visString);
        cher_core_placed.addPhysVolID("core", 1).addPhysVolID("clad", 0);

    } else {
        // Scintillation cladding
        Tube        scin_clad_solid(0.0*mm, m_scin_clad_outer_r, length_rounded_down);
        Volume      scin_clad_volume("scin_clad", scin_clad_solid, m_scin_clad_material);
        if (m_scin_clad_isSensitive) scin_clad_volume.setSensitiveDetector(*m_sens);
        PlacedVolume scin_clad_placed = capillary_volume.placeVolume(scin_clad_volume);
        scin_clad_volume.setVisAttributes(*m_description, m_scin_clad_visString);
        scin_clad_placed.addPhysVolID("clad", 1).addPhysVolID("cherenkov", 0);

        // Scintillation core
        Tube        scin_core_solid(0.0*mm, m_scin_core_outer_r, length_rounded_down);
        Volume      scin_core_volume("DRBT_scin_core", scin_core_solid, m_scin_core_material);
        if (m_scin_core_isSensitive) scin_core_volume.setSensitiveDetector(*m_sens);
        PlacedVolume    scin_core_placed = scin_clad_volume.placeVolume(scin_core_volume);
        scin_core_volume.setVisAttributes(*m_description, m_scin_core_visString);
        scin_core_placed.addPhysVolID("core", 1).addPhysVolID("clad", 0);
    }

    tube_volume_map->insert(std::make_pair(key, capillary_volume));
    
}


// Similar to calculate_tower_width (see abelow), but for the encasing trapezoid support volume
// Width is calculated for any given point in th z-y plane
double DDDRCaloTubes::DRTubesconstructor::calculate_trap_width(double given_y, double given_z, bool backface)
{
    // Calculate width (x_direction) of trapezoid at given y
    // Assuming y=0 corresponds to right angle edge of the tower
    if (given_z>2*m_trap_half_length) throw std::runtime_error("calculate_trap_width: Given z is larger than length of trapezoid");

    // Since the height (in y direction of the trapezoid) changes along z, we need to calculate the maximum y for a given z
    double max_y_for_given_z;
    // given_z can be interpreted as distance from either the front or the back of the tower (in z direction)
    if (backface) max_y_for_given_z = m_trap_backface_y - given_z*m_tower_tan_theta;
    else          max_y_for_given_z = m_trap_frontface_y + given_z*m_tower_tan_theta;
    if (given_y > max_y_for_given_z) throw std::runtime_error("calculate_trap_width: Given y is larger than maximum y for given z");

    // How much the width (x direction) changes due to the y position 
    double delta_y = -2.0*given_y*std::tan(m_angle_edges_x);
    // How much the width (x direction) changes due to the z position
    double delta_z;
    if (backface) delta_z = -2.0*given_z*std::cos(m_covered_theta)*m_tower_tan_half_phi;
    else          delta_z =  2.0*given_z*std::cos(m_covered_theta)*m_tower_tan_half_phi;

    double trap_x;
    if (backface) trap_x = m_trap_backface_rightangleedge_x  + delta_y + delta_z;
    else          trap_x = m_trap_frontface_rightangleedge_x + delta_y + delta_z;

    return trap_x;

}


// Calculate width of tower in x direction at given row of tubes
/*
    Sketch of the tower volume shape in the xy plane
    _____________________
    \                   /
     \-----------------/  Note that the width will not be calculated at the centre of the row
      \               /   but instead, where the tubes first touch the wall
   y^  \             /  
    |   \___________/
    |   
    |______> x
    /
 |/_ z                                          */
double DDDRCaloTubes::DRTubesconstructor::calculate_tower_width(int given_row, bool backface)
{
    // Calculate width (x_direction) of tower at given row
    // Assuming row 0 is at the right angle edge

    // y distance where tube hits side of the wall with given angle between backfaces (y=0 corresponds to m_tower_backface_rightangleedge_x)
    double y = m_capillary_outer_r + given_row*m_V + std::cos(90*deg-m_angle_edges_x)*m_capillary_outer_r;


    double tower_x;
    // Toggle between the back and front in z direction (not shown in sketch)
    if (backface) tower_x = m_tower_backface_rightangleedge_x  - 2.0*y*std::tan(m_angle_edges_x);
    else          tower_x = m_tower_frontface_rightangleedge_x - 2.0*y*std::tan(m_angle_edges_x);

    return tower_x;
}


// Place all tubes which make up the tower
void DDDRCaloTubes::DRTubesconstructor::assemble_tower(Volume& tower_air_volume)
{    
    // Y-distance of rightangle wall from coordinate system origin
    // Used throughout this function
    double tower_centre_r = m_tower_half_length/std::cos(m_tower_polar_angle);
    double tower_centre_half_y = tower_centre_r*std::sin(m_tower_polar_angle)*std::sin(m_tower_azimuthal_angle) + m_tower_frontface_y/2.0;

    
    // Coordinates of corners of trapezoid
    // when looking at the tower from the front with the right angle edge on top
    Position back_upper_right_corner  = Position(m_tower_backface_rightangleedge_x/2.0,  -tower_centre_half_y, m_tower_half_length);
    Position back_lower_right_corner  = Position(m_tower_backface_thetaangleedge_x/2.0,  tower_centre_half_y,  m_tower_half_length);
    Position front_upper_right_corner = Position(m_tower_frontface_rightangleedge_x/2.0, -tower_centre_half_y, -m_tower_half_length);

    // plane equations to calculate length of tubes at the tower sides
    std::vector<double> plane_right_coefficients = get_plane_equation(back_upper_right_corner, back_lower_right_corner, front_upper_right_corner);
    Direction line_direction = Direction(0, 0, -1);

    // width of tower at row 0 (right angle edge)
    // "at row X" meaning the width of the tower at the height of where the tubes first touch the wall in row X
    double tower_x = calculate_tower_width(0);

    // number of total columns of tubes at the back face of the tower
    unsigned int num_back_cols_rightangleedge = 1 + fast_floor((tower_x-2*std::sin(90*deg-m_angle_edges_x)*m_capillary_outer_r)/m_capillary_diameter);

    /*
    Depending on even or odd number of tubes in the row, the column ID looks like this:
         
    Odd: \ O O O O O O O /
        -6 -4 -2 0 2 4 6   central tube has column ID 0, then increasing by two (with negative values on one side)

    Even: \ O O O O O O /
         -5 -3 -1 1 3 5    there is no single "central" tube, so start with one and increase by two

    This way, it's immediately clear if you are in a row with even number of tubes or not.
    Easier for the position reconstruction, since there is an offset between even and odd rows in hexagonal stacking
    Essentially following doubled offset coordinates from https://www.redblobgames.com/grids/hexagons/#coordinates-doubled
    */

    // Safety counter, should be 0 at the end
    int num_bad_rows = 0;

    // How much distance in y is covered already by the tubes (increases in loop)
    // This value is always: how much this row will cover, once it is placed (for checking when we need to start shortening tubes)
    double covered_tower_y = m_capillary_diameter;

    // Number of rows of tubes in the back face of the tower
    unsigned int num_rows = fast_floor((m_tower_backface_y-m_capillary_diameter)/m_V) + 1; 

    std::cout << "TOTAL ROWS = " << num_rows << " COLS = " << num_back_cols_rightangleedge << std::endl;

    // Loop over the rows of tubes in the tower, starting at the right angle edge
    for (unsigned int row = 0; row < num_rows; row++, covered_tower_y+=m_V)
    {

        // Staggering of tubes at the lower edge (theta edge/face)
        double row_staggered_z = 0.0*mm;
        if (covered_tower_y > m_tower_frontface_y) row_staggered_z = (covered_tower_y-m_tower_frontface_y)/m_tower_tan_theta/2.0;

        // Failsafe for tubes which would have 'negative length'
        // Should not happen, but if it does, following rows will also be too short, so can skip the rest
        if (row_staggered_z > m_tower_half_length) {
            num_bad_rows = num_rows-row;
            std::cout << "Encountered bad row at row " << row << std::endl;
            std::cout << "Number of leftover bad rows: " << num_bad_rows << std::endl;
            break;
        }


        // Update the tower width for this row
        tower_x = calculate_tower_width(row);

        // First row is the defining row. It has either even or odd number of tubes.
        // The following rows will need to alternate between odd and even, because of the hexagonal stacking.
        // 
        // Calculate starting column ID based on even or odd row and adapt the covered_tower_x accordingly.
        unsigned int col;
        if (num_back_cols_rightangleedge & 1)   // Uneven number of tubes (so we have a central tube with column ID = 0)
        {
            col = (row&1) ? 1 : 0;              // alternating between 0 and 1 (row index starts at 0, so first is colID = 0)

        } else                                  // Even number of tubes (no central tube, colID starts at 1)
        {
            col = (row&1) ? 0 : 1;
        }

        double covered_tower_x;  // How much the first tube covers in x direction, then increases in loop
        // Start value depends on whether there is a central tube or not:
        covered_tower_x = (col & 1) ? m_capillary_outer_r + m_capillary_outer_r*std::cos(m_angle_edges_x) 
                                    : m_capillary_outer_r*std::cos(m_angle_edges_x);

        // Width of the tower at the front face for this row
        // Used to check when to shorten the tubes, along with covered_tower_x
        double tower_front_x = calculate_tower_width(row, false);

        // We don't calculate how many tubes will fit beforehand, since this varies between rows
        // Instead, place tubes as long as there is space
        // The column placement starts in the middle and goes outwards in both directions
        while (covered_tower_x < tower_x/2.0)
        {
            // TODO: Check what objects can be moved outside of loop (string _name, Tube _solid, etc.)

            // Calculate the position of the tube
            double x = col*m_capillary_outer_r;
            double y = row*m_V + m_capillary_outer_r;

            // To calculate the length of the tubes on the tower sides ("wings"), we use a plane equation to get the point where the tube intersects the wall
            double col_staggered_z = 0.0*mm;
            if (covered_tower_x > tower_front_x/2.0) // Beyond the front face of the tower, the tubes need to be shortened
            {
                // Point of tube where it hits the wall is not the centre, it will be at the outer edge of course
                // But exact position depends on the angle of the walls (m_angle_edges_x)
                Position line_point = Position(x+m_capillary_outer_r*std::cos(m_angle_edges_x), y-tower_centre_half_y+m_capillary_outer_r*std::sin(m_angle_edges_x), m_tower_half_length);
                Position intersection = get_intersection(plane_right_coefficients, line_point, line_direction);
                col_staggered_z = (m_tower_half_length + intersection.z())/2.0;
            }

            // Negative length tubes are not allowed
            // Shouldn't occur, unless I have made a mistake somewhere (this has saved me in the past already)
            if (row_staggered_z > m_tower_half_length) 
            {
                std::cout << "Encountered bad column at (row, col) = (" << row << ", " << col << ")" << std::endl;
                break;
            }

            // If we stagger in both directions, take the shorter length, so the bigger stagger value
            double z = (row_staggered_z > col_staggered_z) ? row_staggered_z : col_staggered_z ;

            double tube_half_length = m_tower_half_length - z;

            // Reference point for tube placement in tower (trapezoid) centre
            auto position = Position(x, y-tower_centre_half_y, z);
            // And mirrored position for the other side of the tower (since the tower is symmetric in phi (left-right), the tubes are identical)
            auto position_mirrored = Position(-x, y-tower_centre_half_y, z);

            // TubeID composed of col in first 16 bits, row in last 16 bits
            int tube_id          = (col << 16) | row;
            int tube_id_mirrored = (-col << 16) | row;


            // Selecting the right fibre to be placed
            bool cher = (row & 1);
            std::unordered_map<int, Volume>* volume_map;
            if (cher) volume_map = &m_cher_tube_volume_map;
            else      volume_map = &m_scin_tube_volume_map;
            // Round length down to next multiple of tolerance
            int key = static_cast<int>(fast_floor(tube_half_length / m_tolerance));
            // Zero or negative length tubes shouldn't occur at this point, but if so, try with the next tube
            if (key < 1) 
            {
                col += 2;
                covered_tower_x += m_capillary_diameter;
                continue;
            }
            this->assert_tube_existence(key, cher);

            // Get the right tube to be placed, including daughters
            Volume capillary_vol_to_be_placed = volume_map->at(key);

            // Place the right side tube
            PlacedVolume    tube_placed = tower_air_volume.placeVolume(capillary_vol_to_be_placed, tube_id, position);
            tube_placed.addPhysVolID("air",0).addPhysVolID("col", col).addPhysVolID("row", row);

            // If column is not the central one, place the mirrored tube on the other side of the tower
            if (col>0)
            {
                PlacedVolume    tube_placed2 = tower_air_volume.placeVolume(capillary_vol_to_be_placed, tube_id_mirrored, position_mirrored);
                tube_placed2.addPhysVolID("air",0).addPhysVolID("col", -col).addPhysVolID("row", row);
            }
            
            col += 2;
            covered_tower_x += m_capillary_diameter;
        }

    }

}


// Function to calculate the position of the tower in stave
void DDDRCaloTubes::DRTubesconstructor::calculate_tower_position()
{
    // Since the Trapezoids are defined including polar and azimuthal angles, we need to convert between cartesian and polar coordinates
    double trap_centre_r = m_trap_half_length/std::cos(m_trap_polar_angle);

    double trap_centre_half_y = trap_centre_r*std::sin(m_trap_polar_angle)*std::sin(m_trap_azimuthal_angle) + m_trap_frontface_y/2.0;

    // distance in radial direction (from the interaction point) along the z direction of the trapezoid
    // Due to the coordinate system of the trapezoid, this is not in the centre of x and y
    double trap_rad_centre = m_calo_inner_r/std::cos(m_covered_theta) + m_back_shift + m_trap_half_length;

    // coordinates of the tower as if it was in global coordinates (was easier to visualise this way during development)
    double stave_x = std::cos(m_covered_theta)*trap_rad_centre - std::sin(m_covered_theta)*trap_centre_half_y;
    double stave_z = std::sin(m_covered_theta)*trap_rad_centre + std::cos(m_covered_theta)*trap_centre_half_y;

    // coordinates of the tower converted to the stave coordinate system
    double tower_x = 0;
    double tower_y = stave_z;
    double tower_z = stave_x-(m_calo_inner_r+m_stave_half_length);

    m_tower_position = dd4hep::Position(tower_x, tower_y, tower_z);
}


// Function to construct the trapezoidal supoprt structure for the tower in which fibres are placed
void DDDRCaloTubes::DRTubesconstructor::construct_tower_trapezoid(Volume& trap_volume)
{

        
        // polar coordinate conversion
        double delta_y = (m_trap_backface_y - m_trap_frontface_y)/2.0;
        double delta_z = 2.0*m_trap_half_length;
        m_trap_polar_angle = std::acos(delta_z/std::sqrt(delta_y*delta_y + delta_z*delta_z));
        m_trap_azimuthal_angle = 90.0*deg;

        Trap trap_solid("trap_solid", m_trap_half_length, m_trap_polar_angle, m_trap_azimuthal_angle, 
                                      m_trap_frontface_y/2.0, m_trap_frontface_rightangleedge_x/2.0, m_trap_frontface_thetaangleedge_x/2.0, 0.,
                                      m_trap_backface_y/2.0,  m_trap_backface_rightangleedge_x/2.0,  m_trap_backface_thetaangleedge_x/2.0,  0.);


        // Air volume in which fibres are placed
        double delta_y_air = (m_tower_backface_y - m_tower_frontface_y)/2.0;
        double delta_z_air = 2.0*m_tower_half_length;
        m_tower_polar_angle = std::acos(delta_z_air/std::sqrt(delta_y_air*delta_y_air + delta_z_air*delta_z_air));
        m_tower_azimuthal_angle = 90.0*deg;

        Trap tower_air_solid("tower_solid", m_tower_half_length, m_tower_polar_angle, m_tower_azimuthal_angle, 
                                      m_tower_frontface_y/2.0, m_tower_frontface_rightangleedge_x/2.0, m_tower_frontface_thetaangleedge_x/2.0, 0.,
                                      m_tower_backface_y/2.0,  m_tower_backface_rightangleedge_x/2.0,  m_tower_backface_thetaangleedge_x/2.0,  0.);

        // Position of the air volume in the trapezoid
        // y coordinate depends on wall thickness at the sides and how much the one wall is rotated in theta
        // z coordinates depedns on how thick the front and back walls are                        
        Position tower_air_pos = Position(0,
                                         (1.0-1.0/std::cos(m_tower_theta))*m_trap_wall_thickness_sides,
                                         (m_trap_wall_thickness_front-m_trap_wall_thickness_back)/2.0/* +10*nm */);


        // Subtraction solid used sometimes for easier visualisation. NOT TO BE USED IN FINAL GEOMETRY
        // SubtractionSolid solid = SubtractionSolid("trap_final", trap_solid, tower_air_solid, tower_air_pos);
        Volume tower_air_volume("tower_air_volume", tower_air_solid, m_air);
        tower_air_volume.setVisAttributes(*m_description, m_air_visString);

        trap_volume.setSolid(trap_solid);
        trap_volume.setVisAttributes(*m_description, m_trap_visString);

        PlacedVolume tower_air_placed = trap_volume.placeVolume(tower_air_volume, tower_air_pos);
        tower_air_placed.addPhysVolID("air", 1);

        // Place all the tubes inside the tower
        this->assemble_tower(tower_air_volume);
    
}


void DDDRCaloTubes::DRTubesconstructor::construct_tower(Volume& trap_volume)
{
    // For each placed tower, recalculate the parameters
    this->calculate_theta_parameters();
    // and construct the tower from this
    this->construct_tower_trapezoid(trap_volume);
}


// Placement of the tower in the stave volume
void DDDRCaloTubes::DRTubesconstructor::place_tower(Volume& stave_volume,
                 Volume& tower_volume,
                 unsigned int tower)
{
    
    double tower_x = m_tower_position.x();
    double tower_y = m_tower_position.y();
    double tower_z = m_tower_position.z();

    // Forward barrel region
    RotationX rot_fwd = RotationX(-m_covered_theta);
    Transform3D tower_fwd_tr(rot_fwd, Position(tower_x, tower_y, tower_z));
    PlacedVolume tower_fwd_placed = stave_volume.placeVolume(tower_volume, tower, tower_fwd_tr);
    tower_fwd_placed.addPhysVolID("tower", tower);

    // Backward barrel region
    Position m_tower_bwd_pos = Position(tower_x, -tower_y, tower_z);
    RotationZ rot_first_bwd = RotationZ(180*deg);
    RotationX rot_second_bwd = RotationX(m_covered_theta);
    Transform3D tower_bwd_tr(rot_second_bwd*rot_first_bwd, m_tower_bwd_pos);
    PlacedVolume tower_bwd_placed = stave_volume.placeVolume(tower_volume, -tower, tower_bwd_tr);
    tower_bwd_placed.addPhysVolID("tower", -tower);

}


void DDDRCaloTubes::DRTubesconstructor::construct_calorimeter(Volume& calorimeter_volume)
{
    // Parameters for stave contruction. Shape is a trapezoid over the full barrel region (forward and backward)
    double dy1 = m_calo_inner_half_z;
    double dy2 = m_calo_inner_half_z+2*m_stave_half_length;
    double dx1 = m_calo_inner_r*m_tower_tan_half_phi;
    double dx2 = m_calo_outer_r*std::sin(m_tower_half_phi);
    Trap stave_solid("stave_solid", m_stave_half_length, 0., 0., 
                     dy1, dx1, dx1, 0.,
                     dy2, dx2, dx2, 0.);
    Volume stave_volume("stave_volume", stave_solid, m_air);
    stave_volume.setVisAttributes(*m_description, "stave_vis");
    RotationZ rot_first = RotationZ(90*deg);
    RotationY rot_second = RotationY(90*deg);
    short int tower = 1;
    // Place towers in theta direection into the stave as long we are in the barrel region
    while (m_covered_theta<m_barrel_endcap_angle)
    {
        std::cout << "tower = " << tower << std::endl;
        Volume trap_volume("tower");
        trap_volume.setMaterial(m_trap_material);
        this->construct_tower(trap_volume);

        this->calculate_tower_position();
        /* if (tower==8)  */this->place_tower(stave_volume, trap_volume, tower);
        this->increase_covered_theta(m_tower_theta);
        
        // if (tower >= 1) break;
        tower++;
    }

    double phi = 0*deg;
    // Variable used to calculate stave position
    double centre_stave_vol = m_calo_inner_r + m_stave_half_length;


    // Placing of the staves
    for (unsigned int stave=0; stave<m_num_phi_towers; stave++, phi+=m_tower_phi)
    {
        RotationZ rot_third = RotationZ(phi);
        // stave position in the calorimeter volume
        double stave_x = centre_stave_vol*std::cos(phi);
        double stave_y = centre_stave_vol*std::sin(phi);
        Transform3D stave_tr(rot_third*rot_second*rot_first, Position(stave_x,stave_y,0));
        PlacedVolume stave_placed = calorimeter_volume.placeVolume(stave_volume, stave, stave_tr);
        stave_placed.addPhysVolID("stave", stave);
    }

    //Print length of tube map m_cher_tube_volume_map and m_scin_tube_volume_map
    std::cout << "Length of C map = " << m_cher_tube_volume_map.size() << std::endl;
    std::cout << "Length of S map = " << m_scin_tube_volume_map.size() << std::endl;

}
