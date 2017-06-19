//==========================================================================
//  Surface installer plugin for generic sliced detector drivers
//--------------------------------------------------------------------------
// DEPRECATED! TO BE REMOVED! USE THE DD4HEP IMPLEMENTATION INSTEAD!
// Author     : N. Nikiforou, adapted from DD4hep/SiTrackerBarrel_surfaces.cpp
//              by M. Frank
//==========================================================================
// Assumes boxes are stacked along one of the primary axes, x,y or z 
// The normal vector (n) must be along one of the axes (i.e. only one non-zero 
// component). The inner/outer thicknesses are calculated according to n.
// Note: Assumes module and component volumes are boxes. For Trapezoids,
// a fitting box is built around the trapezoid which means dx1=dx2=dx1 and
// dy1=dy2=dy. This is in principle fine, since we only access the thicknesses
// (DY in the TrackerEncapSurfacePlugin example) which is supposed to be the same
//==========================================================================

namespace { 
    struct UserData { 
        int dimension ;
        double uvector[3];
        double vvector[3];
        double nvector[3];
        double ovector[3];
        
    }; 
    
}

// Framework include files
#define SURFACEINSTALLER_DATA UserData
#define DD4HEP_USE_SURFACEINSTALL_HELPER GenericSurfaceInstallerPlugin
#include "DD4hep/SurfaceInstaller.h"

// #include <stdio.h>
// #include <string.h>

namespace{
    template <> void Installer<UserData>::handle_arguments(int argc, char** argv)   {

        //Initialize defaults to zero
        data.uvector[0]=0.;
        data.uvector[1]=0.;
        data.uvector[2]=0.;
        data.vvector[0]=0.;
        data.vvector[1]=0.;
        data.vvector[2]=0.;
        data.nvector[0]=0.;
        data.nvector[1]=0.;
        data.nvector[2]=0.;
        data.ovector[0]=0.;
        data.ovector[1]=0.;
        data.ovector[2]=0.;
        std::cout<<"WARNING! The lcgeo implementation of the GenericSurfaceInstallerPlugin is DEPRECATED and will be removed! Use the DD4hep implementation instead (DD4hep_GenericSurfaceInstallerPlugin)."<<std::endl;
        for(int i=0; i<argc; ++i)  {
            double value = -1;
            char* ptr = ::strchr(argv[i],'=');
            if ( ptr )  {
                std::string name( argv[i] , ptr ) ;
                value = dd4hep::_toDouble(++ptr);
                std::cout << "GenericSurfaceInstallerPlugin: argument[" << i << "] = " << name 
                << " = " << value << std::endl;
                if( name=="dimension" ) data.dimension = value ; 
                if( name=="u_x" ) data.uvector[0]=value ; 
                if( name=="u_y" ) data.uvector[1]=value ; 
                if( name=="u_z" ) data.uvector[2]=value ; 
                if( name=="v_x" ) data.vvector[0]=value ; 
                if( name=="v_y" ) data.vvector[1]=value ; 
                if( name=="v_z" ) data.vvector[2]=value ; 
                if( name=="n_x" ) data.nvector[0]=value ; 
                if( name=="n_y" ) data.nvector[1]=value ; 
                if( name=="n_z" ) data.nvector[2]=value ; 
                if( name=="o_x" ) data.ovector[0]=value ; 
                if( name=="o_y" ) data.ovector[1]=value ; 
                if( name=="o_z" ) data.ovector[2]=value ; 
            }
        }
    
        std::cout <<"GenericSurfaceInstallerPlugin: vectors: ";
        std::cout <<"u( "<<data.uvector[0]<<" , "<<data.uvector[1]<<" , "<<data.uvector[2]<<") ";
        std::cout <<"v( "<<data.vvector[0]<<" , "<<data.vvector[1]<<" , "<<data.vvector[2]<<") ";
        std::cout <<"n( "<<data.nvector[0]<<" , "<<data.nvector[1]<<" , "<<data.nvector[2]<<") ";
        std::cout <<"o( "<<data.ovector[0]<<" , "<<data.ovector[1]<<" , "<<data.ovector[2]<<") "<<std::endl;
    
    }  
    
    /// Install measurement surfaces
    template <typename UserData> 
      void Installer<UserData>::install(dd4hep::DetElement component, dd4hep::PlacedVolume pv)   {
        
        
        dd4hep::Volume comp_vol = pv.volume();
        if ( comp_vol.isSensitive() )  {  
            dd4hep::Volume mod_vol  = parentVolume(component);
            //FIXME: WHAT IF TRAPEZOID? Should work if trapezoid since it will fit minimal box and dy1=dy2=dy
            dd4hep::Box mod_shape(mod_vol.solid()), comp_shape(comp_vol.solid());
            
            if ( !comp_shape.isValid() || !mod_shape.isValid() )   {
                invalidInstaller("Components and/or modules are not boxes -- invalid shapes");

            }else if ( !handleUsingCache(component,comp_vol) )  {
                const double* trans = placementTranslation(component);
                double half_module_thickness = 0.;
                double sensitive_z_position  = 0.;

                if (data.nvector[0] !=0 && data.nvector[1] ==0 && data.nvector[2] ==0){
                    half_module_thickness = mod_shape->GetDX();
                    sensitive_z_position = data.nvector[0]>0 ? trans[0] : -trans[0];
                }else if (data.nvector[1] !=0 && data.nvector[0] ==0 && data.nvector[2] ==0){
                    half_module_thickness = mod_shape->GetDY();
                    sensitive_z_position = data.nvector[1]>0 ? trans[1] : -trans[1];

                }else if (data.nvector[2] !=0 && data.nvector[0] ==0 && data.nvector[1] ==0){
                    half_module_thickness = mod_shape->GetDZ();
                    sensitive_z_position = data.nvector[2]>0 ? trans[2] : -trans[2];

                } else {
                    throw std::runtime_error("**** GenericSurfaceInstallerPlugin: normal vector unsupported! It has to be "
                    "perpenidcular to one of the box sides, i.e. only one non-zero component.") ;
                }

                double inner_thickness = half_module_thickness + sensitive_z_position;
                double outer_thickness = half_module_thickness - sensitive_z_position;

                //Surface is placed at the center of the volume, no need to shift origin
                //Make sure u,v,n form a right-handed coordinate system, v along the final z
                Vector3D u(data.uvector[0],data.uvector[1],data.uvector[2]);
                Vector3D v(data.vvector[0],data.vvector[1],data.vvector[2]);
                Vector3D n(data.nvector[0],data.nvector[1],data.nvector[2]);
                Vector3D o(data.ovector[0],data.ovector[1],data.ovector[2]);
                Type type( Type::Sensitive ) ;
                
                if( data.dimension == 1 ) {
                    type.setProperty( Type::Measurement1D , true ) ;
                } else if( data.dimension != 2 ) {
                    throw std::runtime_error("**** GenericSurfaceInstallerPlugin: no or wrong "
                    "'dimension' argument given - has to be 1 or 2") ;
                }
                VolPlane surf(comp_vol, type, inner_thickness, outer_thickness, u, v, n, o);
                addSurface(component,surf);
            }
        }
    }
}// namespace
