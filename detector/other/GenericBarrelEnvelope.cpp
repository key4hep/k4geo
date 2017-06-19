#include "DD4hep/DetFactoryHelper.h"


using namespace std;

using dd4hep::Ref_t;
using dd4hep::Detector;
using dd4hep::SensitiveDetector;
using dd4hep::PolyhedraRegular;
using dd4hep::SubtractionSolid;
using dd4hep::Tube;
using dd4hep::Position;
using dd4hep::RotationZYX;
using dd4hep::DetElement;
using dd4hep::Volume;
using dd4hep::Material;
using dd4hep::Transform3D;
using dd4hep::PlacedVolume;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector /*sens*/) {
    
    xml_det_t Detector = e;
    
    int DetectorID = Detector.id();
    string DetectorName = Detector.nameStr();
    
    xml_comp_t Dimensions = Detector.dimensions();
    
    
    //Get the dimensions of the envelope
    int nsides_inner = Dimensions.nsides_inner();
    int nsides_outer = Dimensions.nsides_outer();
    
    double inner_r = Dimensions.rmin();
    double outer_r = Dimensions.rmax();
    double half_z = Dimensions.zhalf();
    double phi_offset = Dimensions.phi0_offset();

    Volume envelope;
    Material envelope_material = theDetector.material(Dimensions.materialStr());

    if(nsides_inner!=0 && nsides_outer!=0){
        if( (nsides_outer-nsides_inner)==0 && phi_offset==0){
            PolyhedraRegular polyhedron(nsides_outer, inner_r, outer_r, 2*half_z);
            envelope = Volume(DetectorName+"_envelope",polyhedron,envelope_material);
        }else{
            PolyhedraRegular polyhedron_inner (nsides_inner, 0, inner_r, 2.1*half_z); //Here the inner polyhedron hald z is increased for better substraction!
            PolyhedraRegular polyhedron_outer (nsides_outer, phi_offset, 0, outer_r, 2*half_z);
            SubtractionSolid polyhedron(polyhedron_outer,polyhedron_inner);
            envelope = Volume(DetectorName+"_envelope",polyhedron,envelope_material);
        }
    }else{
        if(nsides_inner==0 && nsides_outer==0){
            Tube tube (inner_r,outer_r, half_z);
            envelope = Volume(DetectorName+"_envelope",tube,envelope_material);
        }
        else if(nsides_inner==0 && nsides_outer!=0){
            Tube tube_inner (0, inner_r, 1.1*half_z); //Here the inner polyhedron hald z is increased for better substraction!
            PolyhedraRegular polyhedron_outer (nsides_outer, phi_offset, 0, outer_r, 2*half_z);
            SubtractionSolid polyhedron(polyhedron_outer,tube_inner);
            envelope = Volume(DetectorName+"_envelope",polyhedron,envelope_material);
        }else{
            PolyhedraRegular polyhedron_inner (nsides_inner, 0, inner_r, 2.1*half_z);
            Tube tube_outer (0, outer_r, half_z);
            SubtractionSolid polyhedron(tube_outer,polyhedron_inner);
            envelope = Volume(DetectorName+"_envelope",polyhedron,envelope_material);
        }
    }

    
    xml_comp_t envelope_position = Detector.position();
    xml_comp_t envelope_rotation = Detector.rotation();

    Position position(envelope_position.x(),envelope_position.y(),envelope_position.z());
    RotationZYX rotation(envelope_rotation.z(),envelope_rotation.y(),envelope_rotation.x());

    Transform3D tr(rotation,position);

    
    DetElement sdet (DetectorName,DetectorID);
    Volume motherVol = theDetector.pickMotherVolume(sdet);

    PlacedVolume  env_phv = motherVol.placeVolume(envelope,tr);
    sdet.setPlacement( env_phv );

    envelope.setAttributes(theDetector,Detector.regionStr(),Detector.limitsStr(),Detector.visStr());

    return sdet;
}

DECLARE_DETELEMENT(GenericBarrelEnvelope,create_detector)
