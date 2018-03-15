def writeTopCompactXml( outfile, version='99', name="ILD", Large=True, Option=0, SolenoidMap=False, AntiDID=False, FwdFields=False, Energy=500 ):

    # defaults for option 0, 1
    # si ecal
    ecal_sl1=4
    ecal_sl2=10
    # ahcal
    hcal_sl=3
    if Option==2 or Option==4: # SDHCAL
        hcal_sl=1
    elif Option==3 or Option==4: # ScECAL
        ecal_sl1=3
        ecal_sl2=11
    elif Option>1:
        print 'ERROR: do not know about Option', Option
        return

    if Energy==250:
        enString='250GeV'
    elif Energy==500:
        enString='500GeV'
    else:
        print 'ERROR: do not know about Energy', Energy
        return

    outfile.write('<lccdd xmlns:compact="http://www.lcsim.org/schemas/compact/1.0"\n')
    outfile.write('       xmlns:xs="http://www.w3.org/2001/XMLSchema"\n')
    outfile.write('       xs:noNamespaceSchemaLocation="http://www.lcsim.org/schemas/compact/1.0/compact.xsd">\n')
    outfile.write('  <info name="'+name+'"\n')
    outfile.write('        title="ILD multi-technology model used for the optimisation"\n')
    outfile.write('        author="F. Gaede, S.Lu, D.Jeans"\n')
    outfile.write('        url="http://ilcsoft.desy.de"\n')
    outfile.write('        status="experimental"\n')
    outfile.write('        version="v'+version+'">\n')
    outfile.write('    <comment>ILD simulation models used for detector optimisation </comment>\n')
    outfile.write('  </info>\n')
    outfile.write('  <includes>\n')
    outfile.write('    <gdmlFile  ref="../ILD_common_v02/elements.xml"/>\n')
    outfile.write('    <gdmlFile  ref="../ILD_common_v02/materials.xml"/>\n')
    outfile.write('  </includes>\n')
    outfile.write('  <define>\n')
    if Large:
        outfile.write('    <include ref="../ILD_common_v02/top_defs_ILD_l5_v02.xml"/>\n')
    else:
        outfile.write('    <include ref="../ILD_common_v02/top_defs_ILD_s5_v02.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/top_defs_common_v02.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/basic_defs.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/envelope_defs.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/tube_defs.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/misc_defs.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/tracker_defs.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/fcal_defs.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/ecal_hybrid_defs.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/hcal_defs.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/yoke_defs.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/services_defs.xml"/>\n')
    outfile.write('    <include ref="${DD4hepINSTALL}/DDDetectors/compact/detector_types.xml"/>\n')
    outfile.write('    <include ref="../ILD_common_v02/limits.xml"/>\n')
    outfile.write('    <!-- Readout slice in ecal for reconstruction -->\n')
    outfile.write('    <constant name="Ecal_readout_segmentation_slice0" value="'+str(ecal_sl1)+'"/>\n')
    outfile.write('    <constant name="Ecal_readout_segmentation_slice1" value="'+str(ecal_sl2)+'"/>\n')
    outfile.write('    <!-- Readout slice in hcal for reconstruction -->\n')
    outfile.write('    <constant name="Hcal_readout_segmentation_slice" value="'+str(hcal_sl)+'"/>\n')
    outfile.write('  </define>\n')
    outfile.write('  <limits>\n')
    outfile.write('    <limitset name="cal_limits">\n')
    outfile.write('      <limit name="step_length_max" particles="*" value="cal_steplimit_val" unit="cal_steplimit_unit" />\n')
    outfile.write('    </limitset>\n')
    outfile.write('    <limitset name="TPC_limits">\n')
    outfile.write('      <limit name="step_length_max" particles="*" value="tpc_steplimit_val" unit="tpc_steplimit_unit" />\n')
    outfile.write('    </limitset>\n')
    outfile.write('    <limitset name="Tracker_limits">\n')
    outfile.write('      <limit name="step_length_max" particles="*" value="tracker_steplimit_val" unit="tracker_steplimit_unit" />\n')
    outfile.write('    </limitset>\n')
    outfile.write('  </limits>\n')
    outfile.write('  <include ref="../ILD_common_v02/display.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/Beampipe_o1_v01_01.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/vxd07.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/ftd_simple_staggered_02.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/sit_simple_pixel_sensors_01.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/tpc10_01.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/set_simple_planar_sensors_01.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/SEcal06_hybrid_Barrel.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/SEcal06_hybrid_Endcaps.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/SEcal05_siw_ECRing.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/SHcalSc04_Barrel_v04.xml"/>\n')
    
    if Large:
        outfile.write('  <include ref="../ILD_common_v02/SHcalSc04_Endcaps_v01_LARGE.xml"/>\n')
    else:
        outfile.write('  <include ref="../ILD_common_v02/SHcalSc04_Endcaps_v01_SMALL.xml"/>\n')

    outfile.write('  <include ref="../ILD_common_v02/SHcalSc04_EndcapRing_v01.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/Yoke05_Barrel.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/Yoke06_Endcaps.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/LumiCal.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/LHCal01.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/BeamCal08.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/coil03.xml"/>\n')
    outfile.write('  <include ref="../ILD_common_v02/SServices00.xml"/>\n')
    outfile.write('  <plugins>\n')
    outfile.write('    <plugin name="DD4hepVolumeManager"/>\n')
    outfile.write('    <plugin name="InstallSurfaceManager"/>\n')
    outfile.write('  </plugins>\n')

    if SolenoidMap:
        if Large:
            outfile.write('  <include ref="../ILD_common_v02/Field_Solenoid_Map_l_3.5T.xml"/>\n')
        else:
            outfile.write('  <include ref="../ILD_common_v02/Field_Solenoid_Map_s_4.0T.xml"/>\n')
    else:
        outfile.write('  <include ref="../ILD_common_v02/Field_Solenoid_Ideal.xml"/>\n')

    if AntiDID:
        if Large:
            outfile.write('  <include ref="../ILD_common_v02/Field_AntiDID_Map_l.xml"/>\n')
        else:
            outfile.write('  <include ref="../ILD_common_v02/Field_AntiDID_Map_s.xml"/>\n')

    if FwdFields:
        outfile.write('  <include ref="../ILD_common_v02/Field_FwdMagnets_Ideal_'+enString+'.xml"/>\n')

    outfile.write('</lccdd>\n')

def getVersionNumber(SolenoidMap, AntiDID, Energy):
    version='99'
    if not SolenoidMap:
        if not AntiDID:
            version='02'
        else:
            print 'ERROR; no model defined with ideal solnoid + antiDID!!'
            return
    else: # solenoid map
        if not AntiDID:
            if Energy==250:
                version='03'
            elif Energy==500:
                version='04'
            else:
                print 'ERROR; no model defined for energy', Energy
                return
        else:
            if Energy==250:
                version='05'
            elif Energy==500:
                version='06'
            else:
                print 'ERROR; no model defined for energy', Energy
                return
    return version



#-----------------------------------------------------

import os

prename='ILD' # ILD for final
topdir='../'
mainoutdirname=prename+'_sl5_v02/'
mainoutdir=topdir+mainoutdirname

if not os.path.exists(mainoutdir):
    os.makedirs(mainoutdir)

for Large in (True, False):
    for Option in (0,1,2,3,4):
        for SolenoidMap in (True, False):
            FwdFields=SolenoidMap # always use fwd fields if using solenoid map, don't if not

            ens=[250]
            if FwdFields:
                ens.append(500)

            antis=[False]
            if SolenoidMap:
                antis.append(True)

            for Energy in ens:
                for AntiDID in antis:
                    modelname=prename+'_'
                    if Large:
                        modelname=modelname+'l'
                    else:
                        modelname=modelname+'s'
                    modelname=modelname+'5_'
                    if Option>0:
                        modelname=modelname+'o'+str(Option)+'_'

                    version=getVersionNumber(SolenoidMap, AntiDID, Energy)
                    modelname=modelname+'v'+version
                    
                    outfile=open(mainoutdir+modelname+'.xml','w')
                    writeTopCompactXml( outfile, version=version, name=modelname, Large=Large, Option=Option, SolenoidMap=SolenoidMap, AntiDID=AntiDID, FwdFields=FwdFields, Energy=Energy )
                    outfile.close()

                    if not os.path.exists( topdir+modelname ):
                        print mainoutdirname, modelname
                        os.symlink( mainoutdirname, topdir+modelname )
