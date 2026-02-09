
#--------------------------------------------------
# test(s) for ILD
SET( test_name "test_ILD_l5_v02" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v02/ILD_l5_v02.xml --runType=batch -G -N=1 --outputFile=testILD_l5_v02.slcio )

SET( test_name "test_ILD_s5_v08" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_s5_v08/ILD_s5_v08.xml --runType=batch -G -N=1 --outputFile=testSILD_s5_v08.slcio )

SET( test_name "test_ILD_l5_v09" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v09/ILD_l5_v09.xml --runType=batch -G -N=1 --outputFile=testILD_l5_v09.slcio )

SET( test_name "test_ILD_l5_v10" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v10/ILD_l5_v10.xml --runType=batch -G -N=1 --outputFile=testILD_l5_v10.slcio )

SET( test_name "test_ILD_l5_v11" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v11/ILD_l5_v11.xml --runType=batch -G -N=1 --outputFile=testILD_l5_v11.slcio )

if(INSTALL_BEAMPIPE_STL_FILES)
  SET( test_name "test_ILD_FCCee_v01" )
  add_test(NAME ${test_name} COMMAND sh -c "
    ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ILD_FCCee/compact/ILD_FCCee_v01/ILD_FCCee_v01.xml --runType=batch -G -N=1 --outputFile=testILD_FCCee_v01.slcio 2>&1 |
    python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/filter_G4Tessellated_warnings.py")

  SET( test_name "test_ILD_FCCee_v02" )
  add_test(NAME ${test_name} COMMAND sh -c "
    ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ILD_FCCee/compact/ILD_FCCee_v02/ILD_FCCee_v02.xml --runType=batch -G -N=1 --outputFile=testILD_FCCee_v02.slcio 2>&1 |
    python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/filter_G4Tessellated_warnings.py")
endif()

#--------------------------------------------------
# tests for SiD
SET( test_name "test_SiD_o2_v03" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../SiD/compact/SiD_o2_v03/SiD_o2_v03.xml --runType=batch -G -N=1 --outputFile=testSiD_o2_v03.slcio )
SET( test_name "test_SiD_o2_v04" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../SiD/compact/SiD_o2_v04/SiD_o2_v04.xml --runType=batch -G -N=1 --outputFile=testSiD_o2_v04.slcio )

#--------------------------------------------------
# test for CLIC_o1_v01
SET( test_name "test_CLIC_o1_v01" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o1_v01/CLIC_o1_v01.xml --runType=batch -G -N=1 --outputFile=testCLIC.slcio )

SET( test_name "test_CLIC_o2_v04" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v04.slcio )

SET( test_name "test_CLIC_o3_v15" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o3_v15/CLIC_o3_v15.xml --runType=batch -G -N=1 --outputFile=testCLIC_o3_v15.slcio )

SET( test_name "test_FCCee_o1_v05" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/FCCee_o1_v05/FCCee_o1_v05.xml --runType=batch -G -N=1 --outputFile=testFCCee_o1_v05.slcio )

SET( det_name "CLD_o2_v08" )
add_test(NAME t_test_${det_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/${det_name}/${det_name}.xml --runType=batch -G -N=1 --outputFile=test${det_name}.slcio )

SET( det_name "CLD_o4_v05" )
add_test(NAME t_test_${det_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/${det_name}/${det_name}.xml --runType=batch -G -N=1 --outputFile=test${det_name}.slcio )

SET( test_name "test_steeringFile" )
add_test(NAME ${test_name} COMMAND ddsim --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../example/steeringFile.py --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v04.slcio )

#--------------------------------------------------
# test for IDEA o1 v02
SET( test_name "test_IDEA_o1_v02" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/IDEA/compact/IDEA_o1_v02/IDEA_o1_v02.xml --runType=batch -G -N=51000 --gun.distribution uniform --gun.particle geantino --random.seed 1988301045 --outputFile=testIDEA_o1_v02.slcio )
SET_TESTS_PROPERTIES( ${test_name} PROPERTIES TIMEOUT 360)

#--------------------------------------------------
# test for IDEA o1 v03
if(DCH_INFO_H_EXIST)
SET( test_name "test_IDEA_o1_v03" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml -G -N=1 --gun.distribution uniform --gun.particle geantino --random.seed 1988301045 --steeringFile ${CMAKE_CURRENT_SOURCE_DIR}/../example/SteeringFile_IDEA_o1_v03.py)
SET_TESTS_PROPERTIES( ${test_name} PROPERTIES TIMEOUT 900)
endif()

# test for IDEA o1 v03, with DRC, with electrons
if(DCH_INFO_H_EXIST)
SET( test_name "test_IDEA_with_DRC_o1_v03" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/compact/IDEA_withDRC_o1_v03.xml --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../example/SteeringFile_IDEA_o1_v03.py -G --gun.distribution uniform --gun.particle e- --random.seed 1988301045 )
SET_TESTS_PROPERTIES( ${test_name} PROPERTIES TIMEOUT 900)
endif()

#--------------------------------------------------
# test for IDEA o1 v04
if(DCH_INFO_H_EXIST)
SET( test_name "test_IDEA_o1_v04" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/IDEA/compact/IDEA_o1_v04/IDEA_o1_v04.xml -G -N=1 --gun.distribution uniform --gun.particle geantino --random.seed 1988301045 --steeringFile ${CMAKE_CURRENT_SOURCE_DIR}/../example/SteeringFile_IDEA_o1_v03.py)
SET_TESTS_PROPERTIES( ${test_name} PROPERTIES TIMEOUT 900)
endif()

#--------------------------------------------------
# test for IDEA o2 v01
#if(DCH_INFO_H_EXIST)
#SET( test_name "test_IDEA_o2_v01" )
#add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/IDEA/compact/IDEA_o2_v01/IDEA_o2_v01.xml --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../example/SteeringFile_IDEA_o2_v01.py -N 1 --random.seed 1988301045 --outputFile=testIDEA_o2_v01.root )
#    SET_TESTS_PROPERTIES( ${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  " Exception; EXCEPTION;ERROR;Error" TIMEOUT 600)
#endif()

#--------------------------------------------------
# test for ALLEGRO o1 v02
SET( test_name "test_ALLEGRO_o1_v02" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ALLEGRO/compact/ALLEGRO_o1_v02/ALLEGRO_o1_v02.xml --runType=batch -G -N=1 --gun.distribution uniform --gun.particle geantino --outputFile=testALLEGRO_o1_v02.root )

#--------------------------------------------------
# test for ALLEGRO o1 v03
if(DCH_INFO_H_EXIST)
  SET( test_name "test_ALLEGRO_o1_v03" )
  add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml --runType=batch -G -N=1 --gun.distribution uniform --gun.particle geantino --outputFile=testALLEGRO_o1_v03.root )
  SET_TESTS_PROPERTIES( ${test_name} PROPERTIES TIMEOUT 600)
endif()
#--------------------------------------------------
# test for ALLEGRO o2 v01
SET( test_name "test_ALLEGRO_o2_v01" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ALLEGRO/compact/ALLEGRO_o2_v01/ALLEGRO_o2_v01.xml --runType=batch -G -N=1 --gun.distribution uniform --gun.particle geantino --outputFile=testALLEGRO_o2_v01.root )

#--------------------------------------------------
# test for ARC o1 v01
SET( test_name "test_ARC_o1_v01_run" )
add_test(NAME ${test_name} COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../example/arcfullsim.py --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/compact/ARC_standalone_o1_v01.xml --runType batch  )
SET_TESTS_PROPERTIES( ${test_name} PROPERTIES TIMEOUT 120)

SET( test_name "test_ARC_o1_v01_overlap" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/compact/ARC_standalone_o1_v01.xml --runType run --part.userParticleHandler= --macroFile=${CMAKE_CURRENT_SOURCE_DIR}/../utils/overlap.mac  )
SET_TESTS_PROPERTIES( ${test_name} PROPERTIES TIMEOUT 240)

#--------------------------------------------------
# test for DCH o1 v02
if(DCH_INFO_H_EXIST)
SET( test_name "test_DCH_o1_v02_overlap" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/compact/DCH_standalone_o1_v02.xml --runType run --part.userParticleHandler= --macroFile=${CMAKE_CURRENT_SOURCE_DIR}/../utils/overlap.mac  )
set_tests_properties( ${test_name} PROPERTIES TIMEOUT 400)
endif()

#--------------------------------------------------

ADD_EXECUTABLE( TestSensThickness src/TestSensThickness.cpp )
Target_Link_Libraries( TestSensThickness lcgeo )
INSTALL( TARGETS TestSensThickness DESTINATION bin )

ADD_EXECUTABLE( BeamCalZtest src/BeamCalZtest.cpp )
Target_Link_Libraries( BeamCalZtest lcgeo )
INSTALL( TARGETS BeamCalZtest DESTINATION bin )

add_test(NAME t_SensThickness_Clic_o2_v4 COMMAND ${PROJECT_BINARY_DIR}/bin/TestSensThickness ${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml 300 50 )
add_test(NAME t_SensThickness_Clic_o3_v15 COMMAND ${PROJECT_BINARY_DIR}/bin/TestSensThickness ${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o3_v15/CLIC_o3_v15.xml 100 50 )

#--------------------------------------------------
# Tests for MuonCollider geometries
#--------------------------------------------------
SET( test_name "test_MAIA_v0_run" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml --runType=batch -G -N=1 --outputFile=test_MAIA_v0.edm4hep.root)

SET( test_name "test_MuColl_v1_run" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../MuColl/MuColl/compact/MuColl_v1/MuColl_v1.xml --runType=batch -G -N=1 --outputFile=test_MuColl_v1.edm4hep.root)

SET( test_name "test_MuSIC_v2_run" )
add_test(NAME ${test_name} COMMAND ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../MuColl/MuSIC/compact/MuSIC_v2/MuSIC_v2.xml --runType=batch -G -N=1 --outputFile=test_MuSIC_v2.edm4hep.root)

#--------------------------------------------------
# check if files named the same contain the same in FCCee
add_test(NAME t_test_files_versions_FCCee
    COMMAND "${PROJECT_SOURCE_DIR}/utils/compareIdenticalFiles.sh"
            "${PROJECT_SOURCE_DIR}/FCCee"
            "${PROJECT_SOURCE_DIR}/utils/IdenticalFiles_ToBeIgnored.txt"
)


set(test_environment "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}")
get_property(test_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)
set_tests_properties(${test_names} PROPERTIES ENVIRONMENT "${test_environment}"
                                              FAIL_REGULAR_EXPRESSION " Exception; EXCEPTION;ERROR;Error")
