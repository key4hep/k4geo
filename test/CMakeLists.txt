
#--------------------------------------------------

configure_file( ${DD4hep_ROOT}/cmake/run_test_package.sh run_test_${PackageName}.sh @ONLY)
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/run_test_${PackageName}.sh
  DESTINATION bin )

macro(SET_TEST_ENV test_name)
  # Make sure to point our current installation for the compact files
  SET_TESTS_PROPERTIES( t_${test_name}
    PROPERTIES
      FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error"
      ENVIRONMENT
        K4GEO=${CMAKE_INSTALL_PREFIX}/share/k4geo
    )

endmacro()

#--------------------------------------------------
#ADD_TEST( t_init source "${CMAKE_CURRENT_WORK_DIR}/thisdd4hep.sh" )

#--------------------------------------------------
# test(s) for ILD
SET( test_name "test_ILD_l5_v02" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v02/ILD_l5_v02.xml --runType=batch -G -N=1 --outputFile=testILD_l5_v02.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_ILD_s5_v08" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_s5_v08/ILD_s5_v08.xml --runType=batch -G -N=1 --outputFile=testILD_s5_v08.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_ILD_l5_v09" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v09/ILD_l5_v09.xml --runType=batch -G -N=1 --outputFile=testILD_l5_v09.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_ILD_l5_v10" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v10/ILD_l5_v10.xml --runType=batch -G -N=1 --outputFile=testILD_l5_v10.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_ILD_l5_v11" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v11/ILD_l5_v11.xml --runType=batch -G -N=1 --outputFile=testILD_l5_v11.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_ILD_FCCee_v01" )
ADD_TEST( t_${test_name} sh -c "
 ${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ILD_FCCee/compact/ILD_FCCee_v01/ILD_FCCee_v01.xml --runType=batch -G -N=1 --outputFile=testILD_FCCee_v01.slcio 2>&1 |
  python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/filter_G4Tessellated_warnings.py")
SET_TEST_ENV(${test_name})

SET( test_name "test_ILD_FCCee_v02" )
ADD_TEST( t_${test_name} sh -c "
 ${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ILD_FCCee/compact/ILD_FCCee_v02/ILD_FCCee_v02.xml --runType=batch -G -N=1 --outputFile=testILD_FCCee_v02.slcio 2>&1 |
  python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/filter_G4Tessellated_warnings.py")
SET_TEST_ENV(${test_name})

#--------------------------------------------------
# tests for SiD
SET( test_name "test_SiD_o2_v03" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../SiD/compact/SiD_o2_v03/SiD_o2_v03.xml --runType=batch -G -N=1 --outputFile=testSiD_o2_v03.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_SiD_o2_v04" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../SiD/compact/SiD_o2_v04/SiD_o2_v04.xml --runType=batch -G -N=1 --outputFile=testSiD_o2_v04.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

#--------------------------------------------------
# test for CLIC_o1_v01
SET( test_name "test_CLIC_o1_v01" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o1_v01/CLIC_o1_v01.xml --runType=batch -G -N=1 --outputFile=testCLIC.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o2_v04" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v04.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o3_v15" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
        ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o3_v15/CLIC_o3_v15.xml --runType=batch -G -N=1 --outputFile=testCLIC_o3_v15.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_FCCee_o1_v05" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/FCCee_o1_v05/FCCee_o1_v05.xml --runType=batch -G -N=1 --outputFile=testFCCee_o1_v05.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( det_name "CLD_o2_v08" )
ADD_TEST( t_test_${det_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/${det_name}/${det_name}.xml --runType=batch -G -N=1 --outputFile=test${det_name}.slcio )
SET_TESTS_PROPERTIES( t_test_${det_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( det_name "CLD_o4_v05" )
ADD_TEST( t_test_${det_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/${det_name}/${det_name}.xml --runType=batch -G -N=1 --outputFile=test${det_name}.slcio )
SET_TESTS_PROPERTIES( t_test_${det_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_steeringFile" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../example/steeringFile.py --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v04.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

#--------------------------------------------------
# test for IDEA o1 v02
SET( test_name "test_IDEA_o1_v02" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
    ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/IDEA/compact/IDEA_o1_v02/IDEA_o1_v02.xml --runType=batch -G -N=1 --gun.distribution uniform --gun.particle geantino --random.seed 1988301045 --outputFile=testIDEA_o1_v02.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" TIMEOUT 360)

#--------------------------------------------------
# test for IDEA o1 v03
if(DCH_INFO_H_EXIST)
SET( test_name "test_IDEA_o1_v03" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml -G -N=1 --gun.distribution uniform --gun.particle geantino --random.seed 1988301045 --steeringFile ${CMAKE_CURRENT_SOURCE_DIR}/../example/SteeringFile_IDEA_o1_v03.py)
    SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  " Exception; EXCEPTION;ERROR;Error" TIMEOUT 900)
endif()

#--------------------------------------------------
# test for IDEA o1 v03, with DRC, with Geantinos
if(DCH_INFO_H_EXIST)
SET( test_name "test_IDEA_with_DRC_o1_v03" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/compact/IDEA_withDRC_o1_v03.xml --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../example/SteeringFile_IDEA_o1_v03.py -G -N=1 --gun.distribution uniform --gun.particle geantino --random.seed 1988301045 )
    SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  " Exception; EXCEPTION;ERROR;Error" TIMEOUT 900)
endif()

#--------------------------------------------------
# test for IDEA o1 v03, with DRC, with electrons
if(DCH_INFO_H_EXIST)
SET( test_name "test_IDEA_with_DRC_o1_v03" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
        ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/compact/IDEA_withDRC_o1_v03.xml --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../example/SteeringFile_IDEA_o1_v03.py -G --gun.distribution uniform --gun.particle e- --random.seed 1988301045 )
    SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  " Exception; EXCEPTION;ERROR;Error" TIMEOUT 900)
endif()

#--------------------------------------------------
# test for IDEA o2 v01
#if(DCH_INFO_H_EXIST)
#SET( test_name "test_IDEA_o2_v01" )
#ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
#    ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/IDEA/compact/IDEA_o2_v01/IDEA_o2_v01.xml --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../example/SteeringFile_IDEA_o2_v01.py -N 1 --random.seed 1988301045 --outputFile=testIDEA_o2_v01.root )
#    SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  " Exception; EXCEPTION;ERROR;Error" TIMEOUT 600)
#endif()

#--------------------------------------------------
# test for ALLEGRO o1 v02
SET( test_name "test_ALLEGRO_o1_v02" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ALLEGRO/compact/ALLEGRO_o1_v02/ALLEGRO_o1_v02.xml --runType=batch -G -N=1 --gun.distribution uniform --gun.particle geantino --outputFile=testALLEGRO_o1_v02.root )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

#--------------------------------------------------
# test for ALLEGRO o1 v03
if(DCH_INFO_H_EXIST)
  SET( test_name "test_ALLEGRO_o1_v03" )
  ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
    ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml --runType=batch -G -N=1 --gun.distribution uniform --gun.particle geantino --outputFile=testALLEGRO_o1_v03.root )
  SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  " Exception; EXCEPTION;ERROR;Error" TIMEOUT 600)
endif()
#--------------------------------------------------
# test for ALLEGRO o2 v01
SET( test_name "test_ALLEGRO_o2_v01" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ALLEGRO/compact/ALLEGRO_o2_v01/ALLEGRO_o2_v01.xml --runType=batch -G -N=1 --gun.distribution uniform --gun.particle geantino --outputFile=testALLEGRO_o2_v01.root )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

#--------------------------------------------------
# test for ARC o1 v01
SET( test_name "test_ARC_o1_v01_run" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	python3 ${CMAKE_CURRENT_SOURCE_DIR}/../example/arcfullsim.py --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/compact/ARC_standalone_o1_v01.xml --runType batch  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  " Exception;EXCEPTION;ERROR;Error" )
set_tests_properties( t_${test_name} PROPERTIES TIMEOUT 120) 

SET( test_name "test_ARC_o1_v01_overlap" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/compact/ARC_standalone_o1_v01.xml --runType run --part.userParticleHandler= --macroFile=${CMAKE_CURRENT_SOURCE_DIR}/../utils/overlap.mac  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )
set_tests_properties( t_${test_name} PROPERTIES TIMEOUT 360) 

#--------------------------------------------------
# test for DCH o1 v02
if(DCH_INFO_H_EXIST)
SET( test_name "test_DCH_o1_v02_overlap" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/compact/DCH_standalone_o1_v02.xml --runType run --part.userParticleHandler= --macroFile=${CMAKE_CURRENT_SOURCE_DIR}/../utils/overlap.mac  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )
set_tests_properties( t_${test_name} PROPERTIES TIMEOUT 400)
endif()

#--------------------------------------------------

ADD_EXECUTABLE( TestSensThickness src/TestSensThickness.cpp )
Target_Link_Libraries( TestSensThickness lcgeo )
INSTALL( TARGETS TestSensThickness DESTINATION bin )

ADD_EXECUTABLE( BeamCalZtest src/BeamCalZtest.cpp )
Target_Link_Libraries( BeamCalZtest lcgeo )
INSTALL( TARGETS BeamCalZtest DESTINATION bin )

ADD_TEST( t_SensThickness_Clic_o2_v4 "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
          ${CMAKE_INSTALL_PREFIX}/bin/TestSensThickness ${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml 300 50 )
ADD_TEST( t_SensThickness_CLIC_o3_v15 "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
          ${CMAKE_INSTALL_PREFIX}/bin/TestSensThickness ${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o3_v15/CLIC_o3_v15.xml 100 50 )

#--------------------------------------------------
# Tests for MuonCollider geometries
#--------------------------------------------------
SET( test_name "test_MAIA_v0_run" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml --runType=batch -G -N=1 --outputFile=test_MAIA_v0.edm4hep.root
)
SET_TEST_ENV(${test_name})

SET( test_name "test_MuColl_v1_run" )
ADD_TEST(t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../MuColl/MuColl/compact/MuColl_v1/MuColl_v1.xml --runType=batch -G -N=1 --outputFile=test_MuColl_v1.edm4hep.root
)
SET_TEST_ENV(${test_name})

#--------------------------------------------------
# check if files named the same contain the same in FCCee
ADD_TEST(
    NAME t_test_files_versions_FCCee
    COMMAND "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
            "${PROJECT_SOURCE_DIR}/utils/compareIdenticalFiles.sh"
            "${PROJECT_SOURCE_DIR}/FCCee"
            "${PROJECT_SOURCE_DIR}/utils/IdenticalFiles_ToBeIgnored.txt"
)
