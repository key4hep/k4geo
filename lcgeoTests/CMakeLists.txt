
#--------------------------------------------------

configure_file( ${DD4hep_ROOT}/cmake/run_test_package.sh run_test_${PackageName}.sh @ONLY)
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/run_test_${PackageName}.sh
  DESTINATION bin )

#--------------------------------------------------
#ADD_TEST( t_init source "${CMAKE_CURRENT_WORK_DIR}/thisdd4hep.sh" )

#--------------------------------------------------
# test(s) for ILD
SET( test_name "test_ILD_l5_v02" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v02/ILD_l5_v02.xml --runType=batch -G -N=1 --outputFile=testILD_l5_v02.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_ILD_s5_v02" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_s5_v02/ILD_s5_v02.xml --runType=batch -G -N=1 --outputFile=testILD_s5_v02.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_ILD_l5_v09" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v09/ILD_l5_v09.xml --runType=batch -G -N=1 --outputFile=testILD_l5_v09.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

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

SET( det_name "CLD_o2_v05" )
ADD_TEST( t_test_${det_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/${det_name}/${det_name}.xml --runType=batch -G -N=1 --outputFile=test${det_name}.slcio )
SET_TESTS_PROPERTIES( t_test_${det_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_steeringFile" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../example/steeringFile.py --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v04.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )


#--------------------------------------------------
# test for ARC o1 v01
SET( test_name "test_ARC_o1_v01_run" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	python3 ${CMAKE_CURRENT_SOURCE_DIR}/../example/arcfullsim.py --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/CLD_o3_v01/CLD_o3_v01.xml --runType batch  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  " Exception;EXCEPTION;ERROR;Error" )
set_tests_properties( t_${test_name} PROPERTIES TIMEOUT 120) 

SET( test_name "test_ARC_o1_v01_overlap" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	ddsim --compactFile=${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/CLD_o3_v01/CLD_o3_v01.xml --runType run --part.userParticleHandler= --macroFile=${CMAKE_CURRENT_SOURCE_DIR}/../utils/overlap.mac  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )
set_tests_properties( t_${test_name} PROPERTIES TIMEOUT 240) 

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
