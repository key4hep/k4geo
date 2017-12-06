
#--------------------------------------------------

# include_directories( 
#   )

SET(test_link_libraries DDCore DDRec DDSegmentation ${ROOT_LIBRARIES} Reflex )


configure_file( ${DD4hep_ROOT}/cmake/run_test_package.sh run_test_${PackageName}.sh @ONLY)
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/run_test_${PackageName}.sh
  DESTINATION bin )

#--------------------------------------------------
#ADD_TEST( t_init source "${CMAKE_CURRENT_WORK_DIR}/thisdd4hep.sh" )

#--------------------------------------------------
# test(s) for ILD
SET( test_name "test_ILD_l5_v02" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml --runType=batch -G -N=1 --outputFile=testILD_l5_v02.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_ILD_s5_v02" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/ILD/compact/ILD_s5_v02/ILD_s5_v02.xml --runType=batch -G -N=1 --outputFile=testILD_s5_v02.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

#--------------------------------------------------
# tests for SiD
SET( test_name "test_SiD_o2_v02" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/SiD/compact/SiD_o2_v02/SiD_o2_v02.xml --runType=batch -G -N=1 --outputFile=testSiD_o2_v02.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_SiD_o2_v03" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/SiD/compact/SiD_o2_v03/SiD_o2_v03.xml --runType=batch -G -N=1 --outputFile=testSiD_o2_v03.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

#--------------------------------------------------
# test for CLIC_o1_v01
SET( test_name "test_CLIC_o1_v01" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o1_v01/CLIC_o1_v01.xml --runType=batch -G -N=1 --outputFile=testCLIC.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o2_v04" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v04.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o3_v14" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
        ddsim --compactFile=${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o3_v14/CLIC_o3_v14.xml --runType=batch -G -N=1 --outputFile=testCLIC_o3_v14.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_FCCee_o1_v02" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
	ddsim --compactFile=${CMAKE_SOURCE_DIR}/FCCee/compact/FCCee_o1_v02/FCCee_o1_v02.xml --runType=batch -G -N=1 --outputFile=testFCCee_o1_v02.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )


SET( test_name "test_steeringFile" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --steeringFile=${CMAKE_SOURCE_DIR}/example/steeringFile.py --compactFile=${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v04.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

#--------------------------------------------------

ADD_EXECUTABLE( TestSensThickness src/TestSensThickness.cpp )
Target_Link_Libraries( TestSensThickness ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES} ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} )
INSTALL( TARGETS TestSensThickness DESTINATION bin )

ADD_EXECUTABLE( BeamCalZtest src/BeamCalZtest.cpp )
Target_Link_Libraries( BeamCalZtest ${DD4hep_LIBRARIES} ${ROOT_LIBRARIES} )
INSTALL( TARGETS BeamCalZtest DESTINATION bin )

ADD_TEST( t_SensThickness_Clic_o2_v4 "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
          ${CMAKE_INSTALL_PREFIX}/bin/TestSensThickness ${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml 300 50 )
ADD_TEST( t_SensThickness_CLIC_o3_v14 "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
          ${CMAKE_INSTALL_PREFIX}/bin/TestSensThickness ${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o3_v14/CLIC_o3_v14.xml 100 50 )
