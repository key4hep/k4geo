
#--------------------------------------------------

# include_directories( 
#   )

SET(test_link_libraries DDCore DDRec DDSegmentation ${ROOT_LIBRARIES} Reflex )


configure_file( ${DD4hep_ROOT}/cmake/run_test_package.sh run_test_${PackageName}.sh @ONLY)
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/run_test_${PackageName}.sh
  DESTINATION bin )

#--------------------------------------------------
#ADD_TEST( t_init source "${CMAKE_CURRENT_WORK_DIR}/thisdd4hep.sh" )

SET( test_name "test_ILD_o1_v05" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/ILD/compact/ILD_o1_v05/ILD_o1_v05.xml --runType=batch -G -N=1 --outputFile=testILD.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

#--------------------------------------------------
# test(s) for SiD
SET( test_name "test_SiD_o1_v01" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/SiD/compact/SiD_o1_v01/SiD_o1_v01.xml --runType=batch -G -N=1 --outputFile=testSiD_o1_v01.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

#--------------------------------------------------
# test for CLIC_o1_v01
SET( test_name "test_CLIC_o1_v01" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o1_v01/CLIC_o1_v01.xml --runType=batch -G -N=1 --outputFile=testCLIC.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o2_v01" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o2_v01/CLIC_o2_v01.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o2_v02" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o2_v02/CLIC_o2_v02.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v02.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o2_v03" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o2_v03/CLIC_o2_v03.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v03.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o2_v04" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v04.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o2_v05" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o2_v05/CLIC_o2_v05.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v05.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_steeringFile" )
ADD_TEST( t_${test_name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
  ddsim --steeringFile=${CMAKE_SOURCE_DIR}/example/steeringFile.py --compactFile=${CMAKE_SOURCE_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v04.slcio )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

#--------------------------------------------------

