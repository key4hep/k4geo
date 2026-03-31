
# Helper function to add a ddsim test with a simple particle gun (-G flag).
# The test name and output file are derived from the compact file stem.
# Usage:
#   add_ddsim_gun_test(<compact_file>
#     [NAME <name>]                # override derived test name (test_<stem>)
#     [N_EVENTS <n>]               # number of events (default: 1)
#     [FILTER_SCRIPT <path>]       # pipe output through: 2>&1 | python3 <path>
#     [EXTRA_ARGS <arg1> ...]      # additional ddsim arguments
#     [PROPERTIES <key> <val> ...] # passed verbatim to set_tests_properties
#   )
function(add_ddsim_gun_test compact_file)
  cmake_parse_arguments(ARG
    ""
    "NAME;N_EVENTS;FILTER_SCRIPT"
    "EXTRA_ARGS;PROPERTIES"
    ${ARGN})

  get_filename_component(stem ${compact_file} NAME_WE)

  if(ARG_NAME)
    set(test_name ${ARG_NAME})
  else()
    set(test_name test_${stem})
  endif()

  if(NOT ARG_N_EVENTS)
    set(ARG_N_EVENTS 1)
  endif()

  set(cmd ddsim --compactFile=${compact_file} --runType=batch -G -N=${ARG_N_EVENTS}
    --outputFile=test${stem}.edm4hep.root)

  if(ARG_EXTRA_ARGS)
    list(APPEND cmd ${ARG_EXTRA_ARGS})
  endif()

  if(ARG_FILTER_SCRIPT)
    string(JOIN " " cmd_str ${cmd})
    add_test(NAME ${test_name} COMMAND sh -c "${cmd_str} 2>&1 | python3 ${ARG_FILTER_SCRIPT}")
  else()
    add_test(NAME ${test_name} COMMAND ${cmd})
  endif()

  if(ARG_PROPERTIES)
    set_tests_properties(${test_name} PROPERTIES ${ARG_PROPERTIES})
  endif()
endfunction()

#--------------------------------------------------
# tests for ILD
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v02/ILD_l5_v02.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_s5_v08/ILD_s5_v08.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v09/ILD_l5_v09.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v10/ILD_l5_v10.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_l5_v11/ILD_l5_v11.xml)

if(INSTALL_BEAMPIPE_STL_FILES)
  add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ILD_FCCee/compact/ILD_FCCee_v01/ILD_FCCee_v01.xml
    FILTER_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/scripts/filter_G4Tessellated_warnings.py)

  add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ILD_FCCee/compact/ILD_FCCee_v02/ILD_FCCee_v02.xml
    FILTER_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/scripts/filter_G4Tessellated_warnings.py)
endif()

#--------------------------------------------------
# tests for SiD
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../SiD/compact/SiD_o2_v03/SiD_o2_v03.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../SiD/compact/SiD_o2_v04/SiD_o2_v04.xml)

#--------------------------------------------------
# tests for CLIC
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o1_v01/CLIC_o1_v01.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o3_v15/CLIC_o3_v15.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/FCCee_o1_v05/FCCee_o1_v05.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/CLD_o2_v08/CLD_o2_v08.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/CLD_o4_v05/CLD_o4_v05.xml)

add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml
  NAME test_steeringFile
  EXTRA_ARGS --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../example/steeringFile.py)

#--------------------------------------------------
# test for IDEA o1 v02
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/IDEA/compact/IDEA_o1_v02/IDEA_o1_v02.xml
  N_EVENTS 51000
  EXTRA_ARGS --gun.distribution uniform --gun.particle geantino --random.seed 1988301045
  PROPERTIES TIMEOUT 360)

#--------------------------------------------------
# test for IDEA o1 v03
if(DCH_INFO_H_EXIST)
  add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml
    EXTRA_ARGS --gun.distribution uniform --gun.particle geantino --random.seed 1988301045
      --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../example/SteeringFile_IDEA_o1_v03.py
    PROPERTIES TIMEOUT 900)

  # test for IDEA o1 v03, with DRC, with electrons
  add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/compact/IDEA_withDRC_o1_v03.xml
    EXTRA_ARGS --gun.distribution uniform --gun.particle e- --random.seed 1988301045
      --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../example/SteeringFile_IDEA_o1_v03.py
    PROPERTIES TIMEOUT 900)
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
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ALLEGRO/compact/ALLEGRO_o1_v02/ALLEGRO_o1_v02.xml
  EXTRA_ARGS --gun.distribution uniform --gun.particle geantino)

#--------------------------------------------------
# test for ALLEGRO o1 v03
if(DCH_INFO_H_EXIST)
  add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
    EXTRA_ARGS --gun.distribution uniform --gun.particle geantino
    PROPERTIES TIMEOUT 600)
endif()

#--------------------------------------------------
# test for ALLEGRO o2 v01
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ALLEGRO/compact/ALLEGRO_o2_v01/ALLEGRO_o2_v01.xml
  EXTRA_ARGS --gun.distribution uniform --gun.particle geantino)

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
# test for FCChh
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../FCChh/compact/FCChhBaseline/FCChh_DectMaster.xml
  EXTRA_ARGS --gun.distribution uniform --gun.particle geantino)

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
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../MuColl/MuColl/compact/MuColl_v1/MuColl_v1.xml)
add_ddsim_gun_test(${CMAKE_CURRENT_SOURCE_DIR}/../MuColl/MuSIC/compact/MuSIC_v2/MuSIC_v2.xml)

#--------------------------------------------------
# check if files named the same contain the same in FCCee
add_test(NAME t_test_files_versions_FCCee
    COMMAND "${PROJECT_SOURCE_DIR}/utils/compareIdenticalFiles.sh"
            "${PROJECT_SOURCE_DIR}/FCCee"
            "${PROJECT_SOURCE_DIR}/utils/IdenticalFiles_ToBeIgnored.txt"
)


#--------------------------------------------------
# detector name checks
function(add_detector_name_check compact_file)
  cmake_path(GET compact_file STEM det_name)
  add_test(NAME t_name_check_${det_name}
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../utils/detector_name_check.py
            -c ${compact_file} -n ${det_name})
endfunction()

# Add name checks for all files matching a given GLOB pattern.
function(add_detector_name_checks_for glob_pattern)
  file(GLOB compact_files LIST_DIRECTORIES false "${glob_pattern}")
  foreach(compact_file ${compact_files})
    add_detector_name_check(${compact_file})
  endforeach()
endfunction()

add_detector_name_checks_for(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ILD_FCCee/compact/*/ILD_FCCee_v0?.xml)
add_detector_name_checks_for(${CMAKE_CURRENT_SOURCE_DIR}/../MuColl/MAIA/compact/*/MAIA_v*.xml)
add_detector_name_checks_for(${CMAKE_CURRENT_SOURCE_DIR}/../MuColl/MuSIC/compact/*/MuSIC_v*.xml)
add_detector_name_checks_for(${CMAKE_CURRENT_SOURCE_DIR}/../ILD/compact/ILD_sl5_v02/ILD_l5_o?_v02.xml)
add_detector_name_checks_for(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/CLD/compact/*/CLD_o?_v0?.xml)
# add_detector_name_checks_for(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/ALLEGRO/compact/*/ALLEGRO_o?_v0?.xml)
# add_detector_name_checks_for(${CMAKE_CURRENT_SOURCE_DIR}/../FCCee/IDEA/compact/*/IDEA_o?_v0?.xml)


set(test_environment "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}")
get_property(test_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)
set_tests_properties(${test_names} PROPERTIES ENVIRONMENT "${test_environment}"
                                              FAIL_REGULAR_EXPRESSION " Exception; EXCEPTION;ERROR;Error")
