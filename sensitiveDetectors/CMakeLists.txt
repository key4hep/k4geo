################################################################################
# Package: sensitiveDetectors
################################################################################

file(GLOB sources src/*.cpp)

add_library (sensitiveDetectors SHARED ${sources})
dd4hep_generate_rootmap(sensitiveDetectors)
target_link_libraries(sensitiveDetectors DD4hep::DDCore DD4hep::DDG4 DetCommon)
target_include_directories(sensitiveDetectors
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

file(GLOB headers include/sensitiveDetectors/*.h)
set_target_properties(sensitiveDetectors PROPERTIES PUBLIC_HEADER "${headers}")

# Destination directories are hardcoded because GNUdirectories are not included
install(TARGETS sensitiveDetectors
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sensitiveDetectors"
  )

#include(CTest)
#gaudi_add_test(ExternalTrackerSD
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Detector/sensitiveDetectors/tests/
#               FRAMEWORK tests/options/testSimpleTrackerSD.py)
#gaudi_add_test(ExternalCalorimeterSD
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Detector/sensitiveDetectors/tests/
#               FRAMEWORK tests/options/testSimpleCaloSD.py)
#gaudi_add_test(FullParticleAbsorptionSD
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Detector/sensitiveDetectors/tests/
#               FRAMEWORK tests/options/testDd4hepFullParticleAbsorptionSD.py)
