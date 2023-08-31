################################################################################
# Package: DetSensitive
################################################################################

file(GLOB sources src/*.cpp)

add_library (DetSensitive SHARED ${sources})
dd4hep_generate_rootmap(DetSensitive)
target_link_libraries(DetSensitive DD4hep::DDCore DD4hep::DDG4 DetCommon)
target_include_directories(DetSensitive
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

file(GLOB headers include/DetSensitive/*.h)
set_target_properties(DetSensitive PROPERTIES PUBLIC_HEADER "${headers}")

# Destination directories are hardcoded because GNUdirectories are not included
install(TARGETS DetSensitive
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_PREFIX}/include/DetSensitive"
  )

#include(CTest)
#gaudi_add_test(ExternalTrackerSD
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Detector/DetSensitive/tests/
#               FRAMEWORK tests/options/testSimpleTrackerSD.py)
#gaudi_add_test(ExternalCalorimeterSD
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Detector/DetSensitive/tests/
#               FRAMEWORK tests/options/testSimpleCaloSD.py)
#gaudi_add_test(FullParticleAbsorptionSD
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Detector/DetSensitive/tests/
#               FRAMEWORK tests/options/testDd4hepFullParticleAbsorptionSD.py)
