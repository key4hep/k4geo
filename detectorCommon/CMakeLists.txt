################################################################################
# Package: detectorCommon
################################################################################

file(GLOB sources src/*.cpp)
add_library(detectorCommon SHARED ${sources})
target_link_libraries(detectorCommon DD4hep::DDCore DD4hep::DDG4 detectorSegmentations)
target_include_directories(detectorCommon
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

file(GLOB headers include/detectorCommon/*.h)
set_target_properties(detectorCommon PROPERTIES PUBLIC_HEADER "${headers}")

install(TARGETS detectorCommon
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_PREFIX}/include/detectorCommon"
  COMPONENT dev)


dd4hep_generate_rootmap(detectorCommon)

#include(CTest)
#gaudi_add_test(DumpSimpleBox
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/dumpSimpleBox.py)
