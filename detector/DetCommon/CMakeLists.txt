################################################################################
# Package: DetCommon
################################################################################

file(GLOB sources src/*.cpp)

add_library(DetCommon SHARED ${sources})
dd4hep_generate_rootmap(DetCommon)
target_link_libraries(DetCommon DD4hep::DDCore DD4hep::DDG4 detectorSegmentations)
target_include_directories(DetCommon
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

file(GLOB headers include/DetCommon/*.h)
set_target_properties(DetCommon PROPERTIES PUBLIC_HEADER "${headers}")

# Destination directories are hardcoded because GNUdirectories are not included
install(TARGETS DetCommon
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_PREFIX}/include/DetCommon"
  COMPONENT dev)

# Destination directory is hardcoded because:
# 1. share must be defined because GNUdirectories are not included
# 2. FCCDetectors literal string replaced the labe of project name, kept as literal to keep consistency
install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION share/FCCDetectors/Detector/DetCommon)


#include(CTest)
#gaudi_add_test(DumpSimpleBox
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/dumpSimpleBox.py)
