#################################################################################
##Package : DetSegmentation
#################################################################################

file(GLOB sources src/*.cpp)

add_library(DetSegmentation SHARED ${sources})
dd4hep_generate_rootmap(DetSegmentationPlugin)
target_link_libraries(DetSegmentation DD4hep::DDCore)
target_include_directories(DetSegmentation
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

file(GLOB headers include/DetSegmentation/*.h)
set_target_properties(DetSegmentation PROPERTIES PUBLIC_HEADER "${headers}")

file(GLOB pluginsources src/plugins/*.cpp)
add_library(DetSegmentationPlugin SHARED ${pluginsources})
target_link_libraries(DetSegmentationPlugin DetSegmentation DD4hep::DDCore)
target_include_directories(DetSegmentationPlugin
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)



# Destination directories are hardcoded because GNUdirectories are not included
install(TARGETS DetSegmentation DetSegmentationPlugin
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_PREFIX}/include/DetSegmentation"
  COMPONENT dev)

#
#include(CTest)
#gaudi_add_test(TestSegmentationPhiEta
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/phiEtaSegmentation.py)
