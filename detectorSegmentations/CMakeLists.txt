#################################################################################
##Package : detectorSegmentations
#################################################################################

file(GLOB sources src/*.cpp)

add_library(detectorSegmentations SHARED ${sources} )
add_library(k4geo::detectorSegmentations ALIAS detectorSegmentations)
dd4hep_generate_rootmap(detectorSegmentationsPlugin)
target_link_libraries(detectorSegmentations DD4hep::DDCore)
target_include_directories(detectorSegmentations
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

file(GLOB headers include/detectorSegmentations/*.h)
set_target_properties(detectorSegmentations PROPERTIES PUBLIC_HEADER "${headers}")

file(GLOB pluginsources src/plugins/*.cpp)
add_library(detectorSegmentationsPlugin SHARED ${pluginsources})
target_link_libraries(detectorSegmentationsPlugin detectorSegmentations DD4hep::DDCore)
target_include_directories(detectorSegmentationsPlugin
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)



# Destination directories are hardcoded because GNUdirectories are not included
install(TARGETS detectorSegmentations detectorSegmentationsPlugin
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIDIR} COMPONENT shlib
  PUBLIC_HEADER DESTINATION include/detectorSegmentations COMPONENT dev
)

#
#include(CTest)
#gaudi_add_test(TestSegmentationPhiEta
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/phiEtaSegmentation.py)
