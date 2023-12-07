project(ddDRcalo)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

dd4hep_add_plugin(ddDRcalo SOURCES src/*.cpp USES
  DD4hep::DDCore
  DD4hep::DDCond
  ROOT::Core
  ROOT::Geom
  ROOT::GenVector
  ROOT::MathCore
  DRsegmentation
)

target_include_directories(ddDRcalo PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(ddDRcalo PROPERTIES PUBLIC_HEADER "${headers}")

install(TARGETS ddDRcalo
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)

# install(DIRECTORY compact DESTINATION ${CMAKE_INSTALL_DATADIR})

dd4hep_configure_scripts( ddDRcalo DEFAULT_SETUP WITH_TESTS )
