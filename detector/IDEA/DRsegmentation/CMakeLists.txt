project(DRsegmentation)

file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cpp
)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

add_library(DRsegmentation SHARED ${sources} ${headers})

target_include_directories(DRsegmentation PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(DRsegmentation PROPERTIES PUBLIC_HEADER "${headers}")

target_link_libraries(
  DRsegmentation
  DD4hep::DDCore
)

dd4hep_generate_rootmap(DRsegmentation)

install(TARGETS DRsegmentation EXPORT DetectorTargets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)
