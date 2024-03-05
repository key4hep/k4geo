cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

#---------------------------
set( PackageName k4geo )
#---------------------------

project(${PackageName})

# project version
SET( ${PackageName}_VERSION_MAJOR 0 )
SET( ${PackageName}_VERSION_MINOR 20 )
SET( ${PackageName}_VERSION_PATCH 0 )

SET( ${PackageName}_VERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )
SET( ${PackageName}_SOVERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${PROJECT_SOURCE_DIR}/cmake )
set(LIBRARY_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR lib)
set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

#------------- set the default installation directory to be the source directory

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET( CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR} CACHE PATH
    "install prefix path  - overwrite with -D CMAKE_INSTALL_PREFIX = ..." 
    FORCE )
  MESSAGE(STATUS "CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} - overwrite with -D CMAKE_INSTALL_PREFIX" )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)



#---Options-------------------------------------------------------------------------

option(BUILD_TESTING "Enable and build tests" ON)
option(INSTALL_COMPACT_FILES "Copy compact files to install area" OFF)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


find_package(DD4hep REQUIRED COMPONENTS DDRec DDG4 DDParsers)
# dd4hep_set_compiler_flags() expects DD4hep_SET_RPATH to be set to ON
# otherwise it will not set the rpath when installing
set(DD4HEP_SET_RPATH ON)
dd4hep_set_compiler_flags()

find_package ( ROOT REQUIRED COMPONENTS Geom GenVector)
message ( STATUS "ROOT_VERSION: ${ROOT_VERSION}" )

find_package( Geant4 REQUIRED ) 
OPTION(K4GEO_USE_LCIO "Enable or disable the use of LCIO, which is needed for some detector constructors and plugins" ON)
if(K4GEO_USE_LCIO)
  find_package(LCIO REQUIRED)
  # Shim for older LCIO versions
  if(NOT TARGET LCIO::lcio)
    add_library(LCIO::lcio INTERFACE IMPORTED GLOBAL)
    set_target_properties(LCIO::lcio
      PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${LCIO_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${LCIO_LIBRARIES}"
    )
  endif()
endif()

file(GLOB sources 
  ./detector/tracker/*.cpp
  ./detector/calorimeter/*.cpp
  ./detector/fcal/*.cpp
  ./detector/other/*.cpp
  ./detector/CaloTB/*.cpp
  ./FCalTB/setup/*.cpp
  ./plugins/LinearSortingPolicy.cpp
  ./detector/PID/ARC_geo_o1_v01.cpp
  )

if(NOT K4GEO_USE_LCIO)
  set(lcio_sources # in ./detector/tracker
    TrackerEndcap_o2_v05_geo.cpp
    SiTrackerEndcap_o2_v02ext_geo.cpp
    TrackerBarrel_o1_v03_geo.cpp
    TrackerBarrel_o1_v04_geo.cpp
    TrackerBarrel_o1_v05_geo.cpp
    TrackerEndcap_o1_v05_geo.cpp
    TrackerEndcap_o2_v06_geo.cpp
    VertexBarrel_detailed_o1_v01_geo.cpp
    VertexEndcap_o1_v05_geo.cpp
    ZPlanarTracker_geo.cpp
    )
  foreach(lcio_source ${lcio_sources})
    list(FILTER sources EXCLUDE REGEX "${lcio_source}")
  endforeach()
  message(STATUS "Use of LCIO is DISABLED, some detectors that depend on LCIO will not be built: ${lcio_sources}")
endif()

file(GLOB G4sources
  ./plugins/TPCSDAction.cpp
  ./plugins/CaloPreShowerSDAction.cpp
)

if(DD4HEP_USE_PYROOT)
  ROOT_GENERATE_DICTIONARY(G__ILDEx ${headers} LINKDEF include/ROOT/LinkDef.h)
  list(APPEND sources G__ILDEx.cxx)
endif()

add_dd4hep_plugin(${PackageName} SHARED ${sources})
add_dd4hep_plugin(${PackageName}G4 SHARED ${G4sources})

add_library(lcgeo ALIAS k4geo)

target_include_directories(${PackageName}   PRIVATE ${PROJECT_SOURCE_DIR}/detector/include )
target_include_directories(${PackageName}G4 PRIVATE ${PROJECT_SOURCE_DIR}/detector/include )

target_link_libraries(${PackageName}   DD4hep::DDCore DD4hep::DDRec DD4hep::DDParsers ROOT::Core detectorSegmentations)
target_link_libraries(${PackageName}G4 DD4hep::DDCore DD4hep::DDRec DD4hep::DDParsers DD4hep::DDG4 ROOT::Core ${Geant4_LIBRARIES})

if(K4GEO_USE_LCIO)
  target_link_libraries(${PackageName}   LCIO::lcio)
  target_link_libraries(${PackageName}G4 LCIO::lcio)
endif()

#Create this_package.sh file, and install
dd4hep_instantiate_package(${PackageName})

add_subdirectory(detectorSegmentations)
add_subdirectory(detectorCommon)

# Destination directories are hardcoded because GNUdirectories are not included
install(TARGETS ${PackageName}  ${PackageName}G4
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT shlib
)

#---Testing-------------------------------------------------------------------------
if(BUILD_TESTING)
  include(CTest)
  add_subdirectory(lcgeoTests)
  set(BUILDNAME "${CMAKE_SYSTEM}-${CMAKE_CXX_COMPILER}-${CMAKE_BUILD_TYPE}" CACHE STRING "set build string for cdash")
endif(BUILD_TESTING)

#--- install remaining targets--------------------------
FILE(GLOB hfiles "ILD/include/*.h")
INSTALL(FILES ${hfiles} 
  DESTINATION include/${PackageName} 
  )

#--- install compact files------------------------------
if(INSTALL_COMPACT_FILES)
  INSTALL(DIRECTORY CaloTB CLIC FCalTB FCCee ILD fieldmaps SiD DESTINATION share/k4geo )
endif()

# create k4geoConfig and friends

INSTALL(EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  FILE k4geoConfig-targets.cmake)

include(CMakePackageConfigHelpers)

configure_package_config_file(cmake/k4geoConfig.cmake.in ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_CMAKEDIR}/k4geoConfig.cmake
  PATH_VARS CMAKE_INSTALL_CMAKEDIR
  INSTALL_DESTINATION ${CMAKE_INSTALL_CMAKEDIR} )

write_basic_package_version_file(
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_CMAKEDIR}/k4geoConfigVersion.cmake
  VERSION ${${PackageName}_VERSION}
  COMPATIBILITY SameMajorVersion )

install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_CMAKEDIR}/k4geoConfig.cmake
              ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_CMAKEDIR}/k4geoConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_CMAKEDIR} )
