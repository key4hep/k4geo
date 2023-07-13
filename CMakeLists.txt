cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

#---------------------------
set( PackageName lcgeo )
#---------------------------

project(${PackageName})

# project version
SET( ${PackageName}_VERSION_MAJOR 0 )
SET( ${PackageName}_VERSION_MINOR 18 )
SET( ${PackageName}_VERSION_PATCH 1 )

SET( ${PackageName}_VERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )
SET( ${PackageName}_SOVERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${PROJECT_SOURCE_DIR}/cmake )
set(LIBRARY_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

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
find_package( LCIO REQUIRED)


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

target_include_directories(${PackageName}   PRIVATE ${PROJECT_SOURCE_DIR}/detector/include )
target_include_directories(${PackageName}G4 PRIVATE ${PROJECT_SOURCE_DIR}/detector/include )

target_include_directories(${PackageName}   SYSTEM PUBLIC ${LCIO_INCLUDE_DIRS})
target_include_directories(${PackageName}G4 SYSTEM PUBLIC ${LCIO_INCLUDE_DIRS})

target_link_libraries(${PackageName}   DD4hep::DDCore DD4hep::DDRec DD4hep::DDParsers ${ROOT_LIBRARIES} ${LCIO_LIBRARIES})
target_link_libraries(${PackageName}G4 DD4hep::DDCore DD4hep::DDRec DD4hep::DDParsers DD4hep::DDG4 ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${LCIO_LIBRARIES})

#Create this_package.sh file, and install
dd4hep_instantiate_package(${PackageName})

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
  INSTALL(DIRECTORY CaloTB CLIC FCalTB FCCee ILD fieldmaps SiD DESTINATION share/lcgeo/compact/ )
endif()
