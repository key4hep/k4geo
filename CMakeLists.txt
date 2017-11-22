cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

#---------------------------
set( PackageName lcgeo )
#---------------------------

project(${PackageName})

# project version
SET( ${PackageName}_VERSION_MAJOR 0 )
SET( ${PackageName}_VERSION_MINOR 15 )
SET( ${PackageName}_VERSION_PATCH 3 )

SET( ${PackageName}_VERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )
SET( ${PackageName}_SOVERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${CMAKE_SOURCE_DIR}/cmake ) 
set(LIBRARY_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

#------------- set the default installation directory to be the source directory

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET( CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH  
    "install prefix path  - overwrite with -D CMAKE_INSTALL_PREFIX = ..." 
    FORCE )
  MESSAGE(STATUS "CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} - overwrite with -D CMAKE_INSTALL_PREFIX" )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)



#---Options-------------------------------------------------------------------------

option(BUILD_TESTING "Enable and build tests" ON)
option(CMAKE_MACOSX_RPATH "Build with rpath on macos" ON)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


find_package(DD4hep REQUIRED COMPONENTS DDSegmentation DDRec DDG4 DDParsers)
dd4hep_set_compiler_flags()

SET( DD4hep_REDUCED_COMPONENT ${DD4hep_COMPONENT_LIBRARIES} )
## Drop DDG4 from components for normal geometry plugins
LIST( REMOVE_AT DD4hep_REDUCED_COMPONENT -2 )

find_package ( ROOT REQUIRED COMPONENTS Geom GenVector)
include(${ROOT_USE_FILE})
message ( STATUS "ROOT_VERSION: ${ROOT_VERSION}" )

find_package( Boost REQUIRED ) 
include_directories( SYSTEM ${Boost_INCLUDE_DIRS} )

find_package( Geant4 REQUIRED ) 
include(${Geant4_USE_FILE})

find_package( LCIO REQUIRED)
#find_package( GEAR REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/detector/include )
include_directories( SYSTEM
  ${DD4hep_INCLUDE_DIRS}
  ${DD4hep_COMPONENTS_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${LCIO_INCLUDE_DIRS}
 )
#  ${GEAR_INCLUDE_DIRS}

file(GLOB sources 
  ./detector/tracker/*.cpp 
  ./detector/calorimeter/*.cpp 
  ./detector/fcal/*.cpp 
  ./detector/other/*.cpp 
  ./detector/CaloTB/*.cpp 
  ./FCalTB/setup/*.cpp
  )

file(GLOB G4sources
  ./plugins/TPCSDAction.cpp
  ./plugins/CaloPreShowerSDAction.cpp
)


#file(GLOB headers ./ILD/include/*.h)

if(DD4HEP_USE_PYROOT)
  ROOT_GENERATE_DICTIONARY(G__ILDEx ${headers} LINKDEF include/ROOT/LinkDef.h)
  list(APPEND sources G__ILDEx.cxx)
endif()

#Create target, library, rootmap, install

add_dd4hep_plugin(${PackageName} SHARED ${sources})

add_dd4hep_plugin(${PackageName}G4 SHARED ${G4sources})

target_link_libraries(${PackageName} ${DD4hep_LIBRARIES} ${DD4hep_REDUCED_COMPONENT} ${ROOT_LIBRARIES} ${LCIO_LIBRARIES})

target_link_libraries(${PackageName}G4 ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES} ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${LCIO_LIBRARIES})


#Create this_package.sh file, and install
dd4hep_instantiate_package(${PackageName})

#---Testing-------------------------------------------------------------------------
if(BUILD_TESTING)

  include(CTest)
  enable_testing()
  add_subdirectory(lcgeoTests)
  set(BUILDNAME "${CMAKE_SYSTEM}-${CMAKE_CXX_COMPILER}-${CMAKE_BUILD_TYPE}" CACHE STRING "set build string for cdash")

endif(BUILD_TESTING)

#--- install remaining targets--------------------------
FILE(GLOB hfiles "ILD/include/*.h")
INSTALL(FILES ${hfiles} 
  DESTINATION include/${PackageName} 
  )

INSTALL(DIRECTORY DDSim/DDSim
  DESTINATION lib/python
  FILES_MATCHING PATTERN *.py
  )

INSTALL(PROGRAMS DDSim/bin/ddsim
  DESTINATION bin
  )
