cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

if( POLICY CMP0167 )
  cmake_policy(SET CMP0167 NEW)
endif()

#---------------------------
set( PackageName k4geo )
#---------------------------

project(${PackageName} LANGUAGES CXX)

# project version
SET( ${PackageName}_VERSION_MAJOR 0 )
SET( ${PackageName}_VERSION_MINOR 26 )
SET( ${PackageName}_VERSION_PATCH 0 )

SET( ${PackageName}_VERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )
SET( ${PackageName}_SOVERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${PROJECT_SOURCE_DIR}/cmake )
set(LIBRARY_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

include(cmake/Key4hepConfig.cmake)

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
option(INSTALL_BEAMPIPE_STL_FILES "Download CAD files for building the detailed beampipe" OFF)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

find_package(DD4hep 1.31 REQUIRED COMPONENTS DDRec DDG4 DDParsers)

find_package ( ROOT REQUIRED COMPONENTS Geom GenVector)
message ( STATUS "ROOT_VERSION: ${ROOT_VERSION}" )

find_package( Geant4 REQUIRED )

# Deprecated: k4geo no longer uses LCIO. The option is kept temporarily so that
# existing configurations passing it do not fail outright.
OPTION(K4GEO_USE_LCIO "Deprecated and without effect: k4geo does not depend on LCIO anymore" OFF)
if(K4GEO_USE_LCIO)
  message(DEPRECATION "K4GEO_USE_LCIO is deprecated and has no effect, since k4geo does not depend on LCIO anymore. The option will be removed in a future release, please stop setting it.")
endif()

file(GLOB sources
  ./detector/tracker/*.cpp
  ./detector/calorimeter/*.cpp
  ./detector/calorimeter/dual-readout/src/*.cpp
  ./detector/calorimeter/dual-readout-tubes/src/*.cpp
  ./detector/fcal/*.cpp
  ./detector/muonSystem/*.cpp
  ./detector/other/*.cpp
  ./detector/CaloTB/*.cpp
  ./FCalTB/setup/*.cpp
  ./plugins/LinearSortingPolicy.cpp
  ./detector/PID/ARC_geo_o1_v01.cpp
  )

if(${DD4hep_VERSION} VERSION_LESS 1.29)
  set(FILES_DEPENDINGON_DCH_INFO_H "DriftChamber_o1_v02.cpp" )
  list(FILTER sources EXCLUDE REGEX "${FILES_DEPENDINGON_DCH_INFO_H}")
  message(WARNING "Subdetector ${FILES_DEPENDINGON_DCH_INFO_H} will not be built because the current version of DD4hep does not ship with the header file DDRec/DCH_info.h")
else()
   set(DCH_INFO_H_EXIST "TRUE" CACHE STRING "Indicates whether the drift chamber data extension header file was found (TRUE if found, non set otherwise).")

endif()

find_package(EDM4HEP)
file(GLOB G4sources
  ./plugins/TPCSDAction.cpp
  ./plugins/CaloPreShowerSDAction.cpp
  ./plugins/FiberDRCaloSDAction.h
  ./plugins/FiberDRCaloSDAction.cpp
  ./plugins/Geant4Output2EDM4hep_DRC.cpp
  ./plugins/DRCaloFastSimModel.cpp
  ./plugins/DRCaloFastSimModel.h
  ./plugins/DRTubesSDAction.hh
  ./plugins/DRTubesSDAction.cpp
  ./plugins/SCEPCal_MainSDAction.cpp
  ./plugins/SCEPCal_TimingSDAction.cpp
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

target_include_directories(${PackageName}   PRIVATE ${PROJECT_SOURCE_DIR}/detector/calorimeter/dual-readout/include )
target_include_directories(${PackageName}G4 PRIVATE ${PROJECT_SOURCE_DIR}/detector/calorimeter/dual-readout/include )

target_include_directories(${PackageName}   PRIVATE ${PROJECT_SOURCE_DIR}/detector/calorimeter/dual-readout-tubes/include )
target_include_directories(${PackageName}G4 PRIVATE ${PROJECT_SOURCE_DIR}/detector/calorimeter/dual-readout-tubes/include )

target_link_libraries(${PackageName}   DD4hep::DDCore DD4hep::DDRec DD4hep::DDParsers ROOT::Core detectorSegmentations detectorCommon)
target_link_libraries(${PackageName}G4 DD4hep::DDCore DD4hep::DDRec DD4hep::DDParsers DD4hep::DDG4 ROOT::Core detectorSegmentations detectorCommon podio::podioRootIO EDM4HEP::edm4hep ${Geant4_LIBRARIES})

#Create this_package.sh file, and install
dd4hep_instantiate_package(${PackageName})
if(EXISTS "${PROJECT_BINARY_DIR}/this${PackageName}.sh")
  message(STATUS "this${PackageName}.sh created successfully.")
  # Replace the K4GEO environment variable with the path to the source
  file(APPEND "${PROJECT_BINARY_DIR}/this${PackageName}.sh" "\nexport K4GEO=\"${CMAKE_CURRENT_SOURCE_DIR}\"")
else()
  message(WARNING "this${PackageName}.sh was not created.")
endif()

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
  add_subdirectory(test)
  set(BUILDNAME "${CMAKE_SYSTEM}-${CMAKE_CXX_COMPILER}-${CMAKE_BUILD_TYPE}" CACHE STRING "set build string for cdash")
endif(BUILD_TESTING)

#--- install remaining targets--------------------------
FILE(GLOB hfiles "ILD/include/*.h")
INSTALL(FILES ${hfiles}
  DESTINATION include/${PackageName}
  )

#--- install compact files------------------------------
if(INSTALL_COMPACT_FILES)
  INSTALL(DIRECTORY CaloTB CLIC FCalTB FCCee FCChh ILD fieldmaps MuColl SiD DESTINATION share/k4geo )
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

if(INSTALL_BEAMPIPE_STL_FILES)

    set(STL_FILES
          "AlBeMet162_30042024.stl"
          "Copper_pipe_28092023.stl"
          "Gold_19042024.stl"
          "Paraffine_19042024.stl"
          "Tungsten_mask_02102023.stl"
          "Water_30042024.stl"
        )
    # Set main FCC url
    set(FCC_URL "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco")
    set(STL_PATH "stl_files/Pipe_240430")

    # Set the output directory where the file will be placed
    set(OUTPUT_DIR "${PROJECT_SOURCE_DIR}/FCCee/MDI/compact/${STL_PATH}")
    file(MAKE_DIRECTORY ${OUTPUT_DIR})

    foreach(STL_FILE ${STL_FILES})
        set(FULL_URL "${FCC_URL}/MDI/${STL_PATH}/${STL_FILE}")
        message(DEBUG "Downloading file ${FULL_URL}")
        set(OUTPUT_FILE "${OUTPUT_DIR}/${STL_FILE}")

        if(EXISTS "${OUTPUT_FILE}")
          message(STATUS "File ${STL_FILE} already exists. Skipping download.")
        else()
          # Download the file
          file(DOWNLOAD ${FULL_URL} ${OUTPUT_FILE}
            SHOW_PROGRESS
            STATUS download_status)

          list(GET download_status 0 status_code)
          if(NOT status_code EQUAL 0)
              list(GET download_status 1 error_message)
              message(FATAL_ERROR "Failed to download file: ${error_message}")
          endif()
        endif()
    endforeach()

    file(MAKE_DIRECTORY share/k4geo/FCCee)
    INSTALL(DIRECTORY ${OUTPUT_DIR} DESTINATION share/k4geo/FCCee/MDI/compact/stl_files )

endif()
