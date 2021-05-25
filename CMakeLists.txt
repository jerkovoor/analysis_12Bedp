cmake_minimum_required(VERSION 2.8)
project(12BeAnalysis)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MACOSX_RPATH 1)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O0 -g")

#Find ROOT
find_package(ROOT REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

set(SCRIPTS
    EnergyLoss.dat
    Be_in_Ag.dat
    Be_in_Al.dat
    Be_in_B.dat
    Be_in_D2_density0_201.dat
    Be_in_P.dat
    Be_in_Si.dat
    Be_in_SiO2.dat
    Proton_Aluminum.dat
    Proton_Boron.dat
    Proton_DeuteriumTarget.dat
    geometry_s1506.txt
    Sd1rAlphaCalibration.txt
    Sd2rInBeamCalibration.txt
    )

foreach(_script ${SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
        COPYONLY
    )
endforeach()

ROOT_GENERATE_DICTIONARY(G__Analysis ${CMAKE_CURRENT_SOURCE_DIR}/Analysis.h LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/LinkDefAnalysis.h)
add_library(AnalysisLibrary SHARED Analysis.cpp G__Analysis.cxx)
target_link_libraries(AnalysisLibrary ${ROOT_LIBRARIES} Spectrum)
add_executable(analysis AnalysisRun.cpp)
target_link_libraries(analysis AnalysisLibrary)