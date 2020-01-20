##Test
MESSAGE("Hello, world")
IF(DEFINED APPLE)
    MESSAGE("APPLE")
ENDIF(DEFINED APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")

set(CMAKE_OSX_DEPLOYMENT_TARGET "")
INCLUDE(SetCompileFlag.cmake)

cmake_minimum_required(VERSION 2.6)
project(directPsi CXX Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

## include CERN ROOT
include(FindROOT.cmake)
MESSAGE(${ROOT_LIBRARIES})
include(FindEvtGen.cmake)
include(FindHepMC.cmake)


## set fortran flags
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-fno-automatic")


set(RRF_DIR /Users/luchinsky/Work/rrF)
include_directories(~/local/include)
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${EVTGEN_INCLUDE_DIR}) 
include_directories(${RRF_DIR}/src) 

include_directories(./)
include_directories(src/)




#adding sources for segR_omega.exe
add_executable (segR_omega.exe "src/segR_omega.cc" "src/EvtWHad.cpp" "src/EvtXiccXicNpi.cpp" "src/a.for")
target_link_libraries(segR_omega.exe ${ROOT_LIBRARIES} ${EVTGEN_LIBRARY_DIR}/archive/libEvtGen.a ${HEPMC_LIBRARY_DIR}/libHepMC.a ${RRF_DIR}/build/libsegR.a)
