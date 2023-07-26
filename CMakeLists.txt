cmake_minimum_required (VERSION 3.12)
project(pgenfitpythonized LANGUAGES CXX)

cmake_policy(SET CMP0086 NEW)
cmake_policy(SET CMP0078 NEW)
set (CMAKE_CXX_STANDARD 11)


### import modules
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} /usr/include/python3.6m)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Core RIO RooFit RooFitCore RooStats Hist Tree Matrix Physics MathCore Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})


#include source and include dir
include_directories(${PROJECT_SOURCE_DIR}/include)


#---Create a shared library with geneated dictionary
ROOT_GENERATE_DICTIONARY(G__fitF ${PROJECT_SOURCE_DIR}/include/fitF.hh LINKDEF ${PROJECT_SOURCE_DIR}/include/fitFLinkDef.hh)
add_library(fitF SHARED ${PROJECT_SOURCE_DIR}/src/fitF.cc ${PROJECT_SOURCE_DIR}/src/fitF_auxiliary.cc ${PROJECT_SOURCE_DIR}/src/fitF_cal.cc G__fitF.cxx)
target_link_libraries(fitF ${ROOT_LIBRARIES})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#---Create  a main program using the library
add_executable(main
    ${sources} ${headers}
    main.cc)
target_link_libraries(main ${ROOT_LIBRARIES})
target_link_libraries(main fitF)

add_executable(mainbinfit
    ${sources} ${headers}
    mainbinfit.cc)
target_link_libraries(mainbinfit ${ROOT_LIBRARIES})
target_link_libraries(mainbinfit fitF)

add_executable(makefunction
    ${sources} ${headers}
    makefunction.cc)
target_link_libraries(makefunction ${ROOT_LIBRARIES})
target_link_libraries(makefunction fitF)

add_executable(simulationC
    ${sources} ${headers}
    mainsimulation.cc)
target_link_libraries(simulationC ${ROOT_LIBRARIES})
target_link_libraries(simulationC fitF)

ROOT_GENERATE_DICTIONARY(G__decaypath  ${PROJECT_SOURCE_DIR}/include/decaypath.hh  ${PROJECT_SOURCE_DIR}/include/common.hh LINKDEF  ${PROJECT_SOURCE_DIR}/include/decaypathLinkDef.hh)
add_library(decaypath SHARED  ${PROJECT_SOURCE_DIR}/include/decaypath.hh  ${PROJECT_SOURCE_DIR}/src/decaypath.cc G__decaypath.cxx)
target_link_libraries(decaypath ${ROOT_LIBRARIES})


set_property(SOURCE ${PROJECT_SOURCE_DIR}/include/simulation.i PROPERTY CPLUSPLUS ON)
swig_add_library(simulation
  TYPE SHARED
  LANGUAGE python
  SOURCES ${PROJECT_SOURCE_DIR}/include/simulation.i  ${PROJECT_SOURCE_DIR}/src/simulation.cc
  )
swig_link_libraries(simulation decaypath)
target_include_directories(simulation PUBLIC "")

#----------------------------------------------------------------------------
# Copy all resources including scripts, examples, text files, etc... to the build directory, i.e. the directory in which we
# build project. This is so that we can us the resources because it
# relies on these scripts being in the current working directory.
#
set(PROJ_DATA
  simparmsex.txt
  parmsex.txt
  effparmsex.txt
  testdata.root
  tdiffdataex.root
  FRDM-QRPA12-halflife.txt
  stable.csv
  anasim.C
  genFunctionCopy.sh
  runsinglefit.sh
  runbenchmark.sh
  rununbinfit_multithreads.sh
  )


foreach(_script ${PROJ_DATA})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
