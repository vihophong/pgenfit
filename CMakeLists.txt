# From https://root.cern.ch/faq/can-i-integrate-root-my-cmake-build

# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(main)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Core RIO RooFit RooFitCore RooStats Hist Tree Matrix Physics MathCore Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

#---Create a shared library with geneated dictionary
ROOT_GENERATE_DICTIONARY(G__fitF ${PROJECT_SOURCE_DIR}/include/fitF.hh LINKDEF ${PROJECT_SOURCE_DIR}/include/fitFLinkDef.hh)
add_library(fitF SHARED ${PROJECT_SOURCE_DIR}/src/fitF.cc ${PROJECT_SOURCE_DIR}/src/fitF_auxiliary.cc ${PROJECT_SOURCE_DIR}/src/fitF_cal.cc G__fitF.cxx)
target_link_libraries(fitF ${ROOT_LIBRARIES})


#include source and include dir
include_directories(${PROJECT_SOURCE_DIR}/include)
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

add_executable(simulation
    ${sources} ${headers}
    mainsimulation.cc)
target_link_libraries(simulation ${ROOT_LIBRARIES})
target_link_libraries(simulation fitF)

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
