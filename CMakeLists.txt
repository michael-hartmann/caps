cmake_minimum_required (VERSION 3.0)
project (CaPS)

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
include_directories("src/include/")
include_directories("src/cquadpack/include/")
include_directories("src/libhodlr/include/")

# set standard optimization -O3
if(NOT DEFINED OPT)
    set(OPT "-O3")
endif()

# set CFLAGS
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 ${OPT} -Wall -Wextra -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wno-unused-parameter -fstrict-aliasing")

# set CXXFLAGS
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${OPT} -Wall")

# By default, icc violates strict IEEE floating point behaviour - similar to
# GCC's --fast-math option. The code, however, relies on strict IEEE floating
# point behaviour. The option "-fp-model precise" sets the correct floating
# point behaviour.
# Also, when using the math kernel library (MKL), use the non-threaded library.
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fp-model precise -mkl=sequential")
endif()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model precise -mkl=sequential")
endif()

# libhodlr
add_library(hodlr SHARED src/libhodlr/src/hodlr.cpp src/libhodlr/src/HODLR_Matrix.cpp src/libhodlr/src/HODLR_Node.cpp src/libhodlr/src/HODLR_Tree.cpp src/libhodlr/src/HODLR_Tree_NonSPD.cpp src/libhodlr/src/HODLR_Tree_SPD.cpp src/libhodlr/src/KDTree.cpp)
target_compile_definitions(hodlr PRIVATE -DUSE_DOUBLE)
set_target_properties(hodlr PROPERTIES COMPILE_FLAGS "-Wno-unknown-pragmas")

# libcaps
add_library(caps SHARED src/bessel.c src/cache.c src/fcqs.c src/integration.c src/libcaps.c src/logfac.c src/material.c src/matrix.c src/misc.c src/plm.c src/psd.c src/utils.c src/cquadpack/src/dqage.c src/cquadpack/src/dqagi.c src/cquadpack/src/dqags.c src/cquadpack/src/dqext.c src/cquadpack/src/dqk15.c src/cquadpack/src/dqk15i.c src/cquadpack/src/dqk21.c src/cquadpack/src/dqk31.c src/cquadpack/src/dqk41.c src/cquadpack/src/dqk51.c src/cquadpack/src/dqk61.c src/cquadpack/src/dqsort.c)

target_compile_definitions(caps PRIVATE LIBCAPS_VERSION="0.4.3-dev")

# get machine name
execute_process(
  COMMAND uname -smnr
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE MACHINE
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
target_compile_definitions(caps PRIVATE -DMACHINE="${MACHINE}")

# get GIT head and branch
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  target_compile_definitions(caps PRIVATE -DGIT_HEAD="${GIT_COMMIT_HASH}" -DGIT_BRANCH="${GIT_BRANCH}")
endif(EXISTS "${CMAKE_SOURCE_DIR}/.git")


# caps frontend exectuable
add_executable(caps_frontend src/caps.c)
set_target_properties(caps_frontend PROPERTIES OUTPUT_NAME "caps")

target_link_libraries(caps_frontend m)
target_link_libraries(caps_frontend ${MPI_LIBRARIES})
target_link_libraries(caps_frontend ${BLAS_LIBRARIES})
target_link_libraries(caps_frontend ${LAPACK_LIBRARIES})
target_link_libraries(caps_frontend ${CMAKE_BINARY_DIR}/libcaps.so)
target_link_libraries(caps_frontend ${CMAKE_BINARY_DIR}/libhodlr.so)


# caps_logdetD
add_executable(caps_logdetD src/caps_logdetD.c)

target_link_libraries(caps_logdetD m)
target_link_libraries(caps_logdetD ${BLAS_LIBRARIES})
target_link_libraries(caps_logdetD ${LAPACK_LIBRARIES})
target_link_libraries(caps_logdetD ${CMAKE_BINARY_DIR}/libcaps.so)
target_link_libraries(caps_logdetD ${CMAKE_BINARY_DIR}/libhodlr.so)


# capc
add_executable(capc src/capc.cpp src/argparse.c)
set_target_properties(capc PROPERTIES EXCLUDE_FROM_ALL 1)

target_link_libraries(capc ${BLAS_LIBRARIES})
target_link_libraries(capc ${LAPACK_LIBRARIES})
target_link_libraries(capc ${CMAKE_BINARY_DIR}/libcaps.so)
target_link_libraries(capc ${CMAKE_BINARY_DIR}/libhodlr.so)


# cass
add_executable(cass src/cass.c src/argparse.c)
set_target_properties(cass PROPERTIES EXCLUDE_FROM_ALL 1)
target_link_libraries(cass m)

target_link_libraries(cass ${BLAS_LIBRARIES})
target_link_libraries(cass ${LAPACK_LIBRARIES})
target_link_libraries(cass ${CMAKE_BINARY_DIR}/libcaps.so)
target_link_libraries(cass ${CMAKE_BINARY_DIR}/libhodlr.so)


# tests
add_executable(tests src/tests/test_bessels.c src/tests/test_fresnel.c src/tests/test_Lambda.c src/tests/test_lfac.c src/tests/test_logdetD.c src/tests/test_logi.c src/tests/test_mie.c src/tests/test_mie_drude.c src/tests/test_Plm.c src/tests/tests.c src/tests/unittest.c)
set_target_properties(tests PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(tests PROPERTIES OUTPUT_NAME "caps_tests")

target_link_libraries(tests m)
target_link_libraries(tests ${BLAS_LIBRARIES})
target_link_libraries(tests ${LAPACK_LIBRARIES})
target_link_libraries(tests ${CMAKE_BINARY_DIR}/libcaps.so)
target_link_libraries(tests ${CMAKE_BINARY_DIR}/libhodlr.so)


# examples
# AuAl
add_executable(AuAl examples/AuAl.c)
set_target_properties(AuAl PROPERTIES EXCLUDE_FROM_ALL 1)

target_link_libraries(AuAl ${BLAS_LIBRARIES})
target_link_libraries(AuAl ${LAPACK_LIBRARIES})
target_link_libraries(AuAl ${CMAKE_BINARY_DIR}/libcaps.so)
target_link_libraries(AuAl ${CMAKE_BINARY_DIR}/libhodlr.so)

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("optimization:      " ${OPT})
message("C compiler:        " ${CMAKE_C_COMPILER})
message("C flags:           " ${CMAKE_C_FLAGS})
message("C++ compiler:      " ${CMAKE_CXX_COMPILER})
message("C++ flags:         " ${CMAKE_CXX_FLAGS})
message("blas libraries:    " ${BLAS_LIBRARIES})
message("lapack libraries:  " ${LAPACK_LIBRARIES})
message("latest git commit: " ${GIT_COMMIT_HASH})
message("git branch:        " ${GIT_BRANCH})
message("machine:           " ${MACHINE})
message("build directory:   " ${CMAKE_BINARY_DIR})
message("")
