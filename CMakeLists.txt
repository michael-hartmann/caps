cmake_minimum_required (VERSION 3.0)
project (libcasimir)

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
include_directories("include/")
include_directories("cquadpack/include/")
include_directories("libhodlr/include/")

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
add_library(hodlr SHARED libhodlr/src/hodlr.cpp libhodlr/src/HODLR_Matrix.cpp libhodlr/src/HODLR_Node.cpp libhodlr/src/HODLR_Tree.cpp libhodlr/src/HODLR_Tree_NonSPD.cpp libhodlr/src/HODLR_Tree_SPD.cpp libhodlr/src/KDTree.cpp)
target_compile_definitions(hodlr PRIVATE -DUSE_DOUBLE)
set_target_properties(hodlr PROPERTIES COMPILE_FLAGS "-Wno-unknown-pragmas")

# libcasimir
add_library(casimir SHARED bessel.c cache.c fcqs.c integration.c libcasimir.c logfac.c material.c matrix.c misc.c plm.c psd.c utils.c cquadpack/src/dqage.c cquadpack/src/dqagi.c cquadpack/src/dqags.c cquadpack/src/dqext.c cquadpack/src/dqk15.c cquadpack/src/dqk15i.c cquadpack/src/dqk21.c cquadpack/src/dqk31.c cquadpack/src/dqk41.c cquadpack/src/dqk51.c cquadpack/src/dqk61.c cquadpack/src/dqsort.c)

target_compile_definitions(casimir PRIVATE LIBCASIMIR_VERSION="0.4.3-dev")

# get machine name
execute_process(
  COMMAND uname -smnr
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE MACHINE
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
target_compile_definitions(casimir PRIVATE -DMACHINE="${MACHINE}")

# get GIT head and branch
if(EXISTS "${CMAKE_SOURCE_DIR}/../.git")
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

  target_compile_definitions(casimir PRIVATE -DGIT_HEAD="${GIT_COMMIT_HASH}" -DGIT_BRANCH="${GIT_BRANCH}")
endif(EXISTS "${CMAKE_SOURCE_DIR}/../.git")


# casimir frontend exectuable
add_executable(frontend casimir.c)
set_target_properties(frontend PROPERTIES OUTPUT_NAME "casimir")

target_link_libraries(frontend m)
target_link_libraries(frontend ${MPI_LIBRARIES})
target_link_libraries(frontend ${BLAS_LIBRARIES})
target_link_libraries(frontend ${LAPACK_LIBRARIES})
target_link_libraries(frontend ${CMAKE_BINARY_DIR}/libcasimir.so)
target_link_libraries(frontend ${CMAKE_BINARY_DIR}/libhodlr.so)


# casimir_logdetD
add_executable(casimir_logdetD casimir_logdetD.c)

target_link_libraries(casimir_logdetD m)
target_link_libraries(casimir_logdetD ${BLAS_LIBRARIES})
target_link_libraries(casimir_logdetD ${LAPACK_LIBRARIES})
target_link_libraries(casimir_logdetD ${CMAKE_BINARY_DIR}/libcasimir.so)
target_link_libraries(casimir_logdetD ${CMAKE_BINARY_DIR}/libhodlr.so)


# cylinder
add_executable(cylinder cylinder.cpp argparse.c)
set_target_properties(cylinder PROPERTIES EXCLUDE_FROM_ALL 1)

target_link_libraries(cylinder ${BLAS_LIBRARIES})
target_link_libraries(cylinder ${LAPACK_LIBRARIES})
target_link_libraries(cylinder ${CMAKE_BINARY_DIR}/libcasimir.so)
target_link_libraries(cylinder ${CMAKE_BINARY_DIR}/libhodlr.so)

# tests
add_executable(tests tests/test_bessels.c tests/test_fresnel.c tests/test_Lambda.c tests/test_lfac.c tests/test_logdetD.c tests/test_logi.c tests/test_mie.c tests/test_mie_drude.c tests/test_Plm.c tests/tests.c tests/unittest.c)
set_target_properties(tests PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(tests PROPERTIES OUTPUT_NAME "casimir_tests")

target_link_libraries(tests m)
target_link_libraries(tests ${BLAS_LIBRARIES})
target_link_libraries(tests ${LAPACK_LIBRARIES})
target_link_libraries(tests ${CMAKE_BINARY_DIR}/libcasimir.so)
target_link_libraries(tests ${CMAKE_BINARY_DIR}/libhodlr.so)


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
