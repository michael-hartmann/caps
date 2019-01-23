cmake_minimum_required (VERSION 3.0)
project (libcasimir)

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
include_directories("include/")
include_directories("cquadpack/include/")
include_directories("libhodlr/include/")

# optimization
set(OPT, "-O3")

# set CFLAGS
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 ${OPT} -g -Wall -Wextra -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wno-unused-parameter -fstrict-aliasing")

# set CXXFLAGS
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${OPT} -Wall -g")


# libhodlr
add_library(hodlr SHARED libhodlr/src/hodlr.cpp libhodlr/src/HODLR_Matrix.cpp libhodlr/src/HODLR_Node.cpp libhodlr/src/HODLR_Tree.cpp libhodlr/src/HODLR_Tree_NonSPD.cpp libhodlr/src/HODLR_Tree_SPD.cpp libhodlr/src/KDTree.cpp)
set_target_properties(hodlr PROPERTIES COMPILE_FLAGS "-DUSE_DOUBLE -Wno-unknown-pragmas")


# libcasimir
add_library(casimir SHARED bessel.c cache.c fcqs.c integration.c libcasimir.c logfac.c material.c matrix.c misc.c plm.c psd.c utils.c cquadpack/src/dqage.c cquadpack/src/dqagi.c cquadpack/src/dqags.c cquadpack/src/dqext.c cquadpack/src/dqk15.c cquadpack/src/dqk15i.c cquadpack/src/dqk21.c cquadpack/src/dqk31.c cquadpack/src/dqk41.c cquadpack/src/dqk51.c cquadpack/src/dqk61.c cquadpack/src/dqsort.c)

#execute_process(
#  COMMAND uname -smnr
#  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#  OUTPUT_VARIABLE MACHINE
#  OUTPUT_STRIP_TRAILING_WHITESPACE
#)
#add_definitions("-DMACHINE=${MACHINE}")

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

  message("GIT_HEAD=${GIT_COMMIT_HASH}")
  message("GIT_BRANCH=${GIT_BRANCH}")

  add_definitions(casimir -DGIT_HEAD="${GIT_COMMIT_HASH}" -DGIT_BRANCH="{GIT_BRANCH}")
endif(EXISTS "${CMAKE_SOURCE_DIR}/../.git")


# casimir frontend exectuable
add_executable(frontend casimir.c)
set_target_properties(frontend PROPERTIES OUTPUT_NAME "casimir")

target_link_libraries(frontend m)
target_link_libraries(frontend ${MPI_LIBRARIES})
target_link_libraries(frontend ${BLAS_LIBRARIES})
target_link_libraries(frontend ${LAPACK_LIBRARIES})
target_link_libraries(frontend ${CMAKE_SOURCE_DIR}/libcasimir.so)
target_link_libraries(frontend ${CMAKE_SOURCE_DIR}/libhodlr.so)


# casimir_logdetD
add_executable(casimir_logdetD casimir_logdetD.c)

target_link_libraries(casimir_logdetD m)
target_link_libraries(casimir_logdetD ${BLAS_LIBRARIES})
target_link_libraries(casimir_logdetD ${LAPACK_LIBRARIES})
target_link_libraries(casimir_logdetD ${CMAKE_SOURCE_DIR}/libcasimir.so)
target_link_libraries(casimir_logdetD ${CMAKE_SOURCE_DIR}/libhodlr.so)


# cylinder
#add_executable(cylinder cylinder.cpp argparse.c)

#target_link_libraries(cylinder ${BLAS_LIBRARIES})
#target_link_libraries(cylinder ${LAPACK_LIBRARIES})
#target_link_libraries(cylinder ${CMAKE_SOURCE_DIR}/libcasimir.so)
#target_link_libraries(cylinder ${CMAKE_SOURCE_DIR}/libhodlr.so)
