cmake_minimum_required(VERSION 3.5.1)
project( FiniteElementAnalysis )
find_package( Boost )
find_package( GSL )
find_package( GTest )
find_package( benchmark )

# since 3.1
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_BUILD_TYPE Debug)

# add profiler flags for performance analysis
set (CMAKE_CXX_FLAGS "-Ofast")
#set (CMAKE_CXX_FLAGS "-O3 -Ofast")

#LIST(APPEND CMAKE_MODULE_PATH "{CMAKE_CURRENT_SOURCE_DIR}/mathjoy")
#add_subdirectory( mathjoy )
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})
find_package( SFML 2 COMPONENTS network audio graphics window system )

include_directories( ${Boost_INCLUDE_DIRS})
include_directories( ${GSL_INCLUDE_DIRS} )
include_directories( ${SFML_INCLUDE_DIR} )
include_directories( ${CMAKE_SOURCE_DIR}/include )
#include_directories( ${PROJECT_SOURCE_DIR}/mathjoy/include/solver )

# get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
# foreach(dir ${dirs})
#   message(STATUS "dir='${dir}'")
# endforeach()
message(STATUS ${SFML_INCLUDE_DIR} ${SFML_LIBRARIES} ${SFML_DEPENDENCIES} ${Boost_LIBRARIES} )

LIST(APPEND SOURCES
  src/utils/utils.cpp
  src/utils/utils_function.cpp
  src/mesh_generator/mesh_generator.cpp

  src/matrix_io.cpp
  src/FEM_linear.cpp
  src/FEM_general.cpp
  src/FEM_2dim.cpp
  src/FEM_2dim_triangular.cpp
  src/quadrature.cpp

  src/FVM_2dim.cpp
  )

add_library( FEM_FVM STATIC ${SOURCES} )

add_executable(
  testing_general.out
  "src/test/main.cpp" )

add_executable(
  testing_linear.out
  "src/test/main_linear.cpp" )

add_executable(
  testing.out
  "src/test/test.cpp" )

# add_executable(
#   test_delaunary.out
#   src/test/test_delaunay.cpp )

add_executable(
  test_FVM.out
  src/test/test_FVM.cpp )

add_executable(
  test_running
  test/test_running.cpp )

add_executable(
  run_benchmark
  src/benchmark/miscellaneous)

target_link_libraries( testing_general.out FEM_FVM ${Boost_LIBRARIES} GSL::gsl GSL::gslcblas )
target_link_libraries( testing_linear.out FEM_FVM ${Boost_LIBRARIES} GSL::gsl GSL::gslcblas )
target_link_libraries( testing.out FEM_FVM ${Boost_LIBRARIES} GSL::gsl GSL::gslcblas )
#target_link_libraries( test_delaunary.out FEM_FVM ${SFML_LIBRARIES} ${SFML_DEPENDENCIES} )
target_link_libraries( test_FVM.out FEM_FVM ${SFML_LIBRARIES} ${SFML_DEPENDENCIES} GSL::gsl GSL::gslcblas )
target_link_libraries( test_running FEM_FVM ${GTEST_BOTH_LIBRARIES} )# ${GTEST_MAIN_LIBRARIES} ) #${GTEST_LIBRARIES} )
target_link_libraries( run_benchmark benchmark )
