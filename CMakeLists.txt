cmake_minimum_required(VERSION 3.20)
project(kpm_bdg LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

enable_testing()

# Eigen backs the R1 exact-diagonalization oracle only (DESIGN.md §2). A
# CUDA target exists only if the R5 GPU decision is made (DESIGN.md §3).
find_package(Eigen3 REQUIRED NO_MODULE)

add_library(bdg INTERFACE)
target_include_directories(bdg INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(bdg INTERFACE Eigen3::Eigen)

add_executable(test_stub tests/test_stub.cpp)
target_link_libraries(test_stub PRIVATE bdg)
target_compile_options(test_stub PRIVATE -Wall -Wextra)
add_test(NAME stub COMMAND test_stub)

add_executable(test_r1 tests/test_r1.cpp)
target_link_libraries(test_r1 PRIVATE bdg)
target_compile_options(test_r1 PRIVATE -Wall -Wextra)
# Read-only cross-gate against kpm-tightbinding's committed L=14 spectrum
# when a sibling checkout is present (DESIGN.md §6 decision record).
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../kpm-tightbinding/data)
  target_compile_definitions(test_r1 PRIVATE
    KPM_TB_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../kpm-tightbinding/data")
endif()
add_test(NAME r1 COMMAND test_r1)
