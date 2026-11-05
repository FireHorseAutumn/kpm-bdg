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

add_executable(test_r2 tests/test_r2.cpp)
target_link_libraries(test_r2 PRIVATE bdg)
target_compile_options(test_r2 PRIVATE -Wall -Wextra)
add_test(NAME r2 COMMAND test_r2)

add_executable(test_r3 tests/test_r3.cpp)
target_link_libraries(test_r3 PRIVATE bdg)
target_compile_options(test_r3 PRIVATE -Wall -Wextra)
# Same read-only cross-gate as test_r1: the delta->0 L=14 gate reads
# kpm-tightbinding's committed exact_dos_L14.csv (DESIGN.md section 3 R3).
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../kpm-tightbinding/data)
  target_compile_definitions(test_r3 PRIVATE
    KPM_TB_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../kpm-tightbinding/data")
endif()
add_test(NAME r3 COMMAND test_r3)

add_executable(test_r4 tests/test_r4.cpp)
target_link_libraries(test_r4 PRIVATE bdg)
target_compile_options(test_r4 PRIVATE -Wall -Wextra)
# Same read-only cross-gates: the R1 L=14 spectrum gate and the R3
# delta->0 L=14 DOS gate re-run at the d-wave delta (DESIGN.md §3 R4).
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../kpm-tightbinding/data)
  target_compile_definitions(test_r4 PRIVATE
    KPM_TB_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../kpm-tightbinding/data")
endif()
add_test(NAME r4 COMMAND test_r4)

# The R3 figure generator (DESIGN.md section 3: the figure is the
# deliverable). Not a test; writes data/dos_swave_L128.csv.
add_executable(dos_figure src/dos_figure.cpp)
target_link_libraries(dos_figure PRIVATE bdg)
target_compile_options(dos_figure PRIVATE -Wall -Wextra)
target_compile_definitions(dos_figure PRIVATE
  BDG_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
