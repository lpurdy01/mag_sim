cmake_minimum_required(VERSION 3.16)
project(motor_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_compile_options(-Wall -Wextra -Wpedantic)

find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Project include path
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source list (expand later)
file(GLOB SRC_FILES CONFIGURE_DEPENDS
     "${CMAKE_SOURCE_DIR}/src/*.cpp"
     "${CMAKE_SOURCE_DIR}/src/motorsim/*.cpp")

add_executable(motor_sim ${SRC_FILES})
target_include_directories(motor_sim PRIVATE ${CMAKE_SOURCE_DIR}/external)

# Basic test target (placeholder)
enable_testing()
add_executable(smoke_test tests/smoke_test.cpp)
add_test(NAME smoke_run COMMAND smoke_test)

add_executable(analytic_wire_test
    tests/analytic_wire_test.cpp
    src/motorsim/solver.cpp
    src/motorsim/io_csv.cpp
    src/motorsim/ingest.cpp)
target_include_directories(analytic_wire_test PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/external)
add_test(NAME analytic_wire COMMAND analytic_wire_test)

add_executable(two_wire_cancel_test
    tests/two_wire_cancel_test.cpp
    src/motorsim/solver.cpp
    src/motorsim/ingest.cpp)
target_include_directories(two_wire_cancel_test PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/external)
add_test(NAME two_wire_cancel COMMAND two_wire_cancel_test)

add_executable(line_current_interface_test
    tests/line_current_interface_test.cpp
    src/motorsim/solver.cpp
    src/motorsim/ingest.cpp)
target_include_directories(line_current_interface_test PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/external)
add_test(NAME line_current_interface COMMAND line_current_interface_test)

add_executable(magnet_strip_test
    tests/magnet_strip_test.cpp
    src/motorsim/solver.cpp
    src/motorsim/ingest.cpp)
target_include_directories(magnet_strip_test PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/external)
add_test(NAME magnet_strip COMMAND magnet_strip_test)

add_executable(polygon_region_test
    tests/polygon_region_test.cpp
    src/motorsim/ingest.cpp)
target_include_directories(polygon_region_test PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/external)
add_test(NAME polygon_region COMMAND polygon_region_test)

add_executable(output_quantity_test
    tests/output_quantity_test.cpp
    src/motorsim/ingest.cpp
    src/motorsim/io_csv.cpp)
target_include_directories(output_quantity_test PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/external)
add_test(NAME output_quantity COMMAND output_quantity_test)

add_test(
    NAME scenario_python_emit
    COMMAND
        ${Python3_EXECUTABLE}
        ${CMAKE_SOURCE_DIR}/python/examples/two_wire_cancel.py
        --output
        ${CMAKE_BINARY_DIR}/two_wire_cancel_generated.json)

add_test(
    NAME scenario_json_diff
    COMMAND ${CMAKE_COMMAND} -E compare_files
            ${CMAKE_BINARY_DIR}/two_wire_cancel_generated.json
            ${CMAKE_SOURCE_DIR}/inputs/two_wire_cancel.json)
set_tests_properties(scenario_json_diff PROPERTIES DEPENDS scenario_python_emit)

add_test(
    NAME scenario_cli_solve
    COMMAND
        $<TARGET_FILE:motor_sim>
        --scenario
        ${CMAKE_BINARY_DIR}/two_wire_cancel_generated.json
        --solve
        --outputs
        none)
set_tests_properties(scenario_cli_solve PROPERTIES DEPENDS scenario_json_diff)

add_executable(solver_benchmark
    tools/solver_benchmark.cpp
    src/motorsim/solver.cpp)
target_include_directories(solver_benchmark PRIVATE ${CMAKE_SOURCE_DIR}/include)
