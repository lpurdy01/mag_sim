cmake_minimum_required(VERSION 3.16)
project(motor_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_compile_options(-Wall -Wextra -Wpedantic)

# Project include path
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source list (expand later)
file(GLOB SRC_FILES CONFIGURE_DEPENDS
     "${CMAKE_SOURCE_DIR}/src/*.cpp"
     "${CMAKE_SOURCE_DIR}/src/motorsim/*.cpp")

add_executable(motor_sim ${SRC_FILES})

# Basic test target (placeholder)
enable_testing()
add_executable(smoke_test tests/smoke_test.cpp)
add_test(NAME smoke_run COMMAND smoke_test)

add_executable(analytic_wire_test
    tests/analytic_wire_test.cpp
    src/motorsim/solver.cpp
    src/motorsim/io_csv.cpp)
target_include_directories(analytic_wire_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
add_test(NAME analytic_wire COMMAND analytic_wire_test)

add_executable(solver_benchmark
    tools/solver_benchmark.cpp
    src/motorsim/solver.cpp)
target_include_directories(solver_benchmark PRIVATE ${CMAKE_SOURCE_DIR}/include)
