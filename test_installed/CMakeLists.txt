cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# ---- Include guards ----

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

project(Multiplierless LANGUAGES CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# --- Import tools ----

include(../cmake/tools.cmake)

# ---- Dependencies (all installed by the parent build) ----

include(../cmake/CPM.cmake)

# Find dependencies first — MultiplierlessConfig.cmake references fmt::fmt etc.
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
find_package(EllAlgo REQUIRED)
find_package(Ginger REQUIRED)
find_package(Csd REQUIRED)
find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB LONGDOUBLE_LIB)
find_package(Threads REQUIRED)
CPMAddPackage("gh:nlohmann/json@3.12.0")
find_package(${PROJECT_NAME} REQUIRED)

set(SPECIFIC_LIBS Ginger::Ginger EllAlgo::EllAlgo Csd::Csd
    ${FFTW_LIBRARIES} Threads::Threads fmt::fmt spdlog::spdlog)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../standalone ${CMAKE_BINARY_DIR}/test_installed)
