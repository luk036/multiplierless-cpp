cmake_minimum_required(VERSION 3.14...3.22)

project(
  Multiplierless
  VERSION 1.0.4
  LANGUAGES CXX
)

# ---- Include guards ----

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

# ---------------------------------------------------------------------------
# Compiler warnings (as errors)
# ---------------------------------------------------------------------------

function(target_enable_warnings TARGET)
  if(MSVC)
    target_compile_options(${TARGET} PRIVATE /W4 /WX /utf-8 /wd4819)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(${TARGET} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
endfunction()

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(MULTIPLIERLESS_BUILD_TESTS "Build tests" ON)
option(MULTIPLIERLESS_BUILD_STANDALONE "Build standalone executable" ON)
option(MULTIPLIERLESS_BUILD_DOCS "Build Doxygen documentation" OFF)
option(MULTIPLIERLESS_ENABLE_CLANG_TIDY "Enable clang-tidy analysis" OFF)
option(MULTIPLIERLESS_ENABLE_COVERAGE "Enable coverage reporting" OFF)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------

# Download CPM.cmake for reproducible dependency management
# see https://github.com/cpm-cmake/CPM.cmake for more info
set(CPM_DOWNLOAD_VERSION 0.38.5)
set(CPM_HASH_SUM "192aa0ccdc57dfe75bd9e4b176bf7fb5692fd2b3e3f7b09c74856fc39572b31c")
set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
if(NOT (EXISTS ${CPM_DOWNLOAD_LOCATION}))
  message(STATUS "Downloading CPM.cmake to ${CPM_DOWNLOAD_LOCATION}")
  file(
    DOWNLOAD https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake
             ${CPM_DOWNLOAD_LOCATION}
    EXPECTED_HASH SHA256=${CPM_HASH_SUM}
  )
endif()
include(${CPM_DOWNLOAD_LOCATION})

# PackageProject.cmake will be used to make our target installable
CPMAddPackage("gh:TheLartians/PackageProject.cmake@1.8.0")

# FindFFTW.cmake lives in cmake/; expose it to find_package(FFTW).
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# fmt is a transitive dependency of spdlog; prefer a system-installed fmt when available.
find_package(fmt CONFIG QUIET)
if(fmt_FOUND)
  message(STATUS "Found system fmt: ${fmt_DIR}")
  set(CPM_fmt_ADDED YES)
else()
  CPMAddPackage(
    NAME fmt
    GIT_TAG 12.1.0
    GITHUB_REPOSITORY fmtlib/fmt
    OPTIONS "FMT_INSTALL YES" # create an installable target
  )
endif()

# spdlog is used by the compiled logger source; prefer a system-installed spdlog when available.
# SPDLOG_FMT_EXTERNAL is always set so spdlog uses the external fmt and avoids duplicate symbols.
find_package(spdlog CONFIG QUIET)
if(spdlog_FOUND)
  message(STATUS "Found system spdlog: ${spdlog_DIR}")
  set(CPM_spdlog_ADDED YES)
else()
  CPMAddPackage(
    NAME spdlog
    GIT_TAG v1.17.0
    GITHUB_REPOSITORY gabime/spdlog
    OPTIONS "SPDLOG_INSTALL YES" "SPDLOG_FMT_EXTERNAL YES"
  )
endif()

# EllAlgo provides the oracle/arr types used across the public headers.
CPMAddPackage(
  NAME EllAlgo
  GIT_TAG v1.6.9
  GITHUB_REPOSITORY luk036/ellalgo-cpp
  OPTIONS "ELLALGO_BUILD_TESTS OFF" "ELLALGO_BUILD_STANDALONE OFF"
)

# ellalgo v1.6.9 compiles spdlog headers in its library but does not link fmt. With
# SPDLOG_FMT_EXTERNAL (set on the shared spdlog above), its logger needs fmt's include dir.
if(TARGET EllAlgo)
  get_target_property(FMT_INCLUDE_DIRS fmt::fmt INTERFACE_INCLUDE_DIRECTORIES)
  target_include_directories(EllAlgo PRIVATE ${FMT_INCLUDE_DIRS})
endif()

# Csd provides canonical-signed-digit conversions.
CPMAddPackage(
  NAME Csd
  GIT_TAG v1.1.5
  GITHUB_REPOSITORY luk036/csd-cpp
  OPTIONS "CSD_BUILD_TESTS OFF" "CSD_BUILD_STANDALONE OFF"
)

# Ginger provides polynomial root-finding; installed so consumers can find its config.
CPMAddPackage(
  NAME Ginger
  GIT_TAG v1.1.6
  GITHUB_REPOSITORY luk036/ginger-cpp
  OPTIONS "GINGER_BUILD_TESTS OFF" "GINGER_BUILD_STANDALONE OFF"
)

# nlohmann_json is used by the FIR design CLI; installed so consumers can find its config.
CPMAddPackage(
  NAME nlohmann_json
  GIT_TAG v3.12.0
  GITHUB_REPOSITORY nlohmann/json
  OPTIONS "JSON_BuildTests OFF" "JSON_Install ON"
)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# FFTW for spectral_fact_fft. On Windows with conda, help FindFFTW.cmake locate the library.
if(WIN32
   AND DEFINED ENV{CONDA_PREFIX}
   AND NOT FFTW_ROOT
   AND NOT DEFINED ENV{FFTWDIR}
)
  set(FFTW_ROOT "$ENV{CONDA_PREFIX}/Library")
elseif(DEFINED ENV{FFTWDIR} AND NOT FFTW_ROOT)
  set(FFTW_ROOT "$ENV{FFTWDIR}")
endif()

if(MSVC)
  find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB)
  add_definitions(-DFFTW_NO_LONGDOUBLE)
else()
  find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB LONGDOUBLE_LIB)
endif()

set(SPECIFIC_LIBS
    Ginger::Ginger
    EllAlgo::EllAlgo
    Csd::Csd
    ${FFTW_LIBRARIES}
    Threads::Threads
    fmt::fmt
    spdlog::spdlog
)

# ---------------------------------------------------------------------------
# Library target
# ---------------------------------------------------------------------------

# Note: globbing sources is considered bad practice as CMake's generators may not detect new files
# automatically. Keep that in mind when changing files, or explicitly mention them here.
file(GLOB_RECURSE headers CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp"
)
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp")

add_library(${PROJECT_NAME} ${headers} ${sources})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 20)

# Enable interprocedural optimization (LTO) for better inlining and code generation
set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

# being a cross-platform target, we enforce standards conformance on MSVC
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")

# Public headers include <ellalgo/...> and <fftw3.h>, so these dependencies are public.
target_link_libraries(${PROJECT_NAME} PUBLIC ${SPECIFIC_LIBS})

# Any TU using spdlog must compile with spdlog's own interface definitions (SPDLOG_COMPILED_LIB,
# SPDLOG_FMT_EXTERNAL) or the header-only and compiled registries diverge.
get_target_property(SPDLOG_COMPILE_DEFS spdlog INTERFACE_COMPILE_DEFINITIONS)
if(SPDLOG_COMPILE_DEFS)
  target_compile_definitions(${PROJECT_NAME} PRIVATE ${SPDLOG_COMPILE_DEFS})
endif()

# FFTW is found via include dirs, not a CMake target.
if(FFTW_INCLUDE_DIRS)
  target_include_directories(${PROJECT_NAME} PUBLIC ${FFTW_INCLUDE_DIRS})
endif()

target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

target_enable_warnings(${PROJECT_NAME})

# ---- Create an installable target ----
# this allows users to install and find the library via `find_package()`.

# the location where the project's version header will be placed should match the project's regular
# header paths
string(TOLOWER ${PROJECT_NAME}/version.h VERSION_HEADER_LOCATION)

packageProject(
  NAME ${PROJECT_NAME}
  VERSION ${PROJECT_VERSION}
  NAMESPACE ${PROJECT_NAME}
  BINARY_DIR ${PROJECT_BINARY_DIR}
  INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
  INCLUDE_DESTINATION include/${PROJECT_NAME}-${PROJECT_VERSION}
  VERSION_HEADER "${VERSION_HEADER_LOCATION}"
  COMPATIBILITY SameMajorVersion
  DEPENDENCIES "fmt 12.1.0" "EllAlgo 1.6.9" "Csd 1.1.5" "Ginger 1.1.6"
)

# ---------------------------------------------------------------------------
# Standalone executables
# ---------------------------------------------------------------------------

if(MULTIPLIERLESS_BUILD_STANDALONE)
  # cxxopts is a header-only CLI parser used by the standalone example
  CPMAddPackage(
    GITHUB_REPOSITORY jarro2783/cxxopts
    VERSION 3.2.1
    OPTIONS "CXXOPTS_BUILD_EXAMPLES NO" "CXXOPTS_BUILD_TESTS NO" "CXXOPTS_ENABLE_INSTALL NO"
  )
  # keep third-party headers out of our -Werror scope
  get_target_property(CXXOPTS_INCLUDE_DIRS cxxopts::cxxopts INTERFACE_INCLUDE_DIRECTORIES)

  add_executable(MultiplierlessStandalone standalone/source/main.cpp)
  set_target_properties(MultiplierlessStandalone PROPERTIES OUTPUT_NAME "Multiplierless"
                                                       CXX_STANDARD 20
  )
  target_include_directories(MultiplierlessStandalone SYSTEM PRIVATE ${CXXOPTS_INCLUDE_DIRS})
  target_link_libraries(MultiplierlessStandalone PRIVATE Multiplierless::Multiplierless
                                                  cxxopts::cxxopts ${SPECIFIC_LIBS}
  )
  target_enable_warnings(MultiplierlessStandalone)

  add_executable(MultiplierlessLoggerExample standalone/source/logger_example.cpp)
  set_target_properties(MultiplierlessLoggerExample PROPERTIES OUTPUT_NAME "logger_example"
                                                             CXX_STANDARD 20
  )
  target_link_libraries(MultiplierlessLoggerExample PRIVATE Multiplierless::Multiplierless
                                                      ${SPECIFIC_LIBS}
  )
  target_enable_warnings(MultiplierlessLoggerExample)

  add_executable(MultiplierlessFirDesign standalone/source/fir_design.cpp)
  set_target_properties(MultiplierlessFirDesign PROPERTIES OUTPUT_NAME "FirDesign" CXX_STANDARD 20)
  target_link_libraries(MultiplierlessFirDesign PRIVATE Multiplierless::Multiplierless
                                                  nlohmann_json::nlohmann_json ${SPECIFIC_LIBS}
  )
  target_enable_warnings(MultiplierlessFirDesign)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

enable_testing()

if(MULTIPLIERLESS_BUILD_TESTS)
  # Try system-installed doctest first (Ubuntu: doctest-dev, macOS: brew install doctest, Termux:
  # doctest)
  find_package(doctest CONFIG QUIET)

  if(NOT doctest_FOUND)
    CPMAddPackage(
      NAME doctest
      GITHUB_REPOSITORY doctest/doctest
      VERSION 2.5.2
      OPTIONS "DOCTEST_NO_INSTALL ON" # keep the installed prefix minimal
    )
  endif()

  # Code formatting targets (format / fix-format) via Format.cmake
  CPMAddPackage("gh:TheLartians/Format.cmake@1.7.3")

  # Property-based testing (tests only)
  CPMAddPackage(
    NAME rapidcheck
    GITHUB_REPOSITORY emil-e/rapidcheck
    GIT_TAG master
    OPTIONS "RC_BUILD_TESTS OFF" "RC_BUILD_EXAMPLES OFF" "RC_INSTALL OFF"
  )

  file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "test/source/*.cpp")
  add_executable(MultiplierlessTests ${TEST_SOURCES})
  target_link_libraries(MultiplierlessTests PRIVATE doctest::doctest Multiplierless::Multiplierless
                                               ${SPECIFIC_LIBS} rapidcheck
  )
  target_include_directories(MultiplierlessTests PRIVATE ${rapidcheck_SOURCE_DIR}/include)
  target_compile_definitions(MultiplierlessTests PRIVATE RAPIDCHECK_H)
  set_target_properties(MultiplierlessTests PROPERTIES CXX_STANDARD 20)
  target_enable_warnings(MultiplierlessTests)
  add_test(NAME MultiplierlessAllTests COMMAND MultiplierlessTests)
endif()

# ---------------------------------------------------------------------------
# clang-tidy
# ---------------------------------------------------------------------------

if(MULTIPLIERLESS_ENABLE_CLANG_TIDY)
  find_program(CLANG_TIDY clang-tidy)
  if(CLANG_TIDY)
    add_custom_target(
      clang-tidy
      COMMAND ${CLANG_TIDY} --extra-arg=-std=c++20
              --extra-arg=-I${CMAKE_CURRENT_SOURCE_DIR}/include
              ${CMAKE_CURRENT_SOURCE_DIR}/include/multiplierless/*.hpp
      COMMENT "Running clang-tidy static analysis"
    )
  else()
    message(WARNING "clang-tidy not found - analysis target unavailable")
  endif()
endif()

# ---------------------------------------------------------------------------
# Doxygen documentation
# ---------------------------------------------------------------------------

if(MULTIPLIERLESS_BUILD_DOCS)
  find_package(Doxygen QUIET)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_PROJECT_NAME ${PROJECT_NAME})
    set(DOXYGEN_PROJECT_VERSION ${PROJECT_VERSION})
    set(DOXYGEN_PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doxygen")

    configure_file(documentation/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    add_custom_target(
      GenerateDocs
      COMMAND ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
      COMMAND echo "Docs written to: ${DOXYGEN_OUTPUT_DIRECTORY}/html"
      COMMENT "Generating Doxygen HTML documentation in ${DOXYGEN_OUTPUT_DIRECTORY}/html"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
  else()
    add_custom_target(
      GenerateDocs COMMAND ${CMAKE_COMMAND} -E echo
                           "Doxygen not found - skipping documentation generation"
    )
  endif()
endif()

# ---------------------------------------------------------------------------
# Coverage report target (GCC/Clang only - requires gcovr)
# ---------------------------------------------------------------------------

if(MULTIPLIERLESS_ENABLE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  foreach(COVERAGE_TARGET Multiplierless MultiplierlessTests MultiplierlessStandalone)
    if(TARGET ${COVERAGE_TARGET})
      target_compile_options(${COVERAGE_TARGET} PRIVATE -O0 -g -fprofile-arcs -ftest-coverage)
    endif()
  endforeach()
  foreach(COVERAGE_EXEC MultiplierlessTests MultiplierlessStandalone)
    if(TARGET ${COVERAGE_EXEC})
      target_link_options(${COVERAGE_EXEC} PRIVATE -fprofile-arcs -ftest-coverage)
    endif()
  endforeach()

  find_program(GCOVR_PATH gcovr)
  if(GCOVR_PATH)
    set(extra_gcovr_args "")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang")
      find_program(LLVM_COV_PATH llvm-cov)
      if(LLVM_COV_PATH)
        list(APPEND extra_gcovr_args --llvm-cov-binary ${LLVM_COV_PATH})
      endif()
    endif()
    add_custom_target(
      coverage
      COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/coverage
      COMMAND ${GCOVR_PATH} --root ${PROJECT_SOURCE_DIR} --html-details
              ${CMAKE_BINARY_DIR}/coverage/index.html --exclude test/ --exclude build/
              ${extra_gcovr_args}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Generating coverage report in ${CMAKE_BINARY_DIR}/coverage/index.html"
    )
  else()
    message(WARNING "gcovr not found - coverage target unavailable. Install with: pip install gcovr")
  endif()
endif()
