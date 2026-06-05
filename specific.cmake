set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

CPMAddPackage(
  NAME fmt
  GIT_TAG 12.1.0
  GITHUB_REPOSITORY fmtlib/fmt
  OPTIONS "FMT_INSTALL YES"
)

CPMAddPackage(
  NAME spdlog
  GIT_TAG v1.17.0
  GITHUB_REPOSITORY gabime/spdlog
  OPTIONS "SPDLOG_INSTALL YES"
)

CPMAddPackage(
  NAME EllAlgo
  GIT_TAG 1.6.5
  GITHUB_REPOSITORY luk036/ellalgo-cpp
  OPTIONS "INSTALL_ONLY YES"
)

CPMAddPackage(
  NAME Csd
  GIT_TAG 1.1.2
  GITHUB_REPOSITORY luk036/csd-cpp
  OPTIONS "INSTALL_ONLY YES"
)

CPMAddPackage(
  NAME nlohmann_json
  GIT_TAG v3.11.3
  GITHUB_REPOSITORY nlohmann/json
)

CPMAddPackage(
  NAME Ginger
  GIT_TAG 1.1.2
  GITHUB_REPOSITORY luk036/ginger-cpp
  OPTIONS "INSTALL_ONLY YES"
)

# FFTW for spectral_fact_fft
if(MSVC)
  find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB)
  add_definitions(-DFFTW_NO_LONGDOUBLE)
else()
  find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB LONGDOUBLE_LIB)
endif()

if(FFTW_FOUND)
  include_directories(${FFTW_INCLUDE_DIRS})
endif()

set(SPECIFIC_LIBS Ginger::Ginger EllAlgo::EllAlgo Csd::Csd ${FFTW_LIBRARIES} Threads::Threads fmt::fmt spdlog::spdlog)
