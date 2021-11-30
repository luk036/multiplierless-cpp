cpmaddpackage(
  NAME
  fmt
  GIT_TAG
  7.1.3
  GITHUB_REPOSITORY
  fmtlib/fmt
  OPTIONS
  "FMT_INSTALL YES" # create an installable target
)

cpmaddpackage("gh:xtensor-stack/xtl#0.6.23")
if(xtl_ADDED)
  message(STATUS "Found xtl: ${xtl_SOURCE_DIR}")
  include_directories(${xtl_SOURCE_DIR}/include)
endif(xtl_ADDED)

cpmaddpackage("gh:xtensor-stack/xtensor#0.22.0")
if(xtensor_ADDED)
  message(STATUS "Found xtensor: ${xtensor_SOURCE_DIR}")
  include_directories(${xtensor_SOURCE_DIR}/include)
endif(xtensor_ADDED)

cpmaddpackage(
  NAME
  EllAlgo
  GIT_TAG
  1.1
  GITHUB_REPOSITORY
  luk036/ellalgo-cpp
  OPTIONS
  "INSTALL_ONLY YES" # create an installable target
)

find_package(OpenBLAS QUIET)
if(OpenBLAS_FOUND)
  message(STATUS "Found OpenBLAS: ${OpenBLAS_LIBRARIES}")
  # target_include_directories(OpenBLAS::OpenBLAS SYSTEM INTERFACE
  # ${OpenBLAS_INCLUDE_DIRS})
  include_directories(${OpenBLAS_INCLUDE_DIRS})
endif(OpenBLAS_FOUND)

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
  message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
  include_directories(${LAPACK_INCLUDE_DIRS})
endif(LAPACK_FOUND)

find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
  include_directories(${BLAS_INCLUDE_DIRS})
endif(BLAS_FOUND)

cpmaddpackage("gh:xtensor-stack/xtensor-blas#0.18.0")
if(xtensor-blas_ADDED)
  message(STATUS "Found xtensor-blas: ${xtensor-blas_SOURCE_DIR}")
  include_directories(${xtensor-blas_SOURCE_DIR}/include)
endif(xtensor-blas_ADDED)
# remember to turn off the warnings

if(WIN32)
  add_definitions(-DXTENSOR_USE_FLENS_BLAS)
endif()

# .. fftw
if(MSVC)
  # no long double component, since in the Windows conda-forge build it is not
  # available and the "official" prebuilt long double library can only be used
  # from MinGW
  find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB)
  add_definitions(-DFFTW_NO_LONGDOUBLE)
  # add_compile_definitions(_USE_MATH_DEFINES)
else()
  find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB LONGDOUBLE_LIB)
endif(MSVC)

if(FFTW_FOUND)
  set(LIBS ${LIBS} ${FFTW_LIBRARIES})
  message(STATUS "Found FFTW: ${FFTW_LIBRARIES}")
  include_directories(${FFTW_INCLUDE_DIRS})
endif(FFTW_FOUND)

cpmaddpackage("gh:xtensor-stack/xtensor-fftw#0.2.5")
if(xtensor-fftw_ADDED)
  message(STATUS "Found xtensor-fftw: ${xtensor-fftw_SOURCE_DIR}")
  include_directories(${xtensor-fftw_SOURCE_DIR}/include)
endif(xtensor-fftw_ADDED)

set(SPECIFIC_LIBS EllAlgo::EllAlgo ${OpenBLAS_LIBRARIES} ${LAPACK_LIBRARIES}
                  ${BLAS_LIBRARIES} ${FFTW_LIBRARIES} fmt::fmt)
