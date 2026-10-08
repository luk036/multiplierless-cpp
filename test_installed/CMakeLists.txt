cmake_minimum_required(VERSION 3.14)
project(MultiplierlessInstalledTest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# FindFFTW.cmake lives in the parent's cmake/ dir.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Find dependencies first - the installed config references fmt, EllAlgo, Csd, Ginger, and FFTW.
find_package(Threads REQUIRED)
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
find_package(EllAlgo REQUIRED)
find_package(Csd REQUIRED)
find_package(Ginger REQUIRED)
if(MSVC)
  find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB)
else()
  find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB LONGDOUBLE_LIB)
endif()
find_package(Multiplierless REQUIRED)

add_executable(multiplierless_installed_test main.cpp)
target_link_libraries(
  multiplierless_installed_test PRIVATE Multiplierless::Multiplierless ${FFTW_LIBRARIES}
)
target_include_directories(multiplierless_installed_test PRIVATE ${FFTW_INCLUDE_DIRS})

enable_testing()
add_test(NAME multiplierless_installed_test COMMAND multiplierless_installed_test)
