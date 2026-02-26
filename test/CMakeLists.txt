# ---- Options ----

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
# option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)

# ---- Dependencies ----

CPMAddPackage("gh:doctest/doctest#v2.4.11")
CPMAddPackage("gh:TheLartians/Format.cmake@1.7.3")

# Add RapidCheck for property-based testing
CPMAddPackage(
  NAME rapidcheck
  GITHUB_REPOSITORY emil-e/rapidcheck
  GIT_TAG master
  OPTIONS "RC_BUILD_TESTS OFF" "RC_BUILD_EXAMPLES OFF" "RC_INSTALL OFF"
)

# if(TEST_INSTALLED_VERSION) find_package(${PROJECT_NAME} REQUIRED) else() CPMAddPackage(NAME
# ${PROJECT_NAME} SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..) endif()

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
add_executable(${PROJECT_NAME}Tests ${sources})
target_link_libraries(
  ${PROJECT_NAME}Tests doctest::doctest ${PROJECT_NAME}::${PROJECT_NAME} ${SPECIFIC_LIBS}
)
set_target_properties(${PROJECT_NAME}Tests PROPERTIES CXX_STANDARD 17)

# Link RapidCheck if available
if(rapidcheck_ADDED)
  target_compile_definitions(${PROJECT_NAME}Tests PRIVATE RAPIDCHECK_H)
  if(TARGET rapidcheck::rapidcheck)
    target_link_libraries(${PROJECT_NAME}Tests rapidcheck::rapidcheck)
  else()
    target_include_directories(${PROJECT_NAME}Tests PRIVATE ${rapidcheck_SOURCE_DIR}/include)
    if(TARGET rapidcheck)
      target_link_libraries(${PROJECT_NAME}Tests rapidcheck)
    endif()
  endif()
endif()

# enable compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # target_compile_options(${PROJECT_NAME} PUBLIC -Wall -Wpedantic -Wextra -Werror)
  target_compile_options(
    ${PROJECT_NAME} PUBLIC -Wall -Wextra -Wconversion -Wno-parentheses
                           -Wno-unused-but-set-parameter
  )
elseif(MSVC)
  target_compile_options(
    ${PROJECT_NAME}
    PUBLIC /std:c++latest
           /await
           /W4
           /wd4819
           /wd4146
           /wd4127
           /wd4244
  )
  target_compile_definitions(${PROJECT_NAME}Tests PUBLIC DOCTEST_CONFIG_USE_STD_HEADERS)
endif()

# ---- Add ${PROJECT_NAME}Tests ----

enable_testing()

# Add a single test that runs all tests (both doctest and RapidCheck)
# Note: We don't use doctest_discover_tests because it can conflict with RapidCheck
string(TOLOWER ${PROJECT_NAME} PROJECT_FILE_NAME)
add_test(NAME ${PROJECT_FILE_NAME}Tests COMMAND ${PROJECT_NAME}Tests)

# ---- code coverage ----

if(ENABLE_TEST_COVERAGE)
  target_compile_options(${PROJECT_NAME} PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
  target_link_options(${PROJECT_NAME} PUBLIC -fprofile-arcs -ftest-coverage)
endif()
