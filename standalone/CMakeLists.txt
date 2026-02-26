# ---- Dependencies ----

CPMAddPackage(
  GITHUB_REPOSITORY jarro2783/cxxopts
  VERSION 3.2.1
  OPTIONS "CXXOPTS_BUILD_EXAMPLES NO" "CXXOPTS_BUILD_TESTS NO" "CXXOPTS_ENABLE_INSTALL YES"
)

# ---- Create standalone executable ----

add_executable(${PROJECT_NAME}Standalone source/main.cpp)

set_target_properties(
  ${PROJECT_NAME}Standalone PROPERTIES CXX_STANDARD 17 OUTPUT_NAME ${PROJECT_NAME}
)

target_link_libraries(
  ${PROJECT_NAME}Standalone ${PROJECT_NAME}::${PROJECT_NAME} cxxopts::cxxopts ${SPECIFIC_LIBS}
)

# ---- Create logger example executable ----

add_executable(${PROJECT_NAME}LoggerExample source/logger_example.cpp)

set_target_properties(
  ${PROJECT_NAME}LoggerExample PROPERTIES CXX_STANDARD 17 OUTPUT_NAME logger_example
)

target_link_libraries(
  ${PROJECT_NAME}LoggerExample ${PROJECT_NAME}::${PROJECT_NAME} ${SPECIFIC_LIBS}
)
