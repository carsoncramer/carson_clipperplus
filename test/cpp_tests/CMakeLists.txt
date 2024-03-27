cmake_minimum_required(VERSION 3.10)
project(ClipperPlusTests)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.10.0
)
FetchContent_MakeAvailable(googletest)

# Enable testing functionalities
enable_testing()

# Specify the test source files
set(TEST_SOURCES
    unit/test_clipperplus_clique.cpp
    unit/test_clique_corenumber.cpp
    unit/test_clique_optimization.cpp
    unit/test_utils.cpp
)

# Create a test executable
add_executable(tests ${TEST_SOURCES})

# Link the test executable with GoogleTest and the main project library
target_link_libraries(tests gtest_main gmock clipperplus)

find_package(Threads REQUIRED)
target_link_libraries(tests Threads::Threads)

target_include_directories(tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../include/clipperplus)

# Discover tests
include(GoogleTest)
gtest_discover_tests(tests)
