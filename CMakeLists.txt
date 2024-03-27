cmake_minimum_required(VERSION 3.10)
project(clipperplus VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 14)
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose build type." FORCE)
endif()


###############################################################################
# Options
###############################################################################

option(OPTION_DEBUG_FLAG "Enable debugging output" OFF)
option(OPTION_DEBUG_TIMING_FLAG "Enable timing report" ON)
option(OPTION_DEBUG_OPTIM_FLAG "Enable debugging clipper optimiziation" ON)

option(OPTION_BUILD_PMC_HEU "Enable debugging output" ON)
option(OPTION_BUILD_BINDINGS_MATLAB "Build MATLAB bindings" ON)
# option(OPTION_BUILD_BINDINGS_PYTHON "Build Python bindings" ON)
option(BUILD_TESTS "Build the tests" ON)

###############################################################################
# Google Test Integration
###############################################################################

# Download and unpack googletest at configure time
if(BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.10.0
    )

    FetchContent_MakeAvailable(googletest)

    enable_testing()
    add_subdirectory(test/cpp_tests/)
endif()


###############################################################################
# Dependencies
###############################################################################

include(FetchContent)

find_package(Eigen3 REQUIRED)
message(STATUS "Eigen Version: ${EIGEN3_VERSION_STRING} ${EIGEN3_VERSION}")

find_package(OpenMP REQUIRED)

###############################################################################
# Targets
###############################################################################

set(PMC_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_heu.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_cores.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_graph.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_clique_utils.cpp )

set(CLIPPER_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/external/clipper/clipper.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/clipper/utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/clipper/dsd.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/clipper/invariants/euclidean_distance.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/clipper/invariants/pointnormal_distance.cpp )
        
set(CLLIPPER_PLUS_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/clipperplus_clique.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/clique_corenumber.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/clique_optimization.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/utils.cpp )

# source files
add_library(clipperplus SHARED)
target_sources(clipperplus PRIVATE 
        ${PMC_SOURCE_FILES}
        ${CLIPPER_SOURCE_FILES}        
        ${CLLIPPER_PLUS_SOURCE_FILES} )

# headers:
target_include_directories(clipperplus PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

# link libraries:
target_link_libraries(clipperplus PUBLIC
        OpenMP::OpenMP_CXX
        Eigen3::Eigen)

# debug flag  
if(OPTION_DEBUG_FLAG)
        target_compile_definitions(clipperplus PRIVATE DEBUG)
endif()
if(OPTION_DEBUG_OPTIM_FLAG)
        target_compile_definitions(clipperplus PRIVATE DEBUG_OPTIM)
endif()
if(OPTION_DEBUG_TIMING_FLAG)
        target_compile_definitions(clipperplus PRIVATE DEBUG_TIMING)
endif()

# version
set_target_properties(clipperplus PROPERTIES
  VERSION ${PROJECT_VERSION})
target_compile_definitions(clipperplus PUBLIC CLIPPERPLUS_VERSION="${PROJECT_VERSION}")


###############################################################################
# Extras
###############################################################################

if(OPTION_BUILD_BINDINGS_MATLAB)
    message(STATUS "Attempting to build MATLAB bindings.")
    add_subdirectory(bindings/matlab)
endif()

# if(BUILD_BINDINGS_PYTHON)
#     message(STATUS "Building Python bindings.")
#     add_subdirectory(bindings/python)
# endif()

# building PMC heuristic algorithm
if(OPTION_BUILD_PMC_HEU)
        set(PMC_HUE_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_heuristic.cpp )

        add_library(pmcheu SHARED)
        target_sources(pmcheu PRIVATE 
                ${PMC_SOURCE_FILES}
                ${PMC_HUE_SOURCE_FILES} )

        # headers:
        target_include_directories(pmcheu PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                $<INSTALL_INTERFACE:include>)

        # link libraries:
        target_link_libraries(pmcheu PUBLIC
                OpenMP::OpenMP_CXX
                Eigen3::Eigen)

        # debug flag  
        if(OPTION_DEBUG_FLAG)
                target_compile_definitions(pmcheu PRIVATE DEBUG)
        endif()
endif()



###############################################################################
# Tests
###############################################################################

# test for clipper_plus
add_executable(clipperplus_test src/clipperplus_test.cpp) 
target_link_libraries(clipperplus_test 
                      clipperplus  
                      pmcheu                    
                      Eigen3::Eigen)

# if(BUILD_TESTS)
#     enable_testing()
#     add_subdirectory(test)
# endif()


###############################################################################
# Install
###############################################################################
# include(GNUInstallDirs)

# install(TARGETS pmc
#         EXPORT pmc-export
#         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#         )

# export(
#         EXPORT pmc-export
#         FILE ${CMAKE_CURRENT_BINARY_DIR}/pmcTargets.cmake
#         NAMESPACE pmc::
# )

# install(EXPORT pmc-export
#         DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pmc
#         NAMESPACE pmc::
#         FILE pmcTargets.cmake
#         )

# install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# install(FILES cmake/pmcConfig.cmake
#         DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pmc
#         )
