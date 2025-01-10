cmake_minimum_required(VERSION 3.21)

# Project name and language
project(FortranProject LANGUAGES Fortran)
message(STATUS "building ${PROJECT_NAME}")

# Options
option(ENABLE_TESTING "Enable unit tests using Fortuno" ON)

# CMake modules for configuration and thirdparty projects
include(FetchContent)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# install settings
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "..." FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
include(GNUInstallDirs)


# Add the source directory
add_subdirectory(src)

if(ENABLE_TESTING)
    enable_testing()
    message(STATUS "Unit testing is enabled. Fetching Fortuno for testing.")
    include(fortuno)
    add_subdirectory(testings)
else()
    message(STATUS "Unit testing is disabled. Skipping Fortuno.")
endif(ENABLE_TESTING)

# Application
