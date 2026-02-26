cmake_minimum_required(VERSION 3.20)
project(eer2mrc LANGUAGES CXX)

# C++20 + MSVC extensions (to tolerate legacy code while keeping std::filesystem)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# vcpkg-provided libtiff
find_package(TIFF REQUIRED)

# Find EERReaderLib (already built by the parent library)
# Adjust the path if your EERReaderLib outputs elsewhere.
find_library(EERREADERLIB_LIBRARY
    NAMES eer2mrcLib
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib/Debug
)

if(NOT EERREADERLIB_LIBRARY)
    message(FATAL_ERROR "eer2mrcLib.lib not found in ../build/lib. Build eer2mrcLib first.")
endif()

add_executable(eer2mrc eer2mrc.cpp)

# Include EERReaderLib headers if they are in ../include
target_include_directories(eer2mrc PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_link_libraries(eer2mrc PRIVATE
    ${EERREADERLIB_LIBRARY}
    TIFF::TIFF
)