cmake_minimum_required(VERSION 3.25)
project(anisotropy CXX)

# Prohibit in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(
    FATAL_ERROR
    "In-source build prohibited. \n"
    "Use a separate build directory with:\n  cmake -B build"
  )
endif()

# Set the build type to Debug if unspecified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type" FORCE)
endif()

# Export compile commands for IDEs
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Vendor libassert for some cool assertions
include(FetchContent)
FetchContent_Declare(
  libassert
  GIT_REPOSITORY https://github.com/jeremy-rifkin/libassert.git
  GIT_TAG        v2.2.1
)
FetchContent_MakeAvailable(libassert)

# Executable and link libraries
add_executable(anisotropy)
target_link_libraries(anisotropy libassert::assert)
target_compile_features(anisotropy PRIVATE cxx_std_20)
target_sources(anisotropy PRIVATE main.cpp)
