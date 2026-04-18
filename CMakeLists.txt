cmake_minimum_required(VERSION 3.25)

project(anisotropy CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(anisotropy)

target_compile_features(anisotropy PRIVATE cxx_std_20)

target_sources(anisotropy PRIVATE main.cpp)
