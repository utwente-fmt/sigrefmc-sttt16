cmake_minimum_required(VERSION 2.6)
project(sigrefmc C CXX)

include_directories(sylvan/src)
add_subdirectory(sylvan/src)

include_directories(src)
add_subdirectory(src)
