cmake_minimum_required(VERSION 2.6)
project(sigref C CXX)

set(CMAKE_C_FLAGS "-g -O3 -Wall -Wextra -Werror -std=gnu11")
set(CMAKE_CXX_FLAGS "-g -O3 -Wall -Wextra -Werror -std=gnu++11")

include_directories(.)

add_executable(sigrefmc sigref.cpp sigref.h bisim_ctmc.cpp bisim_lts.cpp bisim_imc.cpp getrss.h getrss.c blocks.h blocks.c inert.h inert.c parse_bdd.hpp parse_bdd.cpp parse_xml.hpp parse_xml.cpp refine_sl.c systems.hpp sigref_util.h sigref_util.cpp)
target_link_libraries(sigrefmc sylvan tinyxml)

add_executable(sigrefmc_ht sigref.cpp sigref.h bisim_ctmc.cpp bisim_lts.cpp bisim_imc.cpp getrss.h getrss.c blocks.h blocks.c inert.h inert.c parse_bdd.hpp parse_bdd.cpp parse_xml.hpp parse_xml.cpp refine_ht.c systems.hpp sigref_util.h sigref_util.cpp)
target_link_libraries(sigrefmc_ht sylvan tinyxml)

include(CheckIncludeFiles)
check_include_files("gperftools/profiler.h" HAVE_PROFILER)

find_path(HAVE_TINYXML_H tinyxml.h)

find_package(Boost REQUIRED)

if(NOT HAVE_TINYXML_H)
    message(FATAL_ERROR "tinyxml not found")
endif()

if(HAVE_PROFILER)
    set_target_properties(sigrefmc PROPERTIES COMPILE_DEFINITIONS "HAVE_PROFILER")
    target_link_libraries(sigrefmc profiler)

    set_target_properties(sigrefmc_ht PROPERTIES COMPILE_DEFINITIONS "HAVE_PROFILER")
    target_link_libraries(sigrefmc_ht profiler)
endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # add argp library for OSX
    target_link_libraries(sigrefmc argp)
    target_link_libraries(sigrefmc_ht argp)
endif()
