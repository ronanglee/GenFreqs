cmake_minimum_required(VERSION 3.0.0)
project(vibrational VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

include_directories(vibrational declare.h)
include_directories(vibrational constnant.h)

add_executable(vibrational main.cc funcs.cc)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
