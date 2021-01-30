cmake_minimum_required(VERSION 3.12)

project(pazi1-2)
add_executable(${PROJECT_NAME} "main.c" "curve.c" "curve.h")

find_package(PkgConfig REQUIRED)
pkg_check_modules(math REQUIRED IMPORTED_TARGET libtommath)
target_link_libraries(${PROJECT_NAME} PkgConfig::math)
