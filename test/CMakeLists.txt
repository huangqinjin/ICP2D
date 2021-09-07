cmake_minimum_required(VERSION 3.21)
project(${PACKAGE_NAME}-test)

file(DOWNLOAD "https://github.com/catchorg/Catch2/releases/download/v2.13.7/catch.hpp"
        "${CMAKE_CURRENT_LIST_DIR}/catch.hpp"
        EXPECTED_HASH SHA512=61da96855603c21b5ff3d354fc96ec5f28a57169460495fc9b584e3d91d0eee444231d2d97a5ad6541f4c33caf9f7d646ec958b9621aea352cc117ca1992ddf8)

file(GLOB SRC_FILES *.cpp)

add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC ${PACKAGE_NAME}::${PACKAGE_NAME})
