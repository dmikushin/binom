cmake_minimum_required(VERSION 3.1)

project(binom CXX)

if (NOT TARGET gtest)
add_subdirectory(ThirdParty/googletest EXCLUDE_FROM_ALL)
endif()

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(${PROJECT_NAME}_test "src/${PROJECT_NAME}_test.cpp")
target_link_libraries(${PROJECT_NAME}_test ${PROJECT_NAME} gtest)

