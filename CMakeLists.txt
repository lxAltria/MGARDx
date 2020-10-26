cmake_minimum_required (VERSION 3.0)
project(MGARDx)
project(First VERSION 1.0.0.0)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_FLAGS "-g -O3")

message(STATUS "${CMAKE_INSTALL_PREFIX}")
message(STATUS "${CMAKE_INSTALL_PREFIX}")
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE include)
install(TARGETS ${PROJECT_NAME}
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        INCLUDES DESTINATION ${CMAKE_INSTALL_PREFIX}/include
        )

add_subdirectory (test)
