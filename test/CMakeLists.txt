add_executable (test_thomas_algorithm test_thomas_algorithm.cpp)

add_executable (test_decompose test_decompose.cpp)
target_include_directories(test_decompose PUBLIC ${PROJECT_SOURCE_DIR}/include)
