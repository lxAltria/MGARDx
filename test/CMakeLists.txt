add_executable (test_thomas_algorithm test_thomas_algorithm.cpp)

add_executable (test_decompose test_decompose.cpp)
target_include_directories(test_decompose PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_executable (test_real_data test_real_data.cpp)
target_include_directories(test_real_data PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_real_data PUBLIC "/Users/xin/github/SZ_refactory/include" "/Users/xin/utils/sz_master/include")
find_library(ZSTD_LIB zstd HINTS "/Users/xin/utils/sz_master/lib/")
find_library(SZ_LIB SZ3 HINTS "/Users/xin/github/SZ_refactory/install/lib")
target_link_libraries(test_real_data ${ZSTD_LIB} ${SZ})

add_executable (test_data_reorder test_data_reorder.cpp)
target_include_directories(test_data_reorder PUBLIC ${PROJECT_SOURCE_DIR}/include)
