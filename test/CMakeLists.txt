find_library(ZSTD_LIB zstd HINTS "/Users/xin/utils/sz_master/lib/")
find_library(SZ_LIB SZ3 HINTS "/Users/xin/github/SZ_refactory/install/lib")

add_executable (test_real_data test_real_data.cpp)
target_include_directories(test_real_data PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_real_data PUBLIC "/Users/xin/github/SZ_refactory/include" "/Users/xin/utils/sz_master/include")
target_link_libraries(test_real_data ${ZSTD_LIB} ${SZ})

add_executable (test_decompose test_decompose.cpp)
target_include_directories(test_decompose PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_decompose PUBLIC "/Users/xin/github/SZ_refactory/include" "/Users/xin/utils/sz_master/include")
target_link_libraries(test_decompose ${ZSTD_LIB} ${SZ})

add_executable (test_compress test_compress.cpp)
target_include_directories(test_compress PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_compress PUBLIC "/Users/xin/github/SZ_refactory/include" "/Users/xin/utils/sz_master/include")
target_link_libraries(test_compress ${ZSTD_LIB} ${SZ})

add_executable (test_decompress test_decompress.cpp)
target_include_directories(test_decompress PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_decompress PUBLIC "/Users/xin/github/SZ_refactory/include" "/Users/xin/utils/sz_master/include")
target_link_libraries(test_decompress ${ZSTD_LIB} ${SZ})
