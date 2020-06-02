add_executable (test_decompose test_decompose.cpp)
target_include_directories(test_decompose PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_decompose PUBLIC ${SZ3_INC} ${ZSTD_INC} ${SZ_CPP_INC})
target_link_libraries(test_decompose ${ZSTD_LIB} ${SZ3_LIB} ${SZ_CPP_LIB})

add_executable (test_compress test_compress.cpp)
target_include_directories(test_compress PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_compress PUBLIC ${SZ3_INC} ${ZSTD_INC} ${SZ_CPP_INC})
target_link_libraries(test_compress ${ZSTD_LIB} ${SZ3_LIB} ${SZ_CPP_LIB})

add_executable (test_decompress test_decompress.cpp)
target_include_directories(test_decompress PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_decompress PUBLIC ${SZ3_INC} ${ZSTD_INC} ${SZ_CPP_INC})
target_link_libraries(test_decompress ${ZSTD_LIB} ${SZ3_LIB} ${SZ_CPP_LIB})

add_executable (test_refactor test_refactor.cpp)
target_include_directories(test_refactor PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_refactor PUBLIC ${SZ3_INC} ${ZSTD_INC} ${SZ_CPP_INC})
target_link_libraries(test_refactor ${ZSTD_LIB} ${SZ3_LIB} ${SZ_CPP_LIB})
