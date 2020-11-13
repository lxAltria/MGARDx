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

find_library(ZFP_LIB zfp HINTS "/Users/xin/github/zfp/build/lib")
set (ZFP_INC "/Users/xin/github/zfp/include")

add_executable (test_refactor test_refactor.cpp)
target_include_directories(test_refactor PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_refactor PUBLIC ${ZFP_INC} ${SZ3_INC} ${ZSTD_INC} ${SZ_CPP_INC})
target_link_libraries(test_refactor ${ZFP_LIB} ${ZSTD_LIB} ${SZ3_LIB} ${SZ_CPP_LIB})

add_executable (test_retrieval test_retrieval.cpp)
target_include_directories(test_retrieval PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_retrieval PUBLIC ${ZFP_INC} ${SZ3_INC} ${ZSTD_INC} ${SZ_CPP_INC})
target_link_libraries(test_retrieval ${ZFP_LIB} ${ZSTD_LIB} ${SZ3_LIB} ${SZ_CPP_LIB})

add_executable (test_bitplane_performance test_bitplane_performance.cpp)
target_include_directories(test_bitplane_performance PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_bitplane_performance PUBLIC ${ZFP_INC} ${SZ3_INC} ${ZSTD_INC} ${SZ_CPP_INC})
target_link_libraries(test_bitplane_performance ${ZFP_LIB} ${ZSTD_LIB} ${SZ3_LIB} ${SZ_CPP_LIB})

add_executable (test_quantization_factor test_quantization_factor.cpp)
target_include_directories(test_quantization_factor PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(test_quantization_factor PUBLIC ${ZFP_INC} ${SZ3_INC} ${ZSTD_INC} ${SZ_CPP_INC})
target_link_libraries(test_quantization_factor ${ZFP_LIB} ${ZSTD_LIB} ${SZ3_LIB} ${SZ_CPP_LIB})

