find_package(MPI REQUIRED)

find_library(ZFP_LIB zfp HINTS "/ccs/home/jieyang/xin/utils/zfp/lib64")
set (ZFP_INC "/ccs/home/jieyang/xin/utils/zfp/include")

add_executable (parallel_zfp parallel_zfp.cpp)
target_include_directories(parallel_zfp PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(parallel_zfp PUBLIC ${ZFP_INC})
target_link_libraries(parallel_zfp ${MPI_LIBRARIES} ${ZFP_LIB})
set_target_properties(parallel_zfp PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(parallel_zfp PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")

find_library(SZ_MASTER_LIB SZ HINTS "/ccs/home/jieyang/xin/codes/SZ/build/sz")
set (SZ_MASTER_INC "/ccs/home/jieyang/xin/codes/SZ/sz/include")
add_executable (parallel_sz parallel_sz.cpp)
target_include_directories(parallel_sz PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(parallel_sz PUBLIC ${SZ_MASTER_INC})
target_link_libraries(parallel_sz ${MPI_LIBRARIES} ${SZ_MASTER_LIB})
set_target_properties(parallel_sz PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(parallel_sz PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")

message(STATUS ${SZ3_INC})
add_executable (parallel_mgard parallel_mgard.cpp)
target_include_directories(parallel_mgard PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(parallel_mgard PUBLIC ${SZ3_INC} ${ZSTD_INC} ${SZ_CPP_INC})
target_link_libraries(parallel_mgard ${MPI_LIBRARIES} ${ZSTD_LIB} ${SZ3_LIB} ${SZ_CPP_LIB})
set_target_properties(parallel_mgard PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(parallel_mgard PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")

