project(RT_cuda)
cmake_minimum_required(VERSION 2.8)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

#find_package(Threads)
#message(STATUS ${CMAKE_THREAD_LIBS_INIT})
#message(STATUS "WIN32: " ${CMAKE_USE_WIN32_THREADS_INIT})

# mandatory library
find_package(MXML REQUIRED)


# common library
set(COMMON_LIB rt_common)
aux_source_directory(./common-lib/ COMMONLIB_SRC_LIST)
add_library(${COMMON_LIB} ${COMMONLIB_SRC_LIST})
target_link_libraries(${COMMON_LIB} ${MXML_LIBRARIES})
set_target_properties(${COMMON_LIB} PROPERTIES INCLUDE_DIRECTORIES ${MXML_INCLUDE_DIR})


# Ray-tracing executable
set(RT_PROG ray-tracing)
add_executable(${RT_PROG} ./ray-tracing/main.cpp)
set_target_properties(${RT_PROG} PROPERTIES COMPILE_FLAGS "-g")
target_link_libraries(${RT_PROG} ${COMMON_LIB})
set_target_properties(${RT_PROG} PROPERTIES INCLUDE_DIRECTORIES ${MXML_INCLUDE_DIR})

# Quantum efficiency executable (surface related only part)
set(QE_PROG qe_comp)
add_executable(${QE_PROG} efficiency/efficiency.cpp)
target_link_libraries(${QE_PROG} ${COMMON_LIB})

target_link_libraries(${QE_PROG} ${MXML_LIBRARIES})
set_target_properties(${QE_PROG} PROPERTIES INCLUDE_DIRECTORIES ${MXML_INCLUDE_DIR})

set(LOSS_PROG loss_comp)
add_executable(${LOSS_PROG} efficiency/losses_comp.cpp)
target_link_libraries(${LOSS_PROG} ${COMMON_LIB})


# Auxiliary tools
set(RT_FILE_INFO rt_file_info)
add_executable(${RT_FILE_INFO} tools/rt_file_info.cpp)
target_link_libraries(${RT_FILE_INFO} ${COMMON_LIB})

set(PRINT_QE_FILE print_qe)
add_executable(${PRINT_QE_FILE} tools/print_QE.cpp)

set(GET_WORK_RANGE create_range_file)
add_executable(${GET_WORK_RANGE} tools/get_working_range.cpp)
target_link_libraries(${GET_WORK_RANGE} ${COMMON_LIB})

set(PRINT_LAMBDA print_lambda)
add_executable(${PRINT_LAMBDA} tools/print_lambda.cpp)
target_link_libraries(${PRINT_LAMBDA} ${COMMON_LIB})

find_package(OpenMP QUIET)
if (OPENMP_FOUND) # add OpenMP compile and link properties to the target
    target_compile_definitions(${RT_PROG} PUBLIC USE_OPENMP)

    set_target_properties(${RT_PROG} PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
    set_target_properties(${RT_PROG} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})

    message(STATUS "OpenMP is detected!")
endif()


find_package(CFITSIO QUIET)
if (CFITSIO_FOUND)
    message(STATUS "CFITSIO library was found!")
    message(STATUS "Add \"ccd_frame\" target")
    set(CCD_FRAME ccd_frame)
    add_executable(${CCD_FRAME} ./ccd/ccd_frame.cpp)
    target_link_libraries(${CCD_FRAME} ${COMMON_LIB})
    set_target_properties(${CCD_FRAME} PROPERTIES LINK_FLAGS ${CFITSIO_LIBRARY})
    set_target_properties(${CCD_FRAME} PROPERTIES INCLUDE_DIRECTORIES ${CFITSIO_INCLUDE_DIR})
else()
    message(STATUS "CFITSIO library was not found! Do not compile ccd_frame target!")
endif()


find_package(CUDA QUIET)
if (CUDA_FOUND)
    target_compile_definitions(${RT_PROG} PUBLIC USE_CUDA)
    message(STATUS "Found CUDA version " ${CUDA_VERSION})
    message(STATUS "Use of CUDA engine")

    if (CUDA_VERSION VERSION_LESS 6.0)
        if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            # require max gcc 4.7 (CUDA v. 5 restriction)
            if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
                message(STATUS "GCC version is ok (" ${CMAKE_CXX_COMPILER_VERSION} ")")
            else()
                message(FATAL_ERROR "For CUDA <=5.5 GCC version must be not greater than 4.7!")
            endif()
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            message(STATUS "Clang is detected")
        else()
            message(WARNING "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
        endif()
    endif()

    set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH}:${CUDA_INCLUDE_DIRS}) # CUDA include paths
    message(STATUS "CUDA include paths: " ${CMAKE_INCLUDE_PATH})

    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-arch=compute_50;-prec-sqrt=true;-prec-div=true )
#    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-arch=compute_20)
#    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-arch=compute_20;-prec-sqrt=true;-prec-div=true )

#     cuda_compile(RT_CUDA_ENGINE_COMMON_O rt_cuda_engine_common.cu)
#     cuda_compile(RT_CUDA_ENGINE_BEAM_O rt_cuda_engine_beam.cu)
#     cuda_compile(RT_CUDA_ENGINE_SPLINE_O rt_cuda_engine_spline.cu)
#     cuda_compile(RT_CUDA_ENGINE_SURFACE_O rt_cuda_engine_surface.cu)
# 
#     cuda_add_library(rt_cuda_engine ${RT_CUDA_ENGINE_COMMON_O} ${RT_CUDA_ENGINE_BEAM_O} ${RT_CUDA_ENGINE_SPLINE_O} ${RT_CUDA_ENGINE_SURFACE_O})
    
    FILE(GLOB CU_SRC ./cuda_engine/*.cu)
#    message(STATUS ${CU_SRC})
    
    cuda_compile(RT_CUDA_ENGINE_O ${CU_SRC})
    cuda_add_library(rt_cuda_engine ${RT_CUDA_ENGINE_O})
    
#    target_link_libraries(${RT_PROG} rt_cuda_engine ${CUDA_curand_LIBRARY})
    target_link_libraries(${COMMON_LIB} rt_cuda_engine ${CUDA_curand_LIBRARY})

#    if (CFITSIO_FOUND)
#        target_link_libraries(${CCD_FRAME} rt_cuda_engine)
#    endif()
else()    
    target_compile_definitions(${RT_PROG} PUBLIC USE_SIMD)
    message(STATUS "Use of Intel SIMD engine")

    # needs to detect SIMD instruction set
    set_target_properties(${RT_PROG} PROPERTIES COMPILE_FLAGS "-march=native")

    aux_source_directory(./simd_engine SIMD_SOURCES)

    add_library(rt_simd_engine STATIC ${SIMD_SOURCES})

#    target_link_libraries(${PROJECT_NAME} rt_simd_engine "-Lrt_simd_engine")
#    target_link_libraries(${RT_PROG} rt_simd_engine)
    target_link_libraries(${COMMON_LIB} rt_simd_engine)

#    if (CFITSIO_FOUND)
#        target_link_libraries(${CCD_FRAME} rt_simd_engine)
#    endif()
endif()
