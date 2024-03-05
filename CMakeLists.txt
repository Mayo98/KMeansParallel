cmake_minimum_required(VERSION 3.22)
project(KMeansParallel LANGUAGES CUDA CXX)
set(CMAKE_CUDA_COMPILER /usr/local/cuda-11.8/bin/nvcc)
set(CMAKE_CUDA_STANDARD 17)
#set(CMAKE_CUDA_COMPILER /usr/local/cuda-11/bin/nvcc)
add_executable(KMeansParallel main.cu
        Point.cpp
        Point.h
        Cluster.cpp
        Cluster.h
        SequentialKMeans.cpp
        SequentialKMeans.h
)
target_include_directories(KMeansParallel PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_compile_options(KMeansParallel PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
target_compile_features(KMeansParallel PUBLIC cxx_std_11)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math")
set_target_properties(KMeansParallel PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
