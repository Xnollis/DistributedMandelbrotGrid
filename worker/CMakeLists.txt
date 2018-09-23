#在执行cmake之前，要先判断环境中的NVCC变量是否设置，以及是否将nvcc所在的路径添加到了PATH中
#export NVCC=/usr/local/cuda/bin/nvcc
#export PATH=$PATH:/usr/local/cuda/bin
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
# 项目名称
PROJECT(DMGWorker)
# cmake寻找cuda
enable_language(CUDA) #将来需要启用CUDA时，取消此行注释
# C++和CUDA的编译参数，可选。
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_30,code=sm_30;-std=c++11;)
# 头文件路径，按需
INCLUDE_DIRECTORIES(
    /usr/local/cuda/include
    ./mpir
    ./mpir_CUDA)
# 库文件路径，按需
LINK_DIRECTORIES(
    /usr/local/cuda/lib
    /usr/lib
    /usr/local/lib)
#嵌套的子目录
add_subdirectory(mpir)
add_subdirectory(mpir_CUDA)
#创建编译目标
ADD_EXECUTABLE(DMGWorker
    main.cpp
)
# 链接外部库，按需
TARGET_LINK_LIBRARIES(DMGWorker
    cudart
    mpir
    mpir_CUDA
)
