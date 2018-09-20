#include "cuda_runtime.h"
#include<stdio.h>
#include"mpir_CUDA.h"
__global__ void helloWorld() //__globa__是关键字 
{ 
printf("Hello CUDA! haha\n"); 
}
void sayhello(void)
{
    helloWorld<<<3,2>>>(); //用3个block,每个block2个线程执行helloWorld() 
    cudaDeviceSynchronize(); //函数执行后必须有这个同步函数
}