#include "cuda_runtime.h"
#include<stdio.h>
#include<memory.h>
#include"mpir_CUDA.h"
#include"mpir.h"
__global__ void helloWorld();
void sayhello(void)
{
    helloWorld<<<3,2>>>(); //用3个block,每个block2个线程执行helloWorld() 
    cudaDeviceSynchronize(); //函数执行后必须有这个同步函数
}
#if 0
double MPIRWorking_DEVICE(double fStartValue)
{
	mpf_t f1, f2, f3, f4, f5;
	mpf_init(f1);
	mpf_init(f2);
	mpf_init(f3);
	mpf_init(f4);
	mpf_init(f5);
	mpf_set_d(f1, fStartValue);
	mpf_mul(f2, f1, f1);// f2=f1 * f1;15241.383936
	mpf_div_ui(f3, f2, 10);//f3=f2/10;1524.1383936
	mpf_set_d(f2, 1000L);
	mpf_mul(f3, f3, f2);//f3=f3*f2;1524138.3936
	mpf_div(f4, f3, f1);//f4=f3/f1;12345.6
	mpf_div_ui(f5, f4, 100L);//f4=f3/f1;123.456
	return mpf_get_d(f5);
}
#endif
#undef HAVE_UNISTD_H
//extern "C"{
//#pragma message(VAR_NAME_VALUE(THE_CLZ_TAB))
/*
#include"div.c"
#include"errno.c"
#include"extract-dbl.c"
#include"invalid.c"
#include"memory.c"
#include"mpir.c"
#include"mpn/add_n.c"
#include"mpn/divrem_1.c"
#include"mpn/divrem_2.c"
#include"mpn/divrem_euclidean_qr_1.c"
#include"mpn/divrem_euclidean_qr_2.c"
#include"mpn/divrem_euclidean_r_1.c"
#include"mpn/get_d.c"
#include"mpn/invert.c"
#include"mpn/mod_1.c"
#include"mpn/mul.c"
#include"mpn/mul_1.c"
#include"mpn/mul_n.c"
#include"mpn/neg_n.c"
#include"mpn/sqr_basecase.c"
#include"mpn/sub_n.c"
#include"mpn/tdiv_qr.c"
#include"mpn/toom3_mul.c"
#include"mpn/toom3_mul_n.c"
#include"mpn/toom4_mul.c"
#include"mpn/toom4_mul_n.c"
#include"mpn/toom8h_mul.c"
#include"mpn/toom8_sqr_n.c"//*/
//}