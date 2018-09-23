#include "cuda_runtime.h"
#include<stdio.h>
#include<memory.h>
#include"mpir_CUDA.h"
#include"mpir.h"
//extern "C"{
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
#include"mpn/toom8_sqr_n.c"
//}*/
__global__ void helloWorld() //__globa__ÊÇ¹Ø¼ü×Ö 
{
#pragma message(VAR_NAME_VALUE(mpf_init))
	mpf_t f1, f2;// , f3, f4, f5;
	mpf_init(f1);
	mpf_init(f2);
	mpf_set_d(f1, 123.456);
	//mpf_mul(f2, f1, f1);
	printf("haha, %f\n", mpf_get_d(f1));//*/
	printf("Hello CUDA!\n");
}