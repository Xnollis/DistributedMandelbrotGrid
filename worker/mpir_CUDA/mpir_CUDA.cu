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
	mpf_t f1, f2, f3;// , f3, f4, f5;
	mpf_init(f1);
	mpf_init(f2);
	mpf_init(f3);
	mpf_set_d(f1, 123.456);
	mpf_mul_ui(f2, f1, 10);
	//mpf_mul(f3, f1, f2); 
	printf("haha, %f,%f,%f\n", mpf_get_d(f1), mpf_get_d(f2), mpf_get_d(f3));//*/
	if (sizeof(void*) == 8)
		printf("MPIR sizes:sizeof(mp_limb_t)=%lld,sizeof(mp_limb_signed_t)=%lld,sizeof(mp_size_t)=%lld,sizeof(mp_exp_t)=%lld\n", sizeof(mp_limb_t), sizeof(mp_limb_signed_t), sizeof(mp_size_t), sizeof(mp_exp_t));
	else
		printf("MPIR sizes:sizeof(mp_limb_t)=%d,sizeof(mp_limb_signed_t)=%d,sizeof(mp_size_t)=%d,sizeof(mp_exp_t)=%d\n", sizeof(mp_limb_t), sizeof(mp_limb_signed_t), sizeof(mp_size_t), sizeof(mp_exp_t));
	printf("Hello CUDA!\n");
}