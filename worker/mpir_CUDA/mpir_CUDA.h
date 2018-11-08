#ifndef __MPIR_CUDA_H__
#define __MPIR_CUDA_H__
#include <stdio.h>
#include <stdlib.h>
#if defined (__cplusplus)
extern "C" {
#endif // __cplusplus

void sayhello(void);

#if defined (__cplusplus)
}
#endif // __cplusplus

#define MPF_TEST_AND_PRINT_UI(funcNameStr,src1,src2,result) funcNameStr(result,src1,src2);printf("%s(%f,%f)=%f\n",#funcNameStr,(double)mpf_get_d(src1),(double)(src2),(double)mpf_get_d(result))
#define MPF_TEST_AND_PRINT(funcNameStr,src1,src2,result) funcNameStr(result,src1,src2);printf("%s(%f,%f)=%f\n",#funcNameStr,(double)mpf_get_d(src1),(double)mpf_get_d(src2),(double)mpf_get_d(result))

//////////////////////////////////
///End of file
#endif // !__MPIR_CUDA_H__