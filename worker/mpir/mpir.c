//#include "cuda_runtime.h"
#include<stdio.h>
#include<stdlib.h>
#include<memory.h>
#include "mpir_inter_decl.h"
/////////////////////////////////////////
const static int __gmp_default_fp_limb_precision=2;
/////////////////////////////////////////
void mpf_init (mpf_ptr r)
{
  mp_size_t prec = __gmp_default_fp_limb_precision;
  r->_mp_size = 0;
  r->_mp_exp = 0;
  r->_mp_prec = prec;
  memset(r->_mp_d, 0, sizeof(r->_mp_d));
}
void mpf_clear (mpf_ptr m)
{
	//Do nothing here
}
void mpf_set_d (mpf_ptr r, double d)
{
  int negative;
/* check if is a valid double value
  DOUBLE_NAN_INF_ACTION (d,
                         __gmp_invalid_operation (),
                         __gmp_invalid_operation ());
*/
  if (UNLIKELY (d == 0))
    {
      SIZ(r) = 0;
      EXP(r) = 0;
      return;
    }
  negative = d < 0;
  d = ABS (d);

  SIZ(r) = negative ? -LIMBS_PER_DOUBLE : LIMBS_PER_DOUBLE;
  EXP(r) = __gmp_extract_double (PTR(r), d);
}
/////////////////////////////////////////
#if 0
int IsCUDA_Supported(int bPrintInfoToConsole)
{
	cudaError_t err;
	int deviceCount;
	int device = 0,iGoodDeviceCount=0;

	err = cudaGetDeviceCount(&deviceCount);
	if (err != cudaSuccess)
		return 0;
	if (deviceCount <= 0)
		return 0;
	printf("找到%d个支持CUDA的显卡。以下是其详细信息：\n",deviceCount);
	struct cudaDeviceProp deviceProp;
	for (device = 0; device < deviceCount; ++device)
	{
		err = cudaGetDeviceProperties(&deviceProp, device);
		if (err == cudaSuccess)
		{
			if (bPrintInfoToConsole)
			{
				printf("显卡设备 %d 详情:\n%s 的计算能力(Compute Capability)为： %d.%d\n", device, deviceProp.name, deviceProp.major, deviceProp.minor);
				//枚举详细信息
				//printf("Identify: %s\n", deviceProp.name);
				//printf("Host Memory: %d\n", deviceProp.canMapHostMemory);
				printf("其他信息：\n显存大小： %zu 字节\n", deviceProp.totalGlobalMem);
				printf("核心处理器个数: %d\n", deviceProp.multiProcessorCount);
				printf("时钟频率 %d khz\n", deviceProp.clockRate);
#ifdef DEBUG
				printf("Compute Mode: %d\n", deviceProp.computeMode);
				printf("Device Overlap: %d\n", deviceProp.deviceOverlap);
				printf("Integrated: %d\n", deviceProp.integrated);
				printf("Kernel Exec Timeout Enabled: %d\n", deviceProp.kernelExecTimeoutEnabled);
				printf("Max Grid Size: %d * %d * %d\n", deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]);
				printf("Max Threads Dim: %d * %d * %d\n", deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]);
				printf("Max Threads per Block: %d\n", deviceProp.maxThreadsPerBlock);
				printf("Maximum Pitch: %d bytes\n", deviceProp.memPitch);
				printf("32bit Registers Availble per Block: %d\n", deviceProp.regsPerBlock);
				printf("Shared Memory Available per Block: %d bytes\n", deviceProp.sharedMemPerBlock);
				printf("Alignment Requirement for Textures: %d\n", deviceProp.textureAlignment);
				printf("Constant Memory Available: %d bytes\n", deviceProp.totalConstMem);
				printf("Warp Size: %d threads\n", deviceProp.warpSize);
#endif
				printf("=========================================\n");
			}
			if (deviceProp.major >= 2)//, deviceProp.minor);
				++iGoodDeviceCount;
		}
	}
	return (iGoodDeviceCount>0);
}
#endif