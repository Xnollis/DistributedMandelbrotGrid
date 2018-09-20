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
double mpf_get_d(mpf_srcptr src)
{
	mp_size_t  size, abs_size;
	long       exp;

	size = SIZ(src);
	if (UNLIKELY(size == 0))
		return 0.0;

	abs_size = ABS(size);
	exp = (EXP(src) - abs_size) * GMP_NUMB_BITS;
	return mpn_get_d(PTR(src), abs_size, size, exp);//32-bit, 参数总大小应该是10h字节，但是现在由于有64位的参数，总大小变为18h了 《---由于注释这里有中文，会导致无法编译通过
}
void mpf_mul_ui(mpf_ptr r, mpf_srcptr u, mpir_ui v)
{
#if 0
	mp_srcptr up;
	mp_size_t usize;
	mp_size_t size;
	mp_size_t prec, excess;
	mp_limb_t cy_limb, vl, cbit, cin;
	mp_ptr rp;

	usize = u->_mp_size;
	if (UNLIKELY (v == 0) || UNLIKELY (usize == 0))
	{
		r->_mp_size = 0;
		r->_mp_exp = 0;
		return;
	}

#if BITS_PER_UI > GMP_NUMB_BITS  /* avoid warnings about shift amount */
	if (v > GMP_NUMB_MAX)
	{
		mpf_t     vf;
		mp_limb_t vp[2];
		vp[0] = v & GMP_NUMB_MASK;
		vp[1] = v >> GMP_NUMB_BITS;
		PTR(vf) = vp;
		SIZ(vf) = 2;
		ASSERT_CODE (PREC(vf) = 2);
		EXP(vf) = 2;
		mpf_mul (r, u, vf);
		return;
	}
#endif

	size = ABS (usize);
	prec = r->_mp_prec;
	up = u->_mp_d;
	vl = v;
	excess = size - prec;
	cin = 0;

	if (excess > 0)
	{
		/* up is bigger than desired rp, shorten it to prec limbs and
		determine a carry-in */

		mp_limb_t  vl_shifted = vl << GMP_NAIL_BITS;
		mp_limb_t  hi, lo, next_lo, sum;
		mp_size_t  i;

		/* high limb of top product */
		i = excess - 1;
		umul_ppmm(cin, lo, up[i], vl_shifted);

		/* and carry bit out of products below that, if any */
		for (;;)
		{
			i--;
			if (i < 0)
				break;

			umul_ppmm(hi, next_lo, up[i], vl_shifted);
			lo >>= GMP_NAIL_BITS;
			ADDC_LIMB(cbit, sum, hi, lo);
			cin += cbit;
			lo = next_lo;

			/* Continue only if the sum is GMP_NUMB_MAX.  GMP_NUMB_MAX is the
			only value a carry from below can propagate across.  If we've
			just seen the carry out (ie. cbit!=0) then sum!=GMP_NUMB_MAX,
			so this test stops us for that case too.  */
			if (LIKELY(sum != GMP_NUMB_MAX))
				break;
		}

		up += excess;
		size = prec;
	}

	rp = r->_mp_d;
#if HAVE_NATIVE_mpn_mul_1c
	cy_limb = mpn_mul_1c(rp, up, size, vl, cin);
#else
	cy_limb = mpn_mul_1(rp, up, size, vl);
	__GMPN_ADD_1(cbit, rp, rp, size, cin);
	cy_limb += cbit;
#endif
	rp[size] = cy_limb;
	cy_limb = cy_limb != 0;
	r->_mp_exp = u->_mp_exp + cy_limb;
	size += cy_limb;
	r->_mp_size = usize >= 0 ? size : -size;
#endif
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