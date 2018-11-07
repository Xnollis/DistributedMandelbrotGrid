#ifndef __GMPN_ADD_H__
#define __GMPN_ADD_H__
/**************** mpn inlines ****************/
/* The comments with __GMPN_ADD_1 below apply here too.
 The test for FUNCTION returning 0 should predict well.  If it's assumed
 {yp,ysize} will usually have a random number of bits then the high limb
 won't be full and a carry out will occur a good deal less than 50% of the
 time.
 ysize==0 isn't a documented feature, but is used internally in a few
 places.
 Producing cout last stops it using up a register during the main part of
 the calculation, though gcc (as of 3.0) on an "if (mpn_add (...))"
 doesn't seem able to move the true and false legs of the conditional up
 to the two places cout is generated.  */
#define __GMPN_AORS(cout, wp, xp, xsize, yp, ysize, FUNCTION, TEST)     \
do {                                                                  \
mp_size_t  __gmp_i;                                                 \
mp_limb_t  __gmp_x;                                                 \
\
/* ASSERT ((ysize) >= 0); */                                        \
/* ASSERT ((xsize) >= (ysize)); */                                  \
/* ASSERT (MPN_SAME_OR_SEPARATE2_P (wp, xsize, xp, xsize)); */      \
/* ASSERT (MPN_SAME_OR_SEPARATE2_P (wp, xsize, yp, ysize)); */      \
\
__gmp_i = (ysize);                                                  \
if (__gmp_i != 0)                                                   \
{                                                                 \
if (FUNCTION (wp, xp, yp, __gmp_i))                             \
{                                                             \
do                                                          \
{                                                         \
if (__gmp_i >= (xsize))                                 \
{                                                     \
(cout) = 1;                                         \
goto __gmp_done;                                    \
}                                                     \
__gmp_x = (xp)[__gmp_i];                                \
}                                                         \
while (TEST);                                               \
}                                                             \
}                                                                 \
if ((wp) != (xp))                                                   \
__GMPN_COPY_REST (wp, xp, xsize, __gmp_i);                        \
(cout) = 0;                                                         \
__gmp_done:                                                           \
;                                                                   \
} while (0)
#define __GMPN_ADD(cout, wp, xp, xsize, yp, ysize)              \
__GMPN_AORS (cout, wp, xp, xsize, yp, ysize, mpn_add_n,       \
(((wp)[__gmp_i++] = (__gmp_x + 1) & GMP_NUMB_MASK) == 0))
#define __GMPN_SUB(cout, wp, xp, xsize, yp, ysize)              \
__GMPN_AORS (cout, wp, xp, xsize, yp, ysize, mpn_sub_n,       \
(((wp)[__gmp_i++] = (__gmp_x - 1) & GMP_NUMB_MASK), __gmp_x == 0))
/* The use of __gmp_i indexing is designed to ensure a compile time src==dst
 remains nice and clear to the compiler, so that __GMPN_COPY_REST can
 disappear, and the load/add/store gets a chance to become a
 read-modify-write on CISC CPUs.
 Alternatives:
 Using a pair of pointers instead of indexing would be possible, but gcc
 isn't able to recognise compile-time src==dst in that case, even when the
 pointers are incremented more or less together.  Other compilers would
 very likely have similar difficulty.
 gcc could use "if (__builtin_constant_p(src==dst) && src==dst)" or
 similar to detect a compile-time src==dst.  This works nicely on gcc
 2.95.x, it's not good on gcc 3.0 where __builtin_constant_p(p==p) seems
 to be always false, for a pointer p.  But the current code form seems
 good enough for src==dst anyway.
 gcc on x86 as usual doesn't give particularly good flags handling for the
 carry/borrow detection.  It's tempting to want some multi instruction asm
 blocks to help it, and this was tried, but in truth there's only a few
 instructions to save and any gain is all too easily lost by register
 juggling setting up for the asm.  */
#if GMP_NAIL_BITS == 0
#define __GMPN_AORS_1(cout, dst, src, n, v, OP, CB)        \
do {                                \
mp_size_t  __gmp_i;                        \
mp_limb_t  __gmp_x, __gmp_r;                                \
\
/* ASSERT ((n) >= 1); */                    \
/* ASSERT (MPN_SAME_OR_SEPARATE_P (dst, src, n)); */    \
\
__gmp_x = (src)[0];                        \
__gmp_r = __gmp_x OP (v);                                   \
(dst)[0] = __gmp_r;                        \
if (CB (__gmp_r, __gmp_x, (v)))                             \
{                                \
(cout) = 1;                        \
for (__gmp_i = 1; __gmp_i < (n);)                       \
{                            \
__gmp_x = (src)[__gmp_i];                           \
__gmp_r = __gmp_x OP 1;                             \
(dst)[__gmp_i] = __gmp_r;                           \
++__gmp_i;                        \
if (!CB (__gmp_r, __gmp_x, 1))                      \
{                            \
if ((src) != (dst))                \
__GMPN_COPY_REST (dst, src, n, __gmp_i);      \
(cout) = 0;                    \
break;                        \
}                            \
}                            \
}                                \
else                            \
{                                \
if ((src) != (dst))                    \
__GMPN_COPY_REST (dst, src, n, 1);            \
(cout) = 0;                        \
}                                \
} while (0)
#endif
#if GMP_NAIL_BITS >= 1
#define __GMPN_AORS_1(cout, dst, src, n, v, OP, CB)        \
do {                                \
mp_size_t  __gmp_i;                        \
mp_limb_t  __gmp_x, __gmp_r;                \
\
/* ASSERT ((n) >= 1); */                    \
/* ASSERT (MPN_SAME_OR_SEPARATE_P (dst, src, n)); */    \
\
__gmp_x = (src)[0];                        \
__gmp_r = __gmp_x OP (v);                    \
(dst)[0] = __gmp_r & GMP_NUMB_MASK;                \
if (__gmp_r >> GMP_NUMB_BITS != 0)                \
{                                \
(cout) = 1;                        \
for (__gmp_i = 1; __gmp_i < (n);)            \
{                            \
__gmp_x = (src)[__gmp_i];                \
__gmp_r = __gmp_x OP 1;                \
(dst)[__gmp_i] = __gmp_r & GMP_NUMB_MASK;        \
++__gmp_i;                        \
if (__gmp_r >> GMP_NUMB_BITS == 0)            \
{                            \
if ((src) != (dst))                \
__GMPN_COPY_REST (dst, src, n, __gmp_i);    \
(cout) = 0;                    \
break;                        \
}                            \
}                            \
}                                \
else                            \
{                                \
if ((src) != (dst))                    \
__GMPN_COPY_REST (dst, src, n, 1);            \
(cout) = 0;                        \
}                                \
} while (0)
#endif
#define __GMPN_ADDCB(r,x,y) ((r) < (y))
#define __GMPN_SUBCB(r,x,y) ((x) < (y))
#define __GMPN_ADD_1(cout, dst, src, n, v)         \
__GMPN_AORS_1(cout, dst, src, n, v, +, __GMPN_ADDCB)
#define __GMPN_SUB_1(cout, dst, src, n, v)         \
__GMPN_AORS_1(cout, dst, src, n, v, -, __GMPN_SUBCB)
/* Compare {xp,size} and {yp,size}, setting "result" to positive, zero or
 negative.  size==0 is allowed.  On random data usually only one limb will
 need to be examined to get a result, so it's worth having it inline.  */
#define __GMPN_CMP(result, xp, yp, size)                                \
do {                                                                  \
mp_size_t  __gmp_i;                                                 \
mp_limb_t  __gmp_x, __gmp_y;                                        \
\
/* ASSERT ((size) >= 0); */                                         \
\
(result) = 0;                                                       \
__gmp_i = (size);                                                   \
while (--__gmp_i >= 0)                                              \
{                                                                 \
__gmp_x = (xp)[__gmp_i];                                        \
__gmp_y = (yp)[__gmp_i];                                        \
if (__gmp_x != __gmp_y)                                         \
{                                                             \
/* Cannot use __gmp_x - __gmp_y, may overflow an "int" */   \
(result) = (__gmp_x > __gmp_y ? 1 : -1);                    \
break;                                                      \
}                                                             \
}                                                                 \
} while (0)
#if defined (__GMPN_COPY) && ! defined (__GMPN_COPY_REST)
#define __GMPN_COPY_REST(dst, src, size, start)                 \
do {                                                          \
/* ASSERT ((start) >= 0); */                                \
/* ASSERT ((start) <= (size)); */                           \
__GMPN_COPY ((dst)+(start), (src)+(start), (size)-(start)); \
} while (0)
#endif
/* Copy {src,size} to {dst,size}, starting at "start".  This is designed to
 keep the indexing dst[j] and src[j] nice and simple for __GMPN_ADD_1,
 __GMPN_ADD, etc.  */
#if ! defined (__GMPN_COPY_REST)
#define __GMPN_COPY_REST(dst, src, size, start)                 \
do {                                                          \
mp_size_t __gmp_j;                                          \
/* ASSERT ((size) >= 0); */                                 \
/* ASSERT ((start) >= 0); */                                \
/* ASSERT ((start) <= (size)); */                           \
/* ASSERT (MPN_SAME_OR_SEPARATE_P (dst, src, size)); */     \
for (__gmp_j = (start); __gmp_j < (size); __gmp_j++)        \
(dst)[__gmp_j] = (src)[__gmp_j];                          \
} while (0)
#endif
/* Enhancement: Use some of the smarter code from gmp-impl.h.  Maybe use
 mpn_copyi if there's a native version, and if we don't mind demanding
 binary compatibility for it (on targets which use it).  */
#if ! defined (__GMPN_COPY)
#define __GMPN_COPY(dst, src, size)   __GMPN_COPY_REST (dst, src, size, 0)
#endif



/* Use a library function for invert_limb, if available. */
#define mpn_invert_limb  __MPN(invert_limb)

#if ! defined (invert_limb) && HAVE_NATIVE_mpn_invert_limb
#define invert_limb(invxl,xl)           \
do {                                  \
(invxl) = mpn_invert_limb (xl);     \
} while (0)
#endif

#ifndef invert_limb
#define invert_limb(invxl,xl)                   \
do {                                          \
mp_limb_t dummy=0;                            \
ASSERT ((xl) != 0);                         \
udiv_qrnnd (invxl, dummy, ~(xl), ~CNST_LIMB(0), xl);  \
} while (0)
#endif



/* ADDC_LIMB sets w=x+y and cout to 0 or 1 for a carry from that addition. */
#if GMP_NAIL_BITS == 0
#define ADDC_LIMB(cout, w, x, y)        \
do {                                  \
mp_limb_t  __x = (x);               \
mp_limb_t  __y = (y);               \
mp_limb_t  __w = __x + __y;         \
(w) = __w;                          \
(cout) = __w < __x;                 \
} while (0)
#else
#define ADDC_LIMB(cout, w, x, y)        \
do {                                  \
mp_limb_t  __w;                     \
ASSERT_LIMB (x);                    \
ASSERT_LIMB (y);                    \
__w = (x) + (y);                    \
(w) = __w & GMP_NUMB_MASK;          \
(cout) = __w >> GMP_NUMB_BITS;      \
} while (0)
#endif

/* SUBC_LIMB sets w=x-y and cout to 0 or 1 for a borrow from that
 subtract.  */
#if GMP_NAIL_BITS == 0
#define SUBC_LIMB(cout, w, x, y)        \
do {                                  \
mp_limb_t  __x = (x);               \
mp_limb_t  __y = (y);               \
mp_limb_t  __w = __x - __y;         \
(w) = __w;                          \
(cout) = __w > __x;                 \
} while (0)
#else
#define SUBC_LIMB(cout, w, x, y)        \
do {                                  \
mp_limb_t  __w = (x) - (y);         \
(w) = __w & GMP_NUMB_MASK;          \
(cout) = __w >> (GMP_LIMB_BITS-1);  \
} while (0)
#endif
//////////////////////////////////////////////////
#define mpn_divmod_1(qp,np,nsize,dlimb) \
mpn_divrem_1 (qp, __GMP_CAST (mp_size_t, 0), np, nsize, dlimb)


#define modlimb_invert(inv,n)                        \
do {                                    \
mp_limb_t  __n = (n);                        \
mp_limb_t  __inv;                            \
ASSERT ((__n & 1) == 1);                        \
\
__inv = getValFrom_invert_table((__n / 2) & 0x7F); /*  8 */        \
if (GMP_NUMB_BITS > 8)   __inv = 2 * __inv - __inv * __inv * __n;    \
if (GMP_NUMB_BITS > 16)  __inv = 2 * __inv - __inv * __inv * __n;    \
if (GMP_NUMB_BITS > 32)  __inv = 2 * __inv - __inv * __inv * __n;    \
\
if (GMP_NUMB_BITS > 64)                        \
{                                    \
int  __invbits = 64;                        \
do {                                \
__inv = 2 * __inv - __inv * __inv * __n;            \
__invbits *= 2;                        \
} while (__invbits < GMP_NUMB_BITS);                \
}                                    \
\
ASSERT ((__inv * __n & GMP_NUMB_MASK) == 1);            \
(inv) = __inv & GMP_NUMB_MASK;                    \
} while (0)

/* On the x86s repe/scasl doesn't seem useful, since it takes many cycles to
start up and would need to strip a lot of zeros before it'd be faster
than a simple cmpl loop.  Here are some times in cycles for
std/repe/scasl/cld and cld/repe/scasl (the latter would be for stripping
low zeros).

std   cld
P5    18    16
P6    46    38
K6    36    13
K7    21    20
*/
#ifndef MPN_NORMALIZE
#define MPN_NORMALIZE(DST, NLIMBS) \
  do {									\
      while ((NLIMBS) > 0)                                                \
	        {									\
    if ((DST)[(NLIMBS) - 1] != 0)					\
      break;							\
    (NLIMBS)--;							\
	        }									\
    } while (0)
#endif
#ifndef MPN_NORMALIZE_NOT_ZERO
#define MPN_NORMALIZE_NOT_ZERO(DST, NLIMBS)     \
  do {                                          \
    ASSERT ((NLIMBS) >= 1);                     \
	    while (1)                                   \
		      {                                         \
    if ((DST)[(NLIMBS) - 1] != 0)           \
      break;                                \
    (NLIMBS)--;                             \
		      }                                         \
    } while (0)
#endif

/* Strip least significant zero limbs from {ptr,size} by incrementing ptr
and decrementing size.  low should be ptr[0], and will be the new ptr[0]
on returning.  The number in {ptr,size} must be non-zero, ie. size!=0 and
somewhere a non-zero limb.  */
#define MPN_STRIP_LOW_ZEROS_NOT_ZERO(ptr, size, low)    \
  do {                                                  \
    ASSERT ((size) >= 1);                               \
    ASSERT ((low) == (ptr)[0]);                         \
                                                        \
														    while ((low) == 0)                                  \
															      {                                                 \
        (size)--;                                       \
        ASSERT ((size) >= 1);                           \
        (ptr)++;                                        \
        (low) = *(ptr);                                 \
															      }                                                 \
    } while (0)

/* Multiplicative inverse of 3, modulo 2^GMP_NUMB_BITS.
 Eg. 0xAAAAAAAB for 32 bits, 0xAAAAAAAAAAAAAAAB for 64 bits.
 GMP_NUMB_MAX/3*2+1 is right when GMP_NUMB_BITS is even, but when it's odd
 we need to start from GMP_NUMB_MAX>>1. */
#define MODLIMB_INVERSE_3 (((GMP_NUMB_MAX >> (GMP_NUMB_BITS % 2)) / 3) * 2 + 1)

/* ceil(GMP_NUMB_MAX/3) and ceil(2*GMP_NUMB_MAX/3).
 These expressions work because GMP_NUMB_MAX%3 != 0 for all GMP_NUMB_BITS. */
#define GMP_NUMB_CEIL_MAX_DIV3   (GMP_NUMB_MAX / 3 + 1)
#define GMP_NUMB_CEIL_2MAX_DIV3  ((GMP_NUMB_MAX>>1) / 3 + 1 + GMP_NUMB_HIGHBIT)


/* For a threshold between algorithms A and B, size>=thresh is where B
 should be used.  Special value MP_SIZE_T_MAX means only ever use A, or
 value 0 means only ever use B.  The tests for these special values will
 be compile-time constants, so the compiler should be able to eliminate
 the code for the unwanted algorithm.  */

#define ABOVE_THRESHOLD(size,thresh)    \
((thresh) == 0                        \
|| ((thresh) != MP_SIZE_T_MAX        \
&& (size) >= (thresh)))
#define BELOW_THRESHOLD(size,thresh)  (! ABOVE_THRESHOLD (size, thresh))

/* If MUL_KARATSUBA_THRESHOLD is not already defined, define it to a
 value which is good on most machines.  */
#ifndef MUL_KARATSUBA_THRESHOLD
#define MUL_KARATSUBA_THRESHOLD 32
#endif

#ifndef SQR_KARATSUBA_THRESHOLD
#define SQR_KARATSUBA_THRESHOLD 32
#endif

/* If MUL_TOOM3_THRESHOLD is not already defined, define it to a
 value which is good on most machines.  */
#ifndef MUL_TOOM3_THRESHOLD
#define MUL_TOOM3_THRESHOLD 128
#endif

#ifndef MUL_TOOM4_THRESHOLD
#define MUL_TOOM4_THRESHOLD 300
#endif

#ifndef MULMID_TOOM42_THRESHOLD
#define MULMID_TOOM42_THRESHOLD 36
#endif

#ifndef MUL_TOOM8H_THRESHOLD
#define MUL_TOOM8H_THRESHOLD 401
#endif

#ifndef SQR_TOOM3_THRESHOLD
#define SQR_TOOM3_THRESHOLD 128
#endif

#ifndef SQR_TOOM4_THRESHOLD
#define SQR_TOOM4_THRESHOLD 300
#endif

#ifndef SQR_TOOM8_THRESHOLD
#define SQR_TOOM8_THRESHOLD 401
#endif

#ifndef MULLOW_BASECASE_THRESHOLD
#define MULLOW_BASECASE_THRESHOLD    8
#endif

#ifndef MULLOW_DC_THRESHOLD
#define MULLOW_DC_THRESHOLD    32
#endif

#ifndef MULLOW_MUL_THRESHOLD
#define MULLOW_MUL_THRESHOLD    8192
#endif

#ifndef MULHIGH_BASECASE_THRESHOLD
#define MULHIGH_BASECASE_THRESHOLD    16
#endif

#ifndef MULHIGH_DC_THRESHOLD
#define MULHIGH_DC_THRESHOLD    32
#endif

#ifndef MULHIGH_MUL_THRESHOLD
#define MULHIGH_MUL_THRESHOLD    8192
#endif

#ifndef MULMOD_2EXPM1_THRESHOLD
#define MULMOD_2EXPM1_THRESHOLD    16
#endif

#ifndef FAC_UI_THRESHOLD
#define FAC_UI_THRESHOLD    8192
#endif

#ifndef ROOTREM_THRESHOLD
#define ROOTREM_THRESHOLD    8
#endif

#ifndef DIVREM_HENSEL_QR_1_THRESHOLD
#define DIVREM_HENSEL_QR_1_THRESHOLD 8
#endif

#ifndef RSH_DIVREM_HENSEL_QR_1_THRESHOLD
#define RSH_DIVREM_HENSEL_QR_1_THRESHOLD 8
#endif

#ifndef DIVREM_EUCLID_HENSEL_THRESHOLD
#define DIVREM_EUCLID_HENSEL_THRESHOLD 32
#endif

#ifndef MOD_1_1_THRESHOLD
#define MOD_1_1_THRESHOLD    16
#endif

#ifndef MOD_1_2_THRESHOLD
#define MOD_1_2_THRESHOLD    32
#endif

#ifndef MOD_1_3_THRESHOLD
#define MOD_1_3_THRESHOLD    64
#endif

/* MUL_KARATSUBA_THRESHOLD_LIMIT is the maximum for MUL_KARATSUBA_THRESHOLD.
 In a normal build MUL_KARATSUBA_THRESHOLD is a constant and we use that.
 In a fat binary or tune program build MUL_KARATSUBA_THRESHOLD is a
 variable and a separate hard limit will have been defined.  Similarly for
 TOOM3.  */
#ifndef MUL_KARATSUBA_THRESHOLD_LIMIT
#define MUL_KARATSUBA_THRESHOLD_LIMIT  MUL_KARATSUBA_THRESHOLD
#endif
#ifndef MUL_TOOM3_THRESHOLD_LIMIT
#define MUL_TOOM3_THRESHOLD_LIMIT  MUL_TOOM3_THRESHOLD
#endif
#ifndef MUL_TOOM4_THRESHOLD_LIMIT
#define MUL_TOOM4_THRESHOLD_LIMIT  MUL_TOOM4_THRESHOLD
#endif
#ifndef MUL_TOOM8H_THRESHOLD_LIMIT
#define MUL_TOOM8H_THRESHOLD_LIMIT  MUL_TOOM8H_THRESHOLD
#endif
#ifndef MULLOW_BASECASE_THRESHOLD_LIMIT
#define MULLOW_BASECASE_THRESHOLD_LIMIT  MULLOW_BASECASE_THRESHOLD
#endif

/* SQR_BASECASE_THRESHOLD is where mpn_sqr_basecase should take over from
mpn_mul_basecase in mpn_sqr_n.  Default is to use mpn_sqr_basecase
always.  (Note that we certainly always want it if there's a native
assembler mpn_sqr_basecase.)

If it turns out that mpn_kara_sqr_n becomes faster than mpn_mul_basecase
before mpn_sqr_basecase does, then SQR_BASECASE_THRESHOLD is the
karatsuba threshold and SQR_KARATSUBA_THRESHOLD is 0.  This oddity arises
more or less because SQR_KARATSUBA_THRESHOLD represents the size up to
which mpn_sqr_basecase should be used, and that may be never.  */

#ifndef SQR_BASECASE_THRESHOLD
#define SQR_BASECASE_THRESHOLD 0
#endif

#ifndef SQR_KARATSUBA_THRESHOLD
#define SQR_KARATSUBA_THRESHOLD (2*MUL_KARATSUBA_THRESHOLD)
#endif

#ifndef SQR_TOOM3_THRESHOLD
#define SQR_TOOM3_THRESHOLD 128
#endif

#ifndef SQR_TOOM4_THRESHOLD
#define SQR_TOOM4_THRESHOLD 300
#endif

#ifndef SQR_TOOM8_THRESHOLD
#define SQR_TOOM8_THRESHOLD 400
#endif

/* See comments above about MUL_TOOM3_THRESHOLD_LIMIT.  */
#ifndef SQR_TOOM3_THRESHOLD_LIMIT
#define SQR_TOOM3_THRESHOLD_LIMIT  SQR_TOOM3_THRESHOLD
#endif

#ifndef SQR_TOOM4_THRESHOLD_LIMIT
#define SQR_TOOM4_THRESHOLD_LIMIT  SQR_TOOM4_THRESHOLD
#endif

#ifndef SQR_TOOM8_THRESHOLD_LIMIT
#define SQR_TOOM8_THRESHOLD_LIMIT  SQR_TOOM8_THRESHOLD
#endif

/* points at which fft is used for mul/sqr and mulmod_Bexp resp. */
#ifndef MUL_FFT_FULL_THRESHOLD
#define MUL_FFT_FULL_THRESHOLD   (MUL_TOOM8H_THRESHOLD * 10)
#endif
#ifndef SQR_FFT_FULL_THRESHOLD
#define SQR_FFT_FULL_THRESHOLD   (SQR_TOOM8_THRESHOLD * 10)
#endif

#ifndef MUL_FFT_THRESHOLD
#define MUL_FFT_THRESHOLD   (MUL_FFT_FULL_THRESHOLD / 2)
#endif
#ifndef SQR_FFT_THRESHOLD
#define SQR_FFT_THRESHOLD   (SQR_FFT_FULL_THRESHOLD / 2)
#endif

#ifndef FFT_MULMOD_2EXPP1_CUTOFF
#define FFT_MULMOD_2EXPP1_CUTOFF 128
#endif

#ifndef DC_DIV_QR_THRESHOLD
#define DC_DIV_QR_THRESHOLD    (3 * MUL_KARATSUBA_THRESHOLD)
#endif

#ifndef DC_DIVAPPR_Q_N_THRESHOLD
#define DC_DIVAPPR_Q_N_THRESHOLD    (3 * MUL_KARATSUBA_THRESHOLD)
#endif

#ifndef DC_BDIV_QR_THRESHOLD
#define DC_BDIV_QR_THRESHOLD    (3 * MUL_KARATSUBA_THRESHOLD)
#endif

#ifndef DC_BDIV_Q_THRESHOLD
#define DC_BDIV_Q_THRESHOLD    (3 * MUL_KARATSUBA_THRESHOLD)
#endif

#ifndef INV_DIV_QR_THRESHOLD
#define INV_DIV_QR_THRESHOLD    (MUL_FFT_THRESHOLD/3)
#endif

#ifndef INV_DIVAPPR_Q_N_THRESHOLD
#define INV_DIVAPPR_Q_N_THRESHOLD    (MUL_FFT_THRESHOLD/3)
#endif

#ifndef DC_DIV_Q_THRESHOLD
#define DC_DIV_Q_THRESHOLD    (3 * MUL_KARATSUBA_THRESHOLD)
#endif

#ifndef INV_DIV_Q_THRESHOLD
#define INV_DIV_Q_THRESHOLD    (MUL_FFT_THRESHOLD/3)
#endif

#ifndef BINV_NEWTON_THRESHOLD
#define BINV_NEWTON_THRESHOLD           300
#endif

#ifndef DC_DIVAPPR_Q_THRESHOLD
#define DC_DIVAPPR_Q_THRESHOLD    (3 * MUL_TOOM3_THRESHOLD)
#endif

#ifndef INV_DIVAPPR_Q_THRESHOLD
#define INV_DIVAPPR_Q_THRESHOLD    (MUL_FFT_THRESHOLD/2)
#endif

#ifndef GET_STR_DC_THRESHOLD
#define GET_STR_DC_THRESHOLD             18
#endif

#ifndef GET_STR_PRECOMPUTE_THRESHOLD
#define GET_STR_PRECOMPUTE_THRESHOLD     35
#endif

#ifndef SET_STR_DC_THRESHOLD
#define SET_STR_DC_THRESHOLD            750
#endif

#ifndef SET_STR_PRECOMPUTE_THRESHOLD
#define SET_STR_PRECOMPUTE_THRESHOLD   2000
#endif

#ifndef FAC_ODD_THRESHOLD
#define FAC_ODD_THRESHOLD    35
#endif

#ifndef FAC_DSC_THRESHOLD
#define FAC_DSC_THRESHOLD   400
#endif

#define mpn_toom42_mulmid_itch(n) (3*(n) + 64)

/* kara uses n+1 limbs of temporary space and then recurses with the balance,
so need (n+1) + (ceil(n/2)+1) + (ceil(n/4)+1) + ...  This can be solved to
2n + o(n).  Since n is very limited, o(n) in practice could be around 15.
For now, assume n is arbitrarily large.  */
#define MPN_KARA_MUL_N_TSIZE(n)   (2*(n) + 2*GMP_LIMB_BITS)
#define MPN_KARA_SQR_N_TSIZE(n)   (2*(n) + 2*GMP_LIMB_BITS)

/* toom3 uses 2n + 2n/3 + o(n) limbs of temporary space if mpn_sublsh1_n is
unavailable, but just 2n + o(n) if mpn_sublsh1_n is available.  It is hard
to pin down the value of o(n), since it is a complex function of
MUL_TOOM3_THRESHOLD and n.  Normally toom3 is used between kara and fft; in
that case o(n) will be really limited.  If toom3 is used for arbitrarily
large operands, o(n) will be larger.  These definitions handle operands of
up to 8956264246117233 limbs.  A single multiplication using toom3 on the
fastest hardware currently (2003) would need 100 million years, which
suggests that these limits are acceptable.  */
#if WANT_FFT
#if HAVE_NATIVE_mpn_sublsh1_n
#define MPN_TOOM3_MUL_N_TSIZE(n)  (2*(n) + 63)
#define MPN_TOOM3_MUL_TSIZE(n)    (3*(n) + 63)
#define MPN_TOOM3_SQR_N_TSIZE(n)  (2*(n) + 63)
#else
#define MPN_TOOM3_MUL_N_TSIZE(n)  (2*(n) + 2*(n/3) + 63)
#define MPN_TOOM3_MUL_TSIZE(n)    (3*(n) + 3*(n/3) + 63)
#define MPN_TOOM3_SQR_N_TSIZE(n)  (2*(n) + 2*(n/3) + 63)
#endif
#else /* WANT_FFT */
#if HAVE_NATIVE_mpn_sublsh1_n
#define MPN_TOOM3_MUL_N_TSIZE(n)  (2*(n) + 255)
#define MPN_TOOM3_MUL_TSIZE(n)    (3*(n) + 255)
#define MPN_TOOM3_SQR_N_TSIZE(n)  (2*(n) + 255)
#else
#define MPN_TOOM3_MUL_N_TSIZE(n)  (2*(n) + 2*(n/3) + 255)
#define MPN_TOOM3_MUL_TSIZE(n)    (3*(n) + 3*(n/3) + 255)
#define MPN_TOOM3_SQR_N_TSIZE(n)  (2*(n) + 2*(n/3) + 255)
#endif
#define MPN_TOOM3_MAX_N 285405
#endif /* WANT_FFT */

/* need 2 so that n2>=1 */
#if defined(HAVE_NATIVE_mpn_karaadd) || defined(HAVE_NATIVE_mpn_karasub)
#define MPN_KARA_MUL_N_MINSIZE    8
#define MPN_KARA_SQR_N_MINSIZE    8
#else
#define MPN_KARA_MUL_N_MINSIZE    2
#define MPN_KARA_SQR_N_MINSIZE    2
#endif

/* Need l>=1, ls>=1, and 2*ls > l (the latter for the tD MPN_INCR_U) */
#define MPN_TOOM3_MUL_N_MINSIZE   17 
#define MPN_TOOM4_MUL_N_MINSIZE   32
#define MPN_TOOM8H_MUL_MINSIZE    86
#define MPN_TOOM3_SQR_N_MINSIZE   17
#define MPN_TOOM4_SQR_N_MINSIZE   32
#define MPN_TOOM8_SQR_N_MINSIZE   58
#define MPN_FFT_MUL_N_MINSIZE     64







#ifndef udiv_qrnnd_preinv
#define udiv_qrnnd_preinv udiv_qrnnd_preinv2
#endif

/* Like udiv_qrnnd_preinv, but branch-free. */
#define udiv_qrnnd_preinv2(q, r, nh, nl, d, di)                \
do {                                    \
mp_limb_t _n2, _n10, _nmask, _nadj, _q1;                \
mp_limb_t _xh, _xl;                            \
_n2 = (nh);                                \
_n10 = (nl);                            \
_nmask = LIMB_HIGHBIT_TO_MASK (_n10);                \
_nadj = _n10 + (_nmask & (d));                    \
umul_ppmm (_xh, _xl, di, _n2 - _nmask);                \
add_ssaaaa (_xh, _xl, _xh, _xl, _n2, _nadj);            \
_q1 = ~_xh;                                \
umul_ppmm (_xh, _xl, _q1, d);                    \
add_ssaaaa (_xh, _xl, _xh, _xl, nh, nl);                \
_xh -= (d);                    /* xh = 0 or -1 */    \
(r) = _xl + ((d) & _xh);                        \
(q) = _xh - _q1;                            \
} while (0)
/***************************************************************/
#define LIMB_HIGHBIT_TO_MASK(n)                                 \
(((mp_limb_signed_t) -1 >> 1) < 0                             \
? (mp_limb_signed_t) (n) >> (GMP_LIMB_BITS - 1)              \
: (n) & GMP_LIMB_HIGHBIT ? MP_LIMB_T_MAX : CNST_LIMB(0))
/***************************************************************/

/* MPN_INCR_U does {ptr,size} += n, MPN_DECR_U does {ptr,size} -= n, both
expecting no carry (or borrow) from that.

The size parameter is only for the benefit of assertion checking.  In a
normal build it's unused and the carry/borrow is just propagated as far
as it needs to go.

On random data, usually only one or two limbs of {ptr,size} get updated,
so there's no need for any sophisticated looping, just something compact
and sensible.

FIXME: Switch all code from mpn_{incr,decr}_u to MPN_{INCR,DECR}_U,
declaring their operand sizes, then remove the former.  This is purely
for the benefit of assertion checking.  */

/* Dummy for non-gcc, code involving it will go dead. */
#if ! defined (__GNUC__) || __GNUC__ < 2
#define __builtin_constant_p(x)   0
#endif
#if GMP_NAIL_BITS == 0
#ifndef mpn_incr_u
#define mpn_incr_u(p,incr)                              \
  do {                                                  \
    mp_limb_t __x;                                      \
    mp_ptr __p = (p);                                   \
    if (__builtin_constant_p (incr) && (incr) == 1)     \
	      {                                                 \
		          while (++(*(__p++)) == 0)                       \
          ;                                             \
	      }                                                 \
	    else                                                \
	        {                                                 \
        __x = *__p + (incr);                            \
        *__p = __x;                                     \
        if (__x < (incr))                               \
		          while (++(*(++__p)) == 0)                     \
            ;                                           \
	        }                                                 \
    } while (0)
#endif
#ifndef mpn_decr_u
#define mpn_decr_u(p,incr)                              \
  do {                                                  \
    mp_limb_t __x;                                      \
    mp_ptr __p = (p);                                   \
    if (__builtin_constant_p (incr) && (incr) == 1)     \
	      {                                                 \
		          while ((*(__p++))-- == 0)                       \
          ;                                             \
	      }                                                 \
	    else                                                \
	        {                                                 \
        __x = *__p;                                     \
        *__p = __x - (incr);                            \
        if (__x < (incr))                               \
		          while ((*(++__p))-- == 0)                     \
            ;                                           \
	        }                                                 \
    } while (0)
#endif
#endif

#if GMP_NAIL_BITS >= 1
#ifndef mpn_incr_u
#define mpn_incr_u(p,incr)                              \
  do {							\
    mp_limb_t __x;					\
    mp_ptr __p = (p);					\
    if (__builtin_constant_p (incr) && (incr) == 1)	\
	      {							\
    do						\
	      {						\
        __x = (*__p + 1) & GMP_NUMB_MASK;		\
        *__p++ = __x;				\
	      }						\
		      while (__x == 0);				\
	      }							\
	    else						\
	        {							\
    __x = (*__p + (incr));				\
    *__p++ = __x & GMP_NUMB_MASK;			\
    if (__x >> GMP_NUMB_BITS != 0)			\
	      {						\
        do						\
		          {						\
        __x = (*__p + 1) & GMP_NUMB_MASK;	\
        *__p++ = __x;				\
		          }						\
				          while (__x == 0);				\
	      }						\
	        }							\
    } while (0)
#endif
#ifndef mpn_decr_u
#define mpn_decr_u(p,incr)				\
  do {							\
    mp_limb_t __x;					\
    mp_ptr __p = (p);					\
    if (__builtin_constant_p (incr) && (incr) == 1)	\
	      {							\
    do						\
	      {						\
        __x = *__p;					\
        *__p++ = (__x - 1) & GMP_NUMB_MASK;		\
	      }						\
		      while (__x == 0);				\
	      }							\
	    else						\
	        {							\
    __x = *__p - (incr);				\
    *__p++ = __x & GMP_NUMB_MASK;			\
    if (__x >> GMP_NUMB_BITS != 0)			\
	      {						\
        do						\
		          {						\
        __x = *__p;				\
        *__p++ = (__x - 1) & GMP_NUMB_MASK;	\
		          }						\
				          while (__x == 0);				\
	      }						\
	        }							\
    } while (0)
#endif
#endif
/***************************************************************/

/* Return non-zero if xp,xsize and yp,ysize overlap.
If xp+xsize<=yp there's no overlap, or if yp+ysize<=xp there's no
overlap.  If both these are false, there's an overlap. */
#define MPN_OVERLAP_P(xp, xsize, yp, ysize) \
  ((xp) + (xsize) > (yp) && (yp) + (ysize) > (xp))
#define MEM_OVERLAP_P(xp, xsize, yp, ysize)     \
  (   (char *) (xp) + (xsize) > (char *) (yp)   \
   && (char *) (yp) + (ysize) > (char *) (xp))

/* Return non-zero if xp,xsize and yp,ysize are either identical or not
overlapping.  Return zero if they're partially overlapping. */
#define MPN_SAME_OR_SEPARATE_P(xp, yp, size)    \
  MPN_SAME_OR_SEPARATE2_P(xp, size, yp, size)
#define MPN_SAME_OR_SEPARATE2_P(xp, xsize, yp, ysize)           \
  ((xp) == (yp) || ! MPN_OVERLAP_P (xp, xsize, yp, ysize))

/* Return non-zero if dst,dsize and src,ssize are either identical or
overlapping in a way suitable for an incrementing/decrementing algorithm.
Return zero if they're partially overlapping in an unsuitable fashion. */
#define MPN_SAME_OR_INCR2_P(dst, dsize, src, ssize)             \
  ((dst) <= (src) || ! MPN_OVERLAP_P (dst, dsize, src, ssize))
#define MPN_SAME_OR_INCR_P(dst, src, size)      \
  MPN_SAME_OR_INCR2_P(dst, size, src, size)
#define MPN_SAME_OR_DECR2_P(dst, dsize, src, ssize)             \
  ((dst) >= (src) || ! MPN_OVERLAP_P (dst, dsize, src, ssize))
#define MPN_SAME_OR_DECR_P(dst, src, size)      \
  MPN_SAME_OR_DECR2_P(dst, size, src, size)

/***************************************************************/

/* DIVEXACT_1_THRESHOLD is at what size to use mpn_divexact_1, as opposed to
plain mpn_divrem_1.  Likewise MODEXACT_1_ODD_THRESHOLD for
mpn_modexact_1_odd against plain mpn_mod_1.  On most CPUs divexact and
modexact are faster at all sizes, so the defaults are 0.  Those CPUs
where this is not right have a tuned threshold.  */
#ifndef DIVEXACT_1_THRESHOLD
#define DIVEXACT_1_THRESHOLD  0
#endif
#ifndef MODEXACT_1_ODD_THRESHOLD
#define MODEXACT_1_ODD_THRESHOLD  0
#endif

#define MPN_DIVREM_OR_DIVEXACT_1(dst, src, size, divisor)                     \
  do {                                                                        \
    if (BELOW_THRESHOLD (size, DIVEXACT_1_THRESHOLD))                         \
      ASSERT_NOCARRY (mpn_divrem_1 (dst, (mp_size_t) 0, src, size, divisor)); \
			    else                                                                      \
				        {                                                                       \
        ASSERT (mpn_mod_1 (src, size, divisor) == 0);                         \
        mpn_divexact_1 (dst, src, size, divisor);                             \
				        }                                                                       \
      } while (0)
/***************************************************************/
/* Use a library function for invert_limb, if available. */
#if ! defined (invert_limb) && HAVE_NATIVE_mpn_invert_limb
#define invert_limb(invxl,xl)           \
  do {                                  \
    (invxl) = mpn_invert_limb (xl);     \
    } while (0)
#endif

#ifndef invert_limb
#define invert_limb(invxl,xl)                   \
  do {                                          \
    mp_limb_t dummy;                            \
    ASSERT ((xl) != 0);                         \
    udiv_qrnnd (invxl, dummy, ~(xl), ~CNST_LIMB(0), xl);  \
    } while (0)
#endif

#define mpir_invert_pi1(dinv, d1, d0)					\
  do {									\
    mp_limb_t _v, _p, _t1, _t0, _mask;					\
    invert_limb (_v, d1);						\
    _p = (d1) * _v;							\
    _p += (d0);								\
    if (_p < (d0))							\
	      {									\
	_v--;								\
	_mask = -(mp_limb_t) (_p >= (d1));				\
	_p -= (d1);							\
	_v += _mask;							\
	_p -= _mask & (d1);						\
	      }									\
    umul_ppmm (_t1, _t0, d0, _v);					\
    _p += _t1;								\
    if (_p < _t1)							\
	      {									\
	_v--;								\
	if (UNLIKELY (_p >= (d1)))					\
		  {								\
	    if (_p > (d1) || _t0 >= (d0))				\
	      _v--;							\
		  }								\
	      }									\
    dinv = _v;							\
    } while (0)

/* For compatibility with GMP only */
#define invert_pi1(dinv, d1, d0)				\
   mpir_invert_pi1((dinv).inv32, d1, d0)

/* Compute quotient the quotient and remainder for n / d. Requires d
>= B^2 / 2 and n < d B. di is the inverse

floor ((B^3 - 1) / (d0 + d1 B)) - B.

NOTE: Output variables are updated multiple times. Only some inputs
and outputs may overlap.
*/
#define udiv_qr_3by2(q, r1, r0, n2, n1, n0, d1, d0, dinv)		\
  do {									\
    mp_limb_t _q0, _t1, _t0;					\
    umul_ppmm ((q), _q0, (n2), (dinv));					\
    add_ssaaaa ((q), _q0, (q), _q0, (n2), (n1));			\
                                    \
    /* Compute the two most significant limbs of n - q'd */		\
    (r1) = (n1) - (d1) * (q);						\
    sub_ddmmss ((r1), (r0), (r1), (n0), (d1), (d0));			\
    umul_ppmm (_t1, _t0, (d0), (q));					\
    sub_ddmmss ((r1), (r0), (r1), (r0), _t1, _t0);			\
    (q)++;								\
                                    \
    /* Conditionally adjust q and the remainders */			\
    if ((r1) >= _q0) {				\
       (q)--;							\
       add_ssaaaa ((r1), (r0), (r1), (r0), (d1), (d0));	} \
    if (UNLIKELY ((r1) >= (d1)))					\
	      {									\
    if ((r1) > (d1) || (r0) >= (d0))				\
	      {								\
        (q)++;							\
        sub_ddmmss ((r1), (r0), (r1), (r0), (d1), (d0));		\
	      }								\
	      }									\
    } while (0)

#ifndef udiv_qrnnd_preinv
#define udiv_qrnnd_preinv udiv_qrnnd_preinv2
#endif

/* Divide the two-limb number in (NH,,NL) by D, with DI being the largest
limb not larger than (2**(2*BITS_PER_MP_LIMB))/D - (2**BITS_PER_MP_LIMB).
If this would yield overflow, DI should be the largest possible number
(i.e., only ones).  For correct operation, the most significant bit of D
has to be set.  Put the quotient in Q and the remainder in R.  */
#define udiv_qrnnd_preinv1(q, r, nh, nl, d, di)				\
  do {									\
    mp_limb_t _q, _ql, _r;						\
    mp_limb_t _xh, _xl;							\
    ASSERT ((d) != 0);							\
    umul_ppmm (_q, _ql, (nh), (di));					\
    _q += (nh);	/* Compensate, di is 2**GMP_LIMB_BITS too small */	\
    umul_ppmm (_xh, _xl, _q, (d));					\
    sub_ddmmss (_xh, _r, (nh), (nl), _xh, _xl);				\
    if (_xh != 0)							\
	      {									\
    sub_ddmmss (_xh, _r, _xh, _r, 0, (d));				\
    _q += 1;							\
    if (_xh != 0)							\
	      {								\
        _r -= (d);							\
        _q += 1;							\
	      }								\
	      }									\
    if (_r >= (d))							\
	      {									\
    _r -= (d);							\
    _q += 1;							\
	      }									\
    (r) = _r;								\
    (q) = _q;								\
    } while (0)

/* Like udiv_qrnnd_preinv, but branch-free. */
#define udiv_qrnnd_preinv2(q, r, nh, nl, d, di)				\
  do {									\
    mp_limb_t _n2, _n10, _nmask, _nadj, _q1;				\
    mp_limb_t _xh, _xl;							\
    _n2 = (nh);								\
    _n10 = (nl);							\
    _nmask = LIMB_HIGHBIT_TO_MASK (_n10);				\
    _nadj = _n10 + (_nmask & (d));					\
    umul_ppmm (_xh, _xl, di, _n2 - _nmask);				\
    add_ssaaaa (_xh, _xl, _xh, _xl, _n2, _nadj);			\
    _q1 = ~_xh;								\
    umul_ppmm (_xh, _xl, _q1, d);					\
    add_ssaaaa (_xh, _xl, _xh, _xl, nh, nl);				\
    _xh -= (d);					/* xh = 0 or -1 */	\
    (r) = _xl + ((d) & _xh);						\
    (q) = _xh - _q1;							\
    } while (0)

/* Like udiv_qrnnd_preinv2, but for for any value D.  DNORM is D shifted left
so that its most significant bit is set.  LGUP is ceil(log2(D)).  */
#define udiv_qrnnd_preinv2gen(q, r, nh, nl, d, di, dnorm, lgup) \
  do {									\
    mp_limb_t _n2, _n10, _nmask, _nadj, _q1;				\
    mp_limb_t _xh, _xl;							\
    _n2 = ((nh) << (BITS_PER_MP_LIMB - (lgup))) + ((nl) >> 1 >> (l - 1));\
    _n10 = (nl) << (BITS_PER_MP_LIMB - (lgup));				\
    _nmask = LIMB_HIGHBIT_TO_MASK (_n10);				\
    _nadj = _n10 + (_nmask & (dnorm));					\
    umul_ppmm (_xh, _xl, di, _n2 - _nmask);				\
    add_ssaaaa (_xh, _xl, _xh, _xl, _n2, _nadj);			\
    _q1 = ~_xh;								\
    umul_ppmm (_xh, _xl, _q1, d);					\
    add_ssaaaa (_xh, _xl, _xh, _xl, nh, nl);				\
    _xh -= (d);								\
    (r) = _xl + ((d) & _xh);						\
    (q) = _xh - _q1;							\
    } while (0)

/***************************************************************/
#endif//__GMPN_ADD_H__
