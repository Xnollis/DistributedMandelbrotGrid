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
mp_limb_t dummy;                            \
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

__GMP_DECLSPEC extern const unsigned char  modlimb_invert_table[128];

#define modlimb_invert(inv,n)                        \
do {                                    \
mp_limb_t  __n = (n);                        \
mp_limb_t  __inv;                            \
ASSERT ((__n & 1) == 1);                        \
\
__inv = modlimb_invert_table[(__n/2) & 0x7F]; /*  8 */        \
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


#endif//__GMPN_ADD_H__
