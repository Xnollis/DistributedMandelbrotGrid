#ifndef __GMP_IMPL_H__
#define __GMP_IMPL_H__
#include "config.h"

enum
{
    GMP_ERROR_NONE = 0,
    GMP_ERROR_UNSUPPORTED_ARGUMENT = 1,
    GMP_ERROR_DIVISION_BY_ZERO = 2,
    GMP_ERROR_SQRT_OF_NEGATIVE = 4,
    GMP_ERROR_INVALID_ARGUMENT = 8
};
__GMP_DECLSPEC_G_VALUE extern int __gmp_junk;
__GMP_DECLSPEC_G_VALUE extern const int __gmp_0;
__GMP_DECLSPEC void __gmp_exception  (int) ;
__GMP_DECLSPEC void __gmp_divide_by_zero  (void) ;
__GMP_DECLSPEC void __gmp_sqrt_of_negative  (void) ;
__GMP_DECLSPEC void __gmp_invalid_operation  (void) ;
#define GMP_ERROR(code)   __gmp_exception (code)
#define DIVIDE_BY_ZERO    __gmp_divide_by_zero ()
#define SQRT_OF_NEGATIVE  __gmp_sqrt_of_negative ()

#ifdef __cplusplus
#define __GMP_CAST(type, expr)  (static_cast<type> (expr))
#else
#define __GMP_CAST(type, expr)  ((type) (expr))
#endif


#define __GMP_ABS(x)   ((x) >= 0 ? (x) : -(x))
#define __GMP_MAX(h,i) ((h) > (i) ? (h) : (i))

/* __GMP_USHRT_MAX is not "~ (unsigned short) 0" because short is promoted
 to int by "~".  */
#define __GMP_UINT_MAX   (~ (unsigned) 0)
#define __GMP_ULONG_MAX  (~ (unsigned long) 0)
#define __GMP_USHRT_MAX  ((unsigned short) ~0)


/* Test for gcc >= maj.min, as per __GNUC_PREREQ in glibc */
#if defined (__GNUC__) && defined (__GNUC_MINOR__)
#define __GMP_GNUC_PREREQ(maj, min) \
  ((__GNUC__ << 16) + __GNUC_MINOR__ >= ((maj) << 16) + (min))
#else
#define __GMP_GNUC_PREREQ(maj, min)  0
#endif
/* __builtin_expect is in gcc 3.0, and not in 2.95. */
#if __GMP_GNUC_PREREQ (3,0)
#define __GMP_LIKELY(cond)    __builtin_expect ((cond) != 0, 1)
#define __GMP_UNLIKELY(cond)  __builtin_expect ((cond) != 0, 0)
#else
#define __GMP_LIKELY(cond)    (cond)
#define __GMP_UNLIKELY(cond)  (cond)
#endif

#ifndef WANT_TMP_DEBUG  /* for TMP_ALLOC_LIMBS_2 and others */
#define WANT_TMP_DEBUG 0
#endif



/* if not provided by gmp-mparam.h */
#ifndef BYTES_PER_MP_LIMB
#define BYTES_PER_MP_LIMB  SIZEOF_MP_LIMB_T
#endif
#ifndef BITS_PER_MP_LIMB
#define BITS_PER_MP_LIMB  (8 * SIZEOF_MP_LIMB_T)
#endif

#define BITS_PER_ULONG   (8 * SIZEOF_UNSIGNED_LONG)
#ifdef HAVE_STDINT_H
#define BITS_PER_UINTMAX (8 * SIZEOF_UINTMAX_T)
#endif
#ifdef _MSC_VER
#ifdef _WIN64
typedef unsigned int uint_least32_t;
#else
typedef unsigned int uint_least32_t;
#endif // _WIN64
#endif // _MSC_VER


/* gmp_uint_least32_t is an unsigned integer type with at least 32 bits. */
#if HAVE_UINT_LEAST32_T
typedef uint_least32_t      gmp_uint_least32_t;
#else
#if SIZEOF_UNSIGNED_SHORT >= 4
typedef unsigned short      gmp_uint_least32_t;
#else
#if SIZEOF_UNSIGNED >= 4
typedef unsigned            gmp_uint_least32_t;
#else
typedef unsigned long       gmp_uint_least32_t;
#endif
#endif
#endif
/////////////////////////////////////////////////////////////

/* The "short" defines are a bit different because shorts are promoted to
   ints by ~ or >> etc.

   #ifndef's are used since on some systems (HP?) header files other than
   limits.h setup these defines.  We could forcibly #undef in that case, but
   there seems no need to worry about that.  */

#ifndef ULONG_MAX
#define ULONG_MAX   __GMP_ULONG_MAX
#endif
#ifndef UINT_MAX
#define UINT_MAX    __GMP_UINT_MAX
#endif
#ifndef USHRT_MAX
#define USHRT_MAX   __GMP_USHRT_MAX
#endif
#define MP_LIMB_T_MAX      (~ (mp_limb_t) 0)

/* Must cast ULONG_MAX etc to unsigned long etc, since they might not be
   unsigned on a K&R compiler.  In particular the HP-UX 10 bundled K&R cc
   treats the plain decimal values in <limits.h> as signed.  */
#define ULONG_HIGHBIT      (ULONG_MAX ^ ((unsigned long) ULONG_MAX >> 1))
#define UINT_HIGHBIT       (UINT_MAX ^ ((unsigned) UINT_MAX >> 1))
#define USHRT_HIGHBIT      ((unsigned short) (USHRT_MAX ^ ((unsigned short) USHRT_MAX >> 1)))
#define GMP_LIMB_HIGHBIT  (MP_LIMB_T_MAX ^ (MP_LIMB_T_MAX >> 1))
#ifdef HAVE_STDINT_H
#define UINTMAX_HIGHBIT   (UINTMAX_MAX ^ (UINTMAX_MAX >> 1))
#endif

#ifndef LONG_MIN
#define LONG_MIN           ((long) ULONG_HIGHBIT)
#endif
#ifndef LONG_MAX
#define LONG_MAX           (-(LONG_MIN+1))
#endif

#ifndef INT_MIN
#define INT_MIN            ((int) UINT_HIGHBIT)
#endif
#ifndef INT_MAX
#define INT_MAX            (-(INT_MIN+1))
#endif

#ifndef SHRT_MIN
#define SHRT_MIN           ((short) USHRT_HIGHBIT)
#endif
#ifndef SHRT_MAX
#define SHRT_MAX           ((short) (-(SHRT_MIN+1)))
#endif
#if defined( _WIN64 )
#define MP_SIZE_T_MAX      _I64_MAX
#define MP_SIZE_T_MIN      _I64_MIN
#elif __GMP_MP_SIZE_T_INT
#define MP_SIZE_T_MAX      INT_MAX
#define MP_SIZE_T_MIN      INT_MIN
#else
#define MP_SIZE_T_MAX      LONG_MAX
#define MP_SIZE_T_MIN      LONG_MIN
#endif
//////////////////////////////////////////////////////////////////////////
/* These constants are generated by gen-fib.c header limbbits nailbits */
#if GMP_NUMB_BITS == 32
#define FIB_TABLE_LIMIT         47
#define FIB_TABLE_LUCNUM_LIMIT  46
#endif /* 32 bits */
#if GMP_NUMB_BITS == 64
#define FIB_TABLE_LIMIT         93
#define FIB_TABLE_LUCNUM_LIMIT  92
#endif /* 64 bits */

/* This constants are generated by gen-bases.c header limbbits nailbits */
#if GMP_NUMB_BITS == 32
#define MP_BASES_CHARS_PER_LIMB_10      9
#define MP_BASES_BIG_BASE_10            CNST_LIMB(0x3b9aca00)
#define MP_BASES_BIG_BASE_INVERTED_10   CNST_LIMB(0x12e0be82)
#define MP_BASES_NORMALIZATION_STEPS_10 2
#endif /* 32 bits */
#if GMP_NUMB_BITS == 64
#define MP_BASES_CHARS_PER_LIMB_10      19
#define MP_BASES_BIG_BASE_10            CNST_LIMB(0x8ac7230489e80000)
#define MP_BASES_BIG_BASE_INVERTED_10   CNST_LIMB(0xd83c94fb6d2ac34a)
#define MP_BASES_NORMALIZATION_STEPS_10 0
#endif /* 64 bits */
#define __GMP_HAVE_TOKEN_PASTE 1
#if defined _LONG_LONG_LIMB
#if __GMP_HAVE_TOKEN_PASTE
#define CNST_LIMB(C) ((mp_limb_t) C##LL)
#else
#define CNST_LIMB(C) ((mp_limb_t) C/**/LL)
#endif
#else /* not _LONG_LONG_LIMB */
#if __GMP_HAVE_TOKEN_PASTE
#define CNST_LIMB(C) ((mp_limb_t) C##L)
#else
#define CNST_LIMB(C) ((mp_limb_t) C/**/L)
#endif
#endif /* _LONG_LONG_LIMB */

/* This constants and defines are generated by gen-psqr limbbits nailbits */
#if GMP_LIMB_BITS == 32 && GMP_NAIL_BITS == 0
/* Non-zero bit indicates a quadratic residue mod 0x100.
This test identifies 82.81% as non-squares (212/256). */
static const mp_limb_t
sq_res_0x100[8] = {
	CNST_LIMB(0x2030213),
	CNST_LIMB(0x2020212),
	CNST_LIMB(0x2020213),
	CNST_LIMB(0x2020212),
	CNST_LIMB(0x2030212),
	CNST_LIMB(0x2020212),
	CNST_LIMB(0x2020212),
	CNST_LIMB(0x2020212),
};

/* 2^24-1 = 3^2 * 5 * 7 * 13 * 17 ... */
#define PERFSQR_MOD_BITS  25

/* This test identifies 95.66% as non-squares. */
#define PERFSQR_MOD_TEST(up, usize) \
  do {                              \
    mp_limb_t  r;                   \
    PERFSQR_MOD_34 (r, up, usize);  \
                                    \
    /* 73.33% */                    \
    PERFSQR_MOD_2 (r, CNST_LIMB(45), CNST_LIMB(0xfa4fa5), \
                   CNST_LIMB(0x920), CNST_LIMB(0x1a442481)); \
                                    \
    /* 47.06% */                    \
    PERFSQR_MOD_1 (r, CNST_LIMB(17), CNST_LIMB(0xf0f0f1), \
                   CNST_LIMB(0x1a317)); \
                                    \
    /* 46.15% */                    \
    PERFSQR_MOD_1 (r, CNST_LIMB(13), CNST_LIMB(0xec4ec5), \
                   CNST_LIMB(0x9e5)); \
                                    \
    /* 42.86% */                    \
    PERFSQR_MOD_1 (r, CNST_LIMB( 7), CNST_LIMB(0xdb6db7), \
                   CNST_LIMB(0x69)); \
    } while (0)

/* Grand total sq_res_0x100 and PERFSQR_MOD_TEST, 99.25% non-squares. */

/* helper for tests/mpz/t-perfsqr.c */
#define PERFSQR_DIVISORS  { 256, 45, 17, 13, 7, }

#elif GMP_LIMB_BITS == 64 && GMP_NAIL_BITS == 0

/* Non-zero bit indicates a quadratic residue mod 0x100.
This test identifies 82.81% as non-squares (212/256). */
static const mp_limb_t
sq_res_0x100[4] = {
	CNST_LIMB(0x202021202030213),
	CNST_LIMB(0x202021202020213),
	CNST_LIMB(0x202021202030212),
	CNST_LIMB(0x202021202020212),
};

/* 2^48-1 = 3^2 * 5 * 7 * 13 * 17 * 97 ... */
#define PERFSQR_MOD_BITS  49

/* This test identifies 97.81% as non-squares. */
#define PERFSQR_MOD_TEST(up, usize) \
  do {                              \
    mp_limb_t  r;                   \
    PERFSQR_MOD_34 (r, up, usize);  \
                                    \
    /* 69.23% */                    \
    PERFSQR_MOD_2 (r, CNST_LIMB(91), CNST_LIMB(0xfd2fd2fd2fd3), \
                   CNST_LIMB(0x2191240), CNST_LIMB(0x8850a206953820e1)); \
                                    \
    /* 68.24% */                    \
    PERFSQR_MOD_2 (r, CNST_LIMB(85), CNST_LIMB(0xfcfcfcfcfcfd), \
                   CNST_LIMB(0x82158), CNST_LIMB(0x10b48c4b4206a105)); \
                                    \
    /* 55.56% */                    \
    PERFSQR_MOD_1 (r, CNST_LIMB( 9), CNST_LIMB(0xe38e38e38e39), \
                   CNST_LIMB(0x93)); \
                                    \
    /* 49.48% */                    \
    PERFSQR_MOD_2 (r, CNST_LIMB(97), CNST_LIMB(0xfd5c5f02a3a1), \
                   CNST_LIMB(0x1eb628b47), CNST_LIMB(0x6067981b8b451b5f)); \
    } while (0)

/* Grand total sq_res_0x100 and PERFSQR_MOD_TEST, 99.62% non-squares. */

/* helper for tests/mpz/t-perfsqr.c */
#define PERFSQR_DIVISORS  { 256, 91, 85, 9, 97, }

//#pragma message(VAR_NAME_VALUE(GMP_LIMB_BITS))
#else

/* Some example here */
#error no data available for this limb size in perfsqr.h
#endif
//////////////////////////////////////////////////////////////////////////


/* mp_exp_t is the same as mp_size_t */
#if defined( _WIN64 )
#define MP_EXP_T_MAX   LONG_MAX
#define MP_EXP_T_MIN   LONG_MIN
#else
#define MP_EXP_T_MAX   MP_SIZE_T_MAX
#define MP_EXP_T_MIN   MP_SIZE_T_MIN
#endif

#define LONG_HIGHBIT       LONG_MIN
#define INT_HIGHBIT        INT_MIN
#define SHRT_HIGHBIT       SHRT_MIN

#define GMP_NUMB_HIGHBIT  (CNST_LIMB(1) << (GMP_NUMB_BITS-1))

#if GMP_NAIL_BITS == 0
#define GMP_NAIL_LOWBIT   CNST_LIMB(0)
#else
#define GMP_NAIL_LOWBIT   (CNST_LIMB(1) << GMP_NUMB_BITS)
#endif

/* Swap macros. */

#define MP_LIMB_T_SWAP(x, y)                    \
  do {                                          \
    mp_limb_t __mp_limb_t_swap__tmp = (x);      \
    (x) = (y);                                  \
    (y) = __mp_limb_t_swap__tmp;                \
      } while (0)
#define MP_SIZE_T_SWAP(x, y)                    \
  do {                                          \
    mp_size_t __mp_size_t_swap__tmp = (x);      \
    (x) = (y);                                  \
    (y) = __mp_size_t_swap__tmp;                \
      } while (0)

#define MP_PTR_SWAP(x, y)               \
  do {                                  \
    mp_ptr __mp_ptr_swap__tmp = (x);    \
    (x) = (y);                          \
    (y) = __mp_ptr_swap__tmp;           \
      } while (0)
#define MP_SRCPTR_SWAP(x, y)                    \
  do {                                          \
    mp_srcptr __mp_srcptr_swap__tmp = (x);      \
    (x) = (y);                                  \
    (y) = __mp_srcptr_swap__tmp;                \
      } while (0)

#define MPN_PTR_SWAP(xp,xs, yp,ys)      \
  do {                                  \
    MP_PTR_SWAP (xp, yp);               \
    MP_SIZE_T_SWAP (xs, ys);            \
      } while(0)
#define MPN_SRCPTR_SWAP(xp,xs, yp,ys)   \
  do {                                  \
    MP_SRCPTR_SWAP (xp, yp);            \
    MP_SIZE_T_SWAP (xs, ys);            \
      } while(0)

/////////////////////////////////////////////////////////////
#define ASSERT(expr)
#define ASSERT_ALWAYS(expr)
#define ASSERT_FAIL(expr)
#define ASSERT_LIMB(limb)       do {} while (0)
#define ASSERT_MPN(ptr, size)   do {} while (0)
#define ASSERT_CARRY(expr)     (expr)
#define ASSERT_NOCARRY(expr)   (expr)

/* Define stuff for longlong.h.  */
#if HAVE_ATTRIBUTE_MODE && defined (__GNUC__)
typedef unsigned int UQItype	__attribute__((mode(QI)));
typedef		 int SItype	__attribute__((mode(SI)));
typedef unsigned int USItype	__attribute__((mode(SI)));
typedef		 int DItype	__attribute__((mode(DI)));
typedef unsigned int UDItype	__attribute__((mode(DI)));
#else
typedef unsigned char UQItype;
typedef		 long SItype;
typedef unsigned long USItype;
#if HAVE_LONG_LONG
typedef	long long int DItype;
typedef unsigned long long int UDItype;
#else /* Assume `long' gives us a wide enough type.  Needed for hppa2.0w.  */
typedef long int DItype;
typedef unsigned long int UDItype;
#endif
#endif

typedef unsigned long long UWtype;
typedef unsigned int UHWtype;
#define W_TYPE_SIZE BITS_PER_MP_LIMB

#if HAVE_DOUBLE_IEEE_LITTLE_ENDIAN
#define _GMP_IEEE_FLOATS 1
union ieee_double_extract
{
  struct
    {
      gmp_uint_least32_t manl:32;
      gmp_uint_least32_t manh:20;
      gmp_uint_least32_t exp:11;
      gmp_uint_least32_t sig:1;
    } s;
  double d;
};
#endif

#if HAVE_DOUBLE_IEEE_BIG_ENDIAN
#define _GMP_IEEE_FLOATS 1
union ieee_double_extract
{
  struct
    {
      gmp_uint_least32_t sig:1;
      gmp_uint_least32_t exp:11;
      gmp_uint_least32_t manh:20;
      gmp_uint_least32_t manl:32;
    } s;
  double d;
};
#endif

/* DOUBLE_NAN_INF_ACTION executes code a_nan if x is a NaN, or executes
   a_inf if x is an infinity.  Both are considered unlikely values, for
   branch prediction.  */

#if _GMP_IEEE_FLOATS
#define DOUBLE_NAN_INF_ACTION(x, a_nan, a_inf)  \
  do {                                          \
    union ieee_double_extract  u;               \
    u.d = (x);                                  \
    if (UNLIKELY (u.s.exp == 0x7FF))            \
      {                                         \
        if (u.s.manl == 0 && u.s.manh == 0)     \
          { a_inf; }                            \
        else                                    \
          { a_nan; }                            \
      }                                         \
  } while (0)
#endif

#ifndef DOUBLE_NAN_INF_ACTION
/* Unknown format, try something generic.
   NaN should be "unordered", so x!=x.
   Inf should be bigger than DBL_MAX.  */
#define DOUBLE_NAN_INF_ACTION(x, a_nan, a_inf)                  \
  do {                                                          \
    {                                                           \
      if (UNLIKELY ((x) != (x)))                                \
        { a_nan; }                                              \
      else if (UNLIKELY ((x) > DBL_MAX || (x) < -DBL_MAX))      \
        { a_inf; }                                              \
    }                                                           \
  } while (0)
#endif



#define mpn_store(dst, n,val)			\
  do {						\
    ASSERT ((n) >= 0);				\
    if ((n) != 0)				\
      {						\
    mp_ptr __dst = (dst);			\
    mp_size_t __n = (n);			\
    do					\
      *__dst++ = val;			\
    while (--__n);				\
      }						\
  } while (0)

#define MPN_ZERO(dst,n)	mpn_store(dst,n,0)

#if ! defined (MPN_COPY_INCR) && HAVE_NATIVE_mpn_copyi
#define MPN_COPY_INCR(dst, src, size)                   \
do {                                                  \
ASSERT ((size) >= 0);                               \
ASSERT (MPN_SAME_OR_INCR_P (dst, src, size));       \
mpn_copyi (dst, src, size);                         \
} while (0)
#endif

/* Copy N limbs from SRC to DST incrementing, N==0 allowed.  */
#if ! defined (MPN_COPY_INCR)
#define MPN_COPY_INCR(dst, src, n)                      \
do {                                                  \
ASSERT ((n) >= 0);                                  \
ASSERT (MPN_SAME_OR_INCR_P (dst, src, n));          \
if ((n) != 0)                                       \
{                                                 \
mp_size_t __n = (n) - 1;                        \
mp_ptr __dst = (dst);                           \
mp_srcptr __src = (src);                        \
mp_limb_t __x;                                  \
__x = *__src++;                                 \
if (__n != 0)                                   \
{                                             \
do                                          \
{                                         \
*__dst++ = __x;                         \
__x = *__src++;                         \
}                                         \
while (--__n);                              \
}                                             \
*__dst++ = __x;                                 \
}                                                 \
} while (0)
#endif

#if ! defined (MPN_COPY_DECR) && HAVE_NATIVE_mpn_copyd
#define MPN_COPY_DECR(dst, src, size)                   \
do {                                                  \
ASSERT ((size) >= 0);                               \
ASSERT (MPN_SAME_OR_DECR_P (dst, src, size));       \
mpn_copyd (dst, src, size);                         \
} while (0)
#endif

/* Copy N limbs from SRC to DST decrementing, N==0 allowed.  */
#if ! defined (MPN_COPY_DECR)
#define MPN_COPY_DECR(dst, src, n)                      \
do {                                                  \
ASSERT ((n) >= 0);                                  \
ASSERT (MPN_SAME_OR_DECR_P (dst, src, n));          \
if ((n) != 0)                                       \
{                                                 \
mp_size_t __n = (n) - 1;                        \
mp_ptr __dst = (dst) + __n;                     \
mp_srcptr __src = (src) + __n;                  \
mp_limb_t __x;                                  \
__x = *__src--;                                 \
if (__n != 0)                                   \
{                                             \
do                                          \
{                                         \
*__dst-- = __x;                         \
__x = *__src--;                         \
}                                         \
while (--__n);                              \
}                                             \
*__dst-- = __x;                                 \
}                                                 \
} while (0)
#endif


#ifndef MPN_COPY
#define MPN_COPY(d,s,n)                         \
do {                                          \
ASSERT (MPN_SAME_OR_SEPARATE_P (d, s, n));  \
MPN_COPY_INCR (d, s, n);                    \
} while (0)
#endif


/* Set {dst,size} to the limbs of {src,size} in reverse order. */
#define MPN_REVERSE(dst, src, size)                     \
do {                                                  \
mp_ptr     __dst = (dst);                           \
mp_size_t  __size = (size);                         \
mp_srcptr  __src = (src) + __size - 1;              \
mp_size_t  __i;                                     \
ASSERT ((size) >= 0);                               \
ASSERT (! MPN_OVERLAP_P (dst, size, src, size));    \
for (__i = 0; __i < __size; __i++)                  \
{                                                 \
*__dst = *__src;                                \
__dst++;                                        \
__src--;                                        \
}                                                 \
} while (0)



#ifndef TMP_DECL
#define TMP_ALLOC_LOCAL_BUF_SIZE (sizeof(mp_limb_t)*MPIR_MAX_LIMB_SIZE*4)
#define TMP_DECL char __tp_temp_buf_zzz[TMP_ALLOC_LOCAL_BUF_SIZE];size_t __tp_temp_buf_zzz_cnt=0;
#define TMP_ALLOC(n) TMP_ALLOC_FUNC(n,__tp_temp_buf_zzz,&__tp_temp_buf_zzz_cnt)
#define TMP_SALLOC(n)        TMP_ALLOC(n)
#define TMP_BALLOC(n)        TMP_ALLOC(n)
#define TMP_MARK
#define TMP_FREE
#define TMP_SDECL TMP_DECL
#define TMP_SMARK TMP_MARK
#define TMP_SFREE TMP_FREE
/* Allocating various types. */
#define TMP_ALLOC_TYPE(n,type)  ((type *) TMP_ALLOC ((n) * sizeof (type)))
#define TMP_SALLOC_TYPE(n,type) ((type *) TMP_SALLOC ((n) * sizeof (type)))
#define TMP_BALLOC_TYPE(n,type) ((type *) TMP_BALLOC ((n) * sizeof (type)))
#define TMP_ALLOC_LIMBS(n)      TMP_ALLOC_TYPE(n,mp_limb_t)
#define TMP_SALLOC_LIMBS(n)     TMP_SALLOC_TYPE(n,mp_limb_t)
#define TMP_BALLOC_LIMBS(n)     TMP_BALLOC_TYPE(n,mp_limb_t)
#define TMP_ALLOC_MP_PTRS(n)    TMP_ALLOC_TYPE(n,mp_ptr)
#define TMP_SALLOC_MP_PTRS(n)   TMP_SALLOC_TYPE(n,mp_ptr)
#define TMP_BALLOC_MP_PTRS(n)   TMP_BALLOC_TYPE(n,mp_ptr)
__GMP_DECLSPEC void *TMP_ALLOC_FUNC(size_t n, char pLocalBuf[TMP_ALLOC_LOCAL_BUF_SIZE], size_t *pCnt);
#endif


/* Enhancement: __gmp_allocate_func could have "__attribute__ ((malloc))",
but current gcc (3.0) doesn't seem to support that.  */
__GMP_DECLSPEC extern void * (*__gmp_allocate_func) (size_t);
__GMP_DECLSPEC extern void * (*__gmp_reallocate_func) (void *, size_t, size_t);
__GMP_DECLSPEC extern void(*__gmp_free_func) (void *, size_t);

__GMP_DECLSPEC void *__gmp_default_allocate (size_t);
__GMP_DECLSPEC void *__gmp_default_reallocate (void *, size_t, size_t);
__GMP_DECLSPEC void __gmp_default_free (void *, size_t);

#define __GMP_ALLOCATE_FUNC_TYPE(n,type) \
  ((type *) (*__gmp_allocate_func) ((n) * sizeof (type)))
#define __GMP_ALLOCATE_FUNC_LIMBS(n)   __GMP_ALLOCATE_FUNC_TYPE (n, mp_limb_t)

#define __GMP_REALLOCATE_FUNC_TYPE(p, old_size, new_size, type) \
  ((type *) (*__gmp_reallocate_func)                            \
   (p, (old_size) * sizeof (type), (new_size) * sizeof (type)))
#define __GMP_REALLOCATE_FUNC_LIMBS(p, old_size, new_size) \
  __GMP_REALLOCATE_FUNC_TYPE(p, old_size, new_size, mp_limb_t)
#define __GMP_FREE_FUNC_TYPE(p,n,type) (*__gmp_free_func) (p, (n) * sizeof (type))
#define __GMP_FREE_FUNC_LIMBS(p,n)     __GMP_FREE_FUNC_TYPE (p, n, mp_limb_t)


#endif//__GMP_IMPL_H__
