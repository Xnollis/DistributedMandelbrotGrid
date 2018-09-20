#ifndef __GMP_IMPL_H__
#define __GMP_IMPL_H__
#include "config.h"



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
/////////////////////////////////////////////////////////////
#define ASSERT(expr)


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
#endif//__GMP_IMPL_H__