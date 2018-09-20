#ifndef __MPIR_H__
#define __MPIR_H__

#include <stdio.h>
#include <stdlib.h>
#if defined (__cplusplus)
extern "C" {
#endif // __cplusplus
#include "gmp_impl.h"
//compile language defines
#define __mpir_const const
//compile data type defines
#define _LONG_LONG_LIMB	1  //x64 SYSTEM sign
typedef unsigned long long int	mp_limb_t;
typedef long long int		mp_limb_signed_t;

#define __GMP_BITS_PER_MP_LIMB             64
#define GMP_LIMB_BITS                      64
#define GMP_NAIL_BITS                      0
#define BITS_PER_UI         BITS_PER_MP_LIMB
typedef mp_limb_t           mpir_ui;
typedef mp_limb_signed_t    mpir_si;
typedef mpir_ui             mp_bitcnt_t;
typedef long long int	mp_size_t;
typedef long int		mp_exp_t;
typedef mp_limb_t *		mp_ptr;
typedef __mpir_const mp_limb_t *	mp_srcptr;
/* pre-inverse types for truncating division and modulo */
typedef struct {mp_limb_t inv32;} gmp_pi1_t;
typedef struct {mp_limb_t inv21, inv32, inv53;} gmp_pi2_t;

#define GMP_UI_MAX          ((mpir_ui)(~(mpir_ui)0))
#define GMP_UI_HIBIT        (GMP_UI_MAX ^ (GMP_UI_MAX >> 1))
#define GMP_SI_MAX          ((mpir_si)(GMP_UI_MAX ^ GMP_UI_HIBIT))
#define GMP_SI_MIN          ((mpir_si)GMP_UI_HIBIT)
#define __GMP_BITCNT_MAX    (~(mp_bitcnt_t)0)
#define GMP_NUMB_BITS     (GMP_LIMB_BITS - GMP_NAIL_BITS)
#define GMP_NUMB_MASK     ((~ __GMP_CAST (mp_limb_t, 0)) >> GMP_NAIL_BITS)
#define GMP_NUMB_MAX      GMP_NUMB_MASK
#define GMP_NAIL_MASK     (~ GMP_NUMB_MASK)
/* Use (4.0 * ...) instead of (2.0 * ...) to work around buggy compilers
   that don't convert ulong->double correctly (eg. SunOS 4 native cc).  */
#define MP_BASE_AS_DOUBLE (4.0 * ((mp_limb_t) 1 << (GMP_NUMB_BITS - 2)))
/* Maximum number of limbs it will take to store any `double'.
   We assume doubles have 53 mantissam bits.  */
#define LIMBS_PER_DOUBLE ((53 + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS + 1)

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#undef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#undef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#define numberof(x)  (sizeof (x) / sizeof ((x)[0]))

/* Field access macros.  */
#define SIZ(x) ((x)->_mp_size)
#define ABSIZ(x) ABS (SIZ (x))
#define PTR(x) ((x)->_mp_d)
#define LIMBS(x) ((x)->_mp_d)
#define EXP(x) ((x)->_mp_exp)
#define PREC(x) ((x)->_mp_prec)
#define ALLOC(x) ((x)->_mp_alloc)

#define MPN_CMP(result, xp, yp, size)  __GMPN_CMP(result, xp, yp, size)
#define LIKELY(cond)                   __GMP_LIKELY(cond)
#define UNLIKELY(cond)                 __GMP_UNLIKELY(cond)
//////////////////////////////////////
// mpf Struct Define
#define MPIR_MAX_LIMB_SIZE 300
typedef struct
{
  int _mp_prec;			/* Max precision, in number of `mp_limb_t's.
				   Set by mpf_init and modified by
				   mpf_set_prec.  The area pointed to by the
				   _mp_d field contains `prec' + 1 limbs.  */
  int _mp_size;			/* abs(_mp_size) is the number of limbs the
				   last field points to.  If _mp_size is
				   negative this is a negative number.  */
  mp_exp_t _mp_exp;		/* Exponent, in the base of `mp_limb_t'.  */
  mp_limb_t _mp_d[MPIR_MAX_LIMB_SIZE];		/* Pointer to the limbs.  */
} __mpf_struct;
typedef __mpf_struct *mpf_ptr;
typedef __mpf_struct mpf_t[1];


void mpf_init (mpf_ptr r);
void mpf_clear (mpf_ptr m);
void mpf_set_d (mpf_ptr r, double d);

int IsCUDA_Supported(int bPrintInfoToConsole);
#if defined (__cplusplus)
}
#endif // __cplusplus
//////////////////////////////////
///End of file
#endif // !__MPIR_H__