/* mpn_mod_1(dividend_ptr, dividend_size, divisor_limb) --
   Divide (DIVIDEND_PTR,,DIVIDEND_SIZE) by DIVISOR_LIMB.
   Return the single-limb remainder.
   There are no constraints on the value of the divisor.

Copyright 1991, 1993, 1994, 1999, 2000, 2002 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include "mpir_inter_decl.h"

/* The size where udiv_qrnnd_preinv should be used rather than udiv_qrnnd,
   meaning the quotient size where that should happen, the quotient size
   being how many udiv divisions will be done.

   The default is to use preinv always, CPUs where this doesn't suit have
   tuned thresholds.  Note in particular that preinv should certainly be
   used if that's the only division available (USE_PREINV_ALWAYS).  */

#ifndef MOD_1_NORM_THRESHOLD
#define MOD_1_NORM_THRESHOLD  0
#endif
#ifndef MOD_1_UNNORM_THRESHOLD
#define MOD_1_UNNORM_THRESHOLD  0
#endif


/* The comments in mpn/generic/divrem_1.c apply here too.

   As noted in the algorithms section of the manual, the shifts in the loop
   for the unnorm case can be avoided by calculating r = a%(d*2^n), followed
   by a final (r*2^n)%(d*2^n).  In fact if it happens that a%(d*2^n) can
   skip a division where (a*2^n)%(d*2^n) can't then there's the same number
   of divide steps, though how often that happens depends on the assumed
   distributions of dividend and divisor.  In any case this idea is left to
   CPU specific implementations to consider.  */
#if 1
__GMP_DECLSPEC
mp_limb_t
mpn_mod_1(mp_srcptr up, mp_size_t un, mp_limb_t d)
{
  mp_size_t i;
  mp_limb_t n1, n0, r;
  mp_limb_t dummy;
  dummy = 0;
  ASSERT(un >= 0);
  ASSERT(d != 0);

  /* Botch: Should this be handled at all?  Rely on callers?
     But note un==0 is currently required by mpz/fdiv_r_ui.c and possibly
     other places.  */
  if (un == 0)
    return 0;

#if HAVE_NATIVE_mpn_divrem_euclidean_r_1
  return mpn_divrem_euclidean_r_1(up, un, d);
#endif

  d <<= GMP_NAIL_BITS;

  if ((d & GMP_LIMB_HIGHBIT) != 0)
  {
    /* High limb is initial remainder, possibly with one subtract of
	 d to get r<d.  */
    r = up[un - 1] << GMP_NAIL_BITS;
    if (r >= d)
      r -= d;
    r >>= GMP_NAIL_BITS;
    un--;
    if (un == 0)
      return r;

    if (BELOW_THRESHOLD(un, MOD_1_NORM_THRESHOLD))
    {
    plain:
      for (i = un - 1; i >= 0; i--)
      {
        n0 = up[i] << GMP_NAIL_BITS;
        udiv_qrnnd(dummy, r, r, n0, d);
        r >>= GMP_NAIL_BITS;
      }
      return r;
    }
    else
    {
      mp_limb_t inv;
      invert_limb(inv, d);
      for (i = un - 1; i >= 0; i--)
      {
        n0 = up[i] << GMP_NAIL_BITS;
        udiv_qrnnd_preinv(dummy, r, r, n0, d, inv);
        r >>= GMP_NAIL_BITS;
      }
      return r;
    }
  }
  else
  {
    int norm;

    /* Skip a division if high < divisor.  Having the test here before
	 normalizing will still skip as often as possible.  */
    r = up[un - 1] << GMP_NAIL_BITS;
    if (r < d)
    {
      r >>= GMP_NAIL_BITS;
      un--;
      if (un == 0)
        return r;
    }
    else
      r = 0;

    /* If udiv_qrnnd doesn't need a normalized divisor, can use the simple
	 code above. */
    if (!UDIV_NEEDS_NORMALIZATION && BELOW_THRESHOLD(un, MOD_1_UNNORM_THRESHOLD))
      goto plain;

    count_leading_zeros(norm, d);
    d <<= norm;

    n1 = up[un - 1] << GMP_NAIL_BITS;
    r = (r << norm) | (n1 >> (GMP_LIMB_BITS - norm));

    if (UDIV_NEEDS_NORMALIZATION && BELOW_THRESHOLD(un, MOD_1_UNNORM_THRESHOLD))
    {
      for (i = un - 2; i >= 0; i--)
      {
        n0 = up[i] << GMP_NAIL_BITS;
        udiv_qrnnd(dummy, r, r,
                   (n1 << norm) | (n0 >> (GMP_NUMB_BITS - norm)),
                   d);
        r >>= GMP_NAIL_BITS;
        n1 = n0;
      }
      udiv_qrnnd(dummy, r, r, n1 << norm, d);
      r >>= GMP_NAIL_BITS;
      return r >> norm;
    }
    else
    {
      mp_limb_t inv;
      invert_limb(inv, d);

      for (i = un - 2; i >= 0; i--)
      {
        n0 = up[i] << GMP_NAIL_BITS;
        udiv_qrnnd_preinv(dummy, r, r,
                          (n1 << norm) | (n0 >> (GMP_NUMB_BITS - norm)),
                          d, inv);
        r >>= GMP_NAIL_BITS;
        n1 = n0;
      }
      udiv_qrnnd_preinv(dummy, r, r, n1 << norm, d, inv);
      r >>= GMP_NAIL_BITS;
      REMOVE_WARNINGS_OF_LOCAL_VAR(dummy);
      return r >> norm;
    }
  }
}
#endif
/* in each round we remove one limb from the body, i.e. k = 1 */
__GMP_DECLSPEC
void mpn_mod_1_1(mp_ptr rem, mp_srcptr xp, mp_size_t xn, mp_srcptr db)
{
    mp_limb_t h, l, sh, sl;
    mp_size_t j;
    
    ASSERT(xn >= 3);
    ASSERT_MPN(xp, xn);
    ASSERT_LIMB(db[0]);
    ASSERT_LIMB(db[1]);
    
    h = xp[xn - 1];
    l = xp[xn - 2];
    
    for (j = xn - 3; j >= 0; j--)
    {
        umul_ppmm(sh, sl, l, db[0]);
        add_ssaaaa(sh, sl, sh, sl, 0, xp[j]);
        umul_ppmm(h, l, h, db[1]);
        add_ssaaaa(h, l, h, l, sh, sl);
    }
    
    umul_ppmm(sh, sl, h, db[0]);
    add_ssaaaa(sh, sl, sh, sl, 0, l);
    
    rem[0] = sl;
    rem[1] = sh;
}

/* in each round we remove one limb from the body, i.e. k = 1 */
__GMP_DECLSPEC
void mpn_mod_1_2(mp_ptr rem, mp_srcptr xp, mp_size_t xn, mp_srcptr db)
{
    mp_limb_t h, l, sh, sl, th, tl;
    mp_size_t j;
    
    ASSERT(xn >= 4);
    ASSERT_MPN(xp, xn);
    ASSERT_LIMB(db[0]);
    ASSERT_LIMB(db[1]);
    ASSERT_LIMB(db[2]);
    
    tl = xp[xn - 2];
    th = xp[xn - 1];
    
    for (j = xn - 4; j >= 0; j -= 2)
    {
        umul_ppmm(sh, sl, xp[j + 1], db[0]);
        add_ssaaaa(sh, sl, sh, sl, 0, xp[j]);
        umul_ppmm(h, l, tl, db[1]);
        add_ssaaaa(sh, sl, sh, sl, h, l);
        umul_ppmm(th, tl, th, db[2]);
        add_ssaaaa(th, tl, th, tl, sh, sl);
    }
    
    if (j > -2) /* we have at least three limbs to do i.e. xp[0], ..., tl, th */
    {
        umul_ppmm(sh, sl, tl, db[0]);
        add_ssaaaa(sh, sl, sh, sl, 0, xp[0]);
        umul_ppmm(th, tl, th, db[1]);
        add_ssaaaa(th, tl, th, tl, sh, sl);
    }
    
    umul_ppmm(h, l, th, db[0]);
    add_ssaaaa(h, l, h, l, 0, tl);
    
    rem[0] = l;
    rem[1] = h;
}

/* in each round we remove one limb from the body, i.e. k = 1 */
__GMP_DECLSPEC
void mpn_mod_1_3(mp_ptr rem, mp_srcptr xp, mp_size_t xn, mp_srcptr db)
{
    mp_limb_t h, l, sh, sl, th, tl;
    mp_size_t j, jj;
    
    ASSERT(xn >= 5);
    ASSERT_MPN(xp, xn);
    ASSERT_LIMB(db[0]);
    ASSERT_LIMB(db[1]);
    ASSERT_LIMB(db[2]);
    ASSERT_LIMB(db[3]);
    
    tl = xp[xn - 2];
    th = xp[xn - 1];
    
    for (j = xn - 5; j >= 0; j -= 3)
    {
        umul_ppmm(sh, sl, xp[j + 1], db[0]);
        add_ssaaaa(sh, sl, sh, sl, 0, xp[j]);
        umul_ppmm(h, l, xp[j + 2], db[1]);
        add_ssaaaa(sh, sl, sh, sl, h, l);
        umul_ppmm(h, l, tl, db[2]);
        add_ssaaaa(sh, sl, sh, sl, h, l);
        umul_ppmm(th, tl, th, db[3]);
        add_ssaaaa(th, tl, th, tl, sh, sl);
    }
    
    if (j > -3) /* we have at least three limbs to do, i.e. xp[0], ..., tl, th */
    {
        sh = 0;
        sl = xp[0];
        jj = 1;
        
        if (j == -1)
        {
            umul_ppmm(sh, sl, xp[1], db[0]);
            add_ssaaaa(sh, sl, sh, sl, 0, xp[0]);
            jj = 2;
        }
        
        umul_ppmm(h, l, tl, db[jj - 1]);
        add_ssaaaa(sh, sl, sh, sl, h, l);
        umul_ppmm(th, tl, th, db[jj]);
        add_ssaaaa(th, tl, th, tl, sh, sl);
    }
    
    umul_ppmm(h, l, th, db[0]);
    add_ssaaaa(h, l, h, l, 0, tl);
    
    rem[0] = l;
    rem[1] = h;
}
