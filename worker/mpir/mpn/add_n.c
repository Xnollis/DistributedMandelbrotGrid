/* mpn_add_n -- Add equal length limb vectors.

Copyright 1992, 1993, 1994, 1996, 2000, 2002 Free Software Foundation, Inc.

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

#if GMP_NAIL_BITS == 0

__GMP_DECLSPEC
mp_limb_t
mpn_add_n (mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n)
{
  mp_limb_t ul, vl, sl, rl, cy, cy1, cy2;

  ASSERT (n >= 1);
  ASSERT (MPN_SAME_OR_SEPARATE_P (rp, up, n));
  ASSERT (MPN_SAME_OR_SEPARATE_P (rp, vp, n));

  cy = 0;
  do
    {
      ul = *up++;
      vl = *vp++;
      sl = ul + vl;
      cy1 = sl < ul;
      rl = sl + cy;
      cy2 = rl < sl;
      cy = cy1 | cy2;
      *rp++ = rl;
    }
  while (--n != 0);

  return cy;
}

#endif

#if GMP_NAIL_BITS >= 1

__GMP_DECLSPEC
mp_limb_t
mpn_add_n (mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n)
{
  mp_limb_t ul, vl, rl, cy;

  ASSERT (n >= 1);
  ASSERT (MPN_SAME_OR_SEPARATE_P (rp, up, n));
  ASSERT (MPN_SAME_OR_SEPARATE_P (rp, vp, n));

  cy = 0;
  do
    {
      ul = *up++;
      vl = *vp++;
      rl = ul + vl + cy;
      cy = rl >> GMP_NUMB_BITS;
      *rp++ = rl & GMP_NUMB_MASK;
    }
  while (--n != 0);

  return cy;
}

#endif

__GMP_DECLSPEC
mp_limb_t mpn_add (mp_ptr __gmp_wp, mp_srcptr __gmp_xp, mp_size_t __gmp_xsize, mp_srcptr __gmp_yp, mp_size_t __gmp_ysize)
{
    mp_limb_t  __gmp_c;
    __GMPN_ADD (__gmp_c, __gmp_wp, __gmp_xp, __gmp_xsize, __gmp_yp, __gmp_ysize);
    return __gmp_c;
}
__GMP_DECLSPEC
mp_limb_t mpn_add_1 (mp_ptr __gmp_dst, mp_srcptr __gmp_src, mp_size_t __gmp_size, mp_limb_t __gmp_n) __GMP_NOTHROW
{
    mp_limb_t  __gmp_c;
    __GMPN_ADD_1 (__gmp_c, __gmp_dst, __gmp_src, __gmp_size, __gmp_n);
    return __gmp_c;
}
__GMP_DECLSPEC
int mpn_cmp (mp_srcptr __gmp_xp, mp_srcptr __gmp_yp, mp_size_t __gmp_size) __GMP_NOTHROW
{
    int __gmp_result;
    __GMPN_CMP (__gmp_result, __gmp_xp, __gmp_yp, __gmp_size);
    return __gmp_result;
}
