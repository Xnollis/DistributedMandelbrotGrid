/* mpn_mul_n and helper function -- Multiply/square natural numbers.

   THE HELPER FUNCTIONS IN THIS FILE (meaning everything except mpn_mul_n)
   ARE INTERNAL FUNCTIONS WITH MUTABLE INTERFACES.  IT IS ONLY SAFE TO REACH
   THEM THROUGH DOCUMENTED INTERFACES.  IN FACT, IT IS ALMOST GUARANTEED
   THAT THEY'LL CHANGE OR DISAPPEAR IN A FUTURE GNU MP RELEASE.


Copyright 1991, 1993, 1994, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003,
2005, Free Software Foundation, Inc.

Copyright 2009 Jason Moxham
Copyright 2009 William Hart
Copyright 2011 The Code Cavern

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


#if ! HAVE_NATIVE_mpn_karasub && HAVE_NATIVE_mpn_addsub_n

static void	mpn_karasub(mp_ptr rp, mp_ptr tp, mp_size_t n)
{
   mp_size_t n2, n3;
   mp_limb_t c1 = 0, c2, c3, top[2];

   n2 = n>>1;
   n3 = n - n2;

   c2 = mpn_addsub_n(tp, rp, rp + 2*n2, tp, 2*n2);
   c3 = mpn_add_n(rp + n2, rp + n2, tp, 2*n2);

   top[1] = rp[2*n2 + 2*n3 - 1];
   top[0] = rp[2*n2 + 2*n3 - 2];

   mpn_incr_u(rp + 3*n2, c3);

   if (c2 == 1) mpn_incr_u(rp + 3*n2, 1);
   if (c2 == -1) mpn_decr_u(rp + 3*n2, 1);

   if (n2 == n3) 
      return;

   c1=mpn_sub_n(rp + 3*n2, rp + 3*n2, tp + 2*n2, 2);
   c2=mpn_add_n(rp + 3*n2, rp + 3*n2, top, 2);

   if(c2 == 1 && c1 == 0) mpn_incr_u(rp + 3*n2 + 2, 1);
   if(c2 == 0 && c1 == 1) mpn_decr_u(rp + 3*n2 + 2, 1);
}
#endif

#if ! HAVE_NATIVE_mpn_karaadd && HAVE_NATIVE_mpn_addadd_n

static void	mpn_karaadd(mp_ptr rp, mp_ptr tp, mp_size_t n)
{
   mp_size_t n2, n3;
   mp_limb_t c1 = 0, c2, c3;

   n2 = n>>1;
   n3 = n - n2;

   c2 = mpn_addadd_n(tp, rp, rp + 2*n2, tp, 2*n2);

   if (n3 != n2) c1 = mpn_add_n(tp + 2*n2, rp + 4*n2, tp + 2*n2, 2);

   c3 = mpn_add_n(rp + n2, rp + n2, tp, 2*n3);

   mpn_incr_u(rp + n2 + 2*n3, c3 + c1);
   mpn_incr_u(rp + n2 + 2*n2, c2);
}
#endif

#if ! HAVE_NATIVE_mpn_karasub && ! HAVE_NATIVE_mpn_addsub_n

static void	mpn_karasub(mp_ptr rp, mp_ptr tp, mp_size_t n)
{
   mp_size_t n2, n3;
   mp_limb_t c1, c2, c3, top[2];

   n2 = n>>1;
   n3 = n - n2;

   c1 = mpn_sub_n(tp, rp + 2*n2, tp, 2*n2);
   c2 = mpn_add_n(tp, tp, rp, 2*n2);
   c3 = mpn_add_n(rp + n2, rp + n2, tp, 2*n2);

   top[1] = rp[2*n2 + 2*n3 - 1];
   top[0] = rp[2*n2 + 2*n3 - 2];

   mpn_incr_u(rp + 3*n2, c3);
   mpn_incr_u(rp + 3*n2, c2);
   mpn_decr_u(rp + 3*n2, c1);

   if(n2 == n3)
      return;

   c1 = mpn_sub_n(rp + 3*n2, rp + 3*n2, tp + 2*n2, 2);
   c2 = mpn_add_n(rp + 3*n2, rp + 3*n2, top, 2);

   if(c2 == 1 && c1 == 0) mpn_incr_u(rp + 3*n2 + 2, 1);
   if(c2 == 0 && c1 == 1) mpn_decr_u(rp + 3*n2 + 2, 1);
}
#endif

#if ! HAVE_NATIVE_mpn_karaadd && ! HAVE_NATIVE_mpn_addadd_n

static void	mpn_karaadd(mp_ptr rp, mp_ptr tp, mp_size_t n)
{
   mp_size_t n2, n3;
   mp_limb_t c1, c2, c3;

   n2 = n>>1;
   n3 = n - n2;

   c1 = mpn_add_n(tp, rp + 2*n2, tp, 2*n3);
   c2 = mpn_add_n(tp, tp, rp, 2*n2);
   c3 = mpn_add_n(rp + n2, rp + n2, tp, 2*n3);
   
   mpn_incr_u(rp + n2 + 2*n3, c3 + c1);
   mpn_incr_u(rp + n2 + 2*n2, c2);
}
#endif

/* (rp, 2n) = (xp, n)*(yp, n) with temp space (tp, 2*n + C) */
void mpn_kara_mul_n(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n, mp_ptr tp)
{
   mp_size_t n2, n3;
   mp_srcptr xl, yl, xh, yh;
   mp_ptr dx, dy;
   int suboradd;
   mp_limb_t c;

   n2 = n>>1;
   suboradd = -1;
   xl = xp;
   xh = xp + n2;
   yl = yp;
   yh = yp + n2;
   n3 = n - n2;
   dx = rp + 2*n2;
   dy = dx + n3;

   if ((n&1) == 0)
   {
      if (mpn_cmp(xh, xl, n2) >= 0)
         mpn_sub_n(dx, xh, xl, n2);
      else
      {
         mpn_sub_n(dx, xl, xh, n2);
         suboradd = -suboradd;
      }
      
      if (mpn_cmp(yh, yl, n2) >= 0)
         mpn_sub_n(dy, yh, yl, n2);
      else
      {
         mpn_sub_n(dy, yl, yh, n2);
         suboradd = -suboradd;
      }
   }
   else
   {
      if (xh[n2] !=0 || mpn_cmp(xh, xl, n2) >= 0)
      {
         c = mpn_sub_n(dx, xh, xl, n2);
         dx[n2] = xh[n2] - c;
      }
      else
      {
         mpn_sub_n(dx, xl, xh, n2);
         dx[n2] = 0;
         suboradd = -suboradd;
      }
      
      if (yh[n2] != 0 || mpn_cmp(yh, yl, n2) >= 0)
      {
         c = mpn_sub_n(dy, yh, yl, n2);
         dy[n2] = yh[n2] - c;
      }
      else
      {
         mpn_sub_n(dy, yl, yh, n2);
         dy[n2] = 0;
         suboradd = -suboradd;
      }
   }

   if (BELOW_THRESHOLD(n3, MUL_KARATSUBA_THRESHOLD))
   {
      mpn_mul_basecase(rp, xl, n2, yl, n2);
      mpn_mul_basecase(tp, dx, n3, dy, n3);
      mpn_mul_basecase(rp + 2*n2, xh, n3, yh, n3);
   }
   else
   {
      mpn_kara_mul_n(rp, xl, yl, n2, tp + 2*n3);
      mpn_kara_mul_n(tp, dx, dy, n3, tp + 2*n3);   
      mpn_kara_mul_n(rp + 2*n2, xh, yh, n3, tp + 2*n3);
   }

   if (suboradd == -1)
      mpn_karasub(rp, tp, n);
   else
      mpn_karaadd(rp, tp, n);
}

/* (rp, 2n) = (xp, n)^2 with temp space (tp, 2*n + C) */
void mpn_kara_sqr_n(mp_ptr rp, mp_srcptr xp, mp_size_t n, mp_ptr tp)
{
   mp_size_t n2, n3;
   mp_srcptr xl, xh;
   mp_ptr dx;
   mp_limb_t c;

   n2 = n>>1;
   xl = xp;
   xh = xp + n2;
   n3 = n - n2;
   dx = rp + 2*n2;

   if ((n&1) == 0)
   {
      if (mpn_cmp(xh, xl, n2) >=0)
         mpn_sub_n(dx, xh, xl, n2);
      else
         mpn_sub_n(dx, xl, xh, n2);
   }
   else
   {
      if (xh[n2] != 0 || mpn_cmp(xh, xl, n2) >= 0)
      {
         c = mpn_sub_n(dx, xh, xl, n2);
         dx[n2] = xh[n2] - c;
      }
      else
      {
         mpn_sub_n(dx, xl, xh, n2);
         dx[n2] = 0;
      }
   }

   if (BELOW_THRESHOLD(n3, SQR_BASECASE_THRESHOLD))
   {
      mpn_mul_basecase(rp, xl, n2, xl, n2);
      mpn_mul_basecase(tp, dx, n3, dx, n3);
      mpn_mul_basecase(rp + 2*n2, xh, n3, xh, n3);
   }
   else if (BELOW_THRESHOLD(n3, SQR_KARATSUBA_THRESHOLD))
   {
      mpn_sqr_basecase(rp, xl, n2);
      mpn_sqr_basecase(tp, dx, n3);
      mpn_sqr_basecase(rp + 2*n2, xh, n3);
   }
   else
   {
      mpn_kara_sqr_n(rp, xl, n2, tp + 2*n3);
      mpn_kara_sqr_n(tp, dx, n3, tp + 2*n3);   
      mpn_kara_sqr_n(rp + 2*n2, xh, n3, tp + 2*n3);
   }

   mpn_karasub(rp, tp, n);
}

void
mpn_mul_n (mp_ptr p, mp_srcptr a, mp_srcptr b, mp_size_t n)
{
  ASSERT (n >= 1);
  ASSERT (! MPN_OVERLAP_P (p, 2 * n, a, n));
  ASSERT (! MPN_OVERLAP_P (p, 2 * n, b, n));

  if (BELOW_THRESHOLD (n, MUL_KARATSUBA_THRESHOLD))
    {
      mpn_mul_basecase (p, a, n, b, n);
    }
  else if (BELOW_THRESHOLD (n, MUL_TOOM3_THRESHOLD))
    {
      /* Allocate workspace of fixed size on stack: fast! */
      mp_limb_t ws[MPN_KARA_MUL_N_TSIZE (MUL_TOOM3_THRESHOLD_LIMIT-1)];
      ASSERT (MUL_TOOM3_THRESHOLD <= MUL_TOOM3_THRESHOLD_LIMIT);
      mpn_kara_mul_n (p, a, b, n, ws);
    }
  else if (BELOW_THRESHOLD (n, MUL_TOOM4_THRESHOLD))
    {
      mp_ptr ws;
      TMP_SDECL;
      TMP_SMARK;
      ws = TMP_SALLOC_LIMBS (MPN_TOOM3_MUL_N_TSIZE (n));
      mpn_toom3_mul_n (p, a, b, n, ws);
      TMP_SFREE;
    }
  else if (BELOW_THRESHOLD (n, MUL_TOOM8H_THRESHOLD))
    {
       mpn_toom4_mul_n (p, a, b, n);
    }
#if WANT_FFT || TUNE_PROGRAM_BUILD
  else if (BELOW_THRESHOLD (n, MUL_FFT_FULL_THRESHOLD))
    {
       mpn_toom8h_mul (p, a, n, b, n);
    }
#endif
  else
#if WANT_FFT || TUNE_PROGRAM_BUILD
    {
       mpn_mul_fft_main(p, a, n, b, n); 
    }
#else
    {
      /* Toom8 for large operands. */
      mpn_toom8h_mul (p, a, n, b, n);
    }
#endif
}

void
mpn_sqr (mp_ptr p, mp_srcptr a, mp_size_t n)
{
  ASSERT (n >= 1);
  ASSERT (! MPN_OVERLAP_P (p, 2 * n, a, n));

#if 0
  /* FIXME: Can this be removed? */
  if (n == 0)
    return;
#endif

  if (BELOW_THRESHOLD (n, SQR_BASECASE_THRESHOLD))
    { 
      /* mul_basecase is faster than sqr_basecase on small sizes sometimes */
      mpn_mul_basecase (p, a, n, a, n);
    }
  else if (BELOW_THRESHOLD (n, SQR_KARATSUBA_THRESHOLD))
    {
      mpn_sqr_basecase (p, a, n);
    }
  else if (BELOW_THRESHOLD (n, SQR_TOOM3_THRESHOLD))
    {
      /* Allocate workspace of fixed size on stack: fast! */
      mp_limb_t ws[MPN_KARA_SQR_N_TSIZE (SQR_TOOM3_THRESHOLD_LIMIT-1)];
      ASSERT (SQR_TOOM3_THRESHOLD <= SQR_TOOM3_THRESHOLD_LIMIT);
      mpn_kara_sqr_n (p, a, n, ws);
    }
  else if (BELOW_THRESHOLD (n, SQR_TOOM4_THRESHOLD))
    {
      mp_ptr ws;
      TMP_SDECL;
      TMP_SMARK;
      ws = TMP_SALLOC_LIMBS (MPN_TOOM3_SQR_N_TSIZE (n));
      mpn_toom3_sqr_n (p, a, n, ws);
      TMP_SFREE;
    }
  else if (BELOW_THRESHOLD (n, SQR_TOOM8_THRESHOLD))
    {
       mpn_toom4_sqr_n (p, a, n);
    }
#if WANT_FFT || TUNE_PROGRAM_BUILD
  else if (BELOW_THRESHOLD (n, SQR_FFT_FULL_THRESHOLD))
#else
  else 
#endif
    {
       mpn_toom8_sqr_n (p, a, n);
    }
#if WANT_FFT || TUNE_PROGRAM_BUILD
  else
    {
       mpn_mul_fft_main(p, a, n, a, n); 
    }
#endif
}

#if GMP_NAIL_BITS == 0

mp_limb_t
mpn_addmul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
	mp_limb_t ul, cl, hpl, lpl, rl;

	ASSERT(n >= 1);
	ASSERT(MPN_SAME_OR_SEPARATE_P(rp, up, n));

	cl = 0;
	do
	{
		ul = *up++;
		umul_ppmm(hpl, lpl, ul, vl);

		lpl += cl;
		cl = (lpl < cl) + hpl;

		rl = *rp;
		lpl = rl + lpl;
		cl += lpl < rl;
		*rp++ = lpl;
	} while (--n != 0);

	return cl;
}

#endif

#if GMP_NAIL_BITS == 1

mp_limb_t
mpn_addmul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
	mp_limb_t shifted_vl, ul, rl, lpl, hpl, prev_hpl, cl, xl, c1, c2, c3;

	ASSERT(n >= 1);
	ASSERT(MPN_SAME_OR_SEPARATE_P(rp, up, n));
	ASSERT_MPN(rp, n);
	ASSERT_MPN(up, n);
	ASSERT_LIMB(vl);

	shifted_vl = vl << GMP_NAIL_BITS;
	cl = 0;
	prev_hpl = 0;
	do
	{
		ul = *up++;
		rl = *rp;
		umul_ppmm(hpl, lpl, ul, shifted_vl);
		lpl >>= GMP_NAIL_BITS;
		ADDC_LIMB(c1, xl, prev_hpl, lpl);
		ADDC_LIMB(c2, xl, xl, rl);
		ADDC_LIMB(c3, xl, xl, cl);
		cl = c1 + c2 + c3;
		*rp++ = xl;
		prev_hpl = hpl;
	} while (--n != 0);

	return prev_hpl + cl;
}

#endif

#if GMP_NAIL_BITS >= 2

mp_limb_t
mpn_addmul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
	mp_limb_t shifted_vl, ul, rl, lpl, hpl, prev_hpl, xw, cl, xl;

	ASSERT(n >= 1);
	ASSERT(MPN_SAME_OR_SEPARATE_P(rp, up, n));
	ASSERT_MPN(rp, n);
	ASSERT_MPN(up, n);
	ASSERT_LIMB(vl);

	shifted_vl = vl << GMP_NAIL_BITS;
	cl = 0;
	prev_hpl = 0;
	do
	{
		ul = *up++;
		rl = *rp;
		umul_ppmm(hpl, lpl, ul, shifted_vl);
		lpl >>= GMP_NAIL_BITS;
		xw = prev_hpl + lpl + rl + cl;
		cl = xw >> GMP_NUMB_BITS;
		xl = xw & GMP_NUMB_MASK;
		*rp++ = xl;
		prev_hpl = hpl;
	} while (--n != 0);

	return prev_hpl + cl;
}

#endif
mp_limb_t
mpn_addmul_2(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_srcptr vp)
{
	rp[n] = mpn_addmul_1(rp, up, n, vp[0]);
	return mpn_addmul_1(rp + 1, up, n, vp[1]);
}

#if GMP_NAIL_BITS == 0

mp_limb_t
mpn_submul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
	mp_limb_t ul, cl, hpl, lpl, rl;

	ASSERT(n >= 1);
	ASSERT(MPN_SAME_OR_SEPARATE_P(rp, up, n));

	cl = 0;
	do
	{
		ul = *up++;
		umul_ppmm(hpl, lpl, ul, vl);

		lpl += cl;
		cl = (lpl < cl) + hpl;

		rl = *rp;
		lpl = rl - lpl;
		cl += lpl > rl;
		*rp++ = lpl;
	} while (--n != 0);

	return cl;
}

#endif

#if GMP_NAIL_BITS == 1

mp_limb_t
mpn_submul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
	mp_limb_t shifted_vl, ul, rl, lpl, hpl, prev_hpl, cl, xl, c1, c2, c3;

	ASSERT(n >= 1);
	ASSERT(MPN_SAME_OR_SEPARATE_P(rp, up, n));
	ASSERT_MPN(rp, n);
	ASSERT_MPN(up, n);
	ASSERT_LIMB(vl);

	shifted_vl = vl << GMP_NAIL_BITS;
	cl = 0;
	prev_hpl = 0;
	do
	{
		ul = *up++;
		rl = *rp;
		umul_ppmm(hpl, lpl, ul, shifted_vl);
		lpl >>= GMP_NAIL_BITS;
		SUBC_LIMB(c1, xl, rl, prev_hpl);
		SUBC_LIMB(c2, xl, xl, lpl);
		SUBC_LIMB(c3, xl, xl, cl);
		cl = c1 + c2 + c3;
		*rp++ = xl;
		prev_hpl = hpl;
	} while (--n != 0);

	return prev_hpl + cl;
}

#endif

#if GMP_NAIL_BITS >= 2

mp_limb_t
mpn_submul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
	mp_limb_t shifted_vl, ul, rl, lpl, hpl, prev_hpl, xw, cl, xl;

	ASSERT(n >= 1);
	ASSERT(MPN_SAME_OR_SEPARATE_P(rp, up, n));
	ASSERT_MPN(rp, n);
	ASSERT_MPN(up, n);
	ASSERT_LIMB(vl);

	shifted_vl = vl << GMP_NAIL_BITS;
	cl = 0;
	prev_hpl = 0;
	do
	{
		ul = *up++;
		rl = *rp;
		umul_ppmm(hpl, lpl, ul, shifted_vl);
		lpl >>= GMP_NAIL_BITS;
		xw = rl - (prev_hpl + lpl) + cl;
		cl = (mp_limb_signed_t)xw >> GMP_NUMB_BITS; /* FIXME: non-portable */
		xl = xw & GMP_NUMB_MASK;
		*rp++ = xl;
		prev_hpl = hpl;
	} while (--n != 0);

	return prev_hpl - cl;
}

#endif
mp_limb_t
mpn_rshift(mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt)
{
	mp_limb_t high_limb, low_limb;
	unsigned int tnc;
	mp_size_t i;
	mp_limb_t retval;

	ASSERT(n >= 1);
	ASSERT(cnt >= 1);
	ASSERT(cnt < GMP_NUMB_BITS);
	ASSERT(MPN_SAME_OR_INCR_P(rp, up, n));

	tnc = GMP_NUMB_BITS - cnt;
	high_limb = *up++;
	retval = (high_limb << tnc) & GMP_NUMB_MASK;
	low_limb = high_limb >> cnt;

	for (i = n - 1; i != 0; i--)
	{
		high_limb = *up++;
		*rp++ = low_limb | ((high_limb << tnc) & GMP_NUMB_MASK);
		low_limb = high_limb >> cnt;
	}
	*rp = low_limb;

	return retval;
}

mp_limb_t
mpn_lshift(mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt)
{
	mp_limb_t high_limb, low_limb;
	unsigned int tnc;
	mp_size_t i;
	mp_limb_t retval;

	ASSERT(n >= 1);
	ASSERT(cnt >= 1);
	ASSERT(cnt < GMP_NUMB_BITS);
	ASSERT(MPN_SAME_OR_DECR_P(rp, up, n));

	up += n;
	rp += n;

	tnc = GMP_NUMB_BITS - cnt;
	low_limb = *--up;
	retval = low_limb >> tnc;
	high_limb = (low_limb << cnt) & GMP_NUMB_MASK;

	for (i = n - 1; i != 0; i--)
	{
		low_limb = *--up;
		*--rp = high_limb | (low_limb >> tnc);
		high_limb = (low_limb << cnt) & GMP_NUMB_MASK;
	}
	*--rp = high_limb;

	return retval;
}

void
mpn_mul_basecase(mp_ptr rp,
mp_srcptr up, mp_size_t un,
mp_srcptr vp, mp_size_t vn)
{
	ASSERT(un >= vn);
	ASSERT(vn >= 1);
	ASSERT(!MPN_OVERLAP_P(rp, un + vn, up, un));
	ASSERT(!MPN_OVERLAP_P(rp, un + vn, vp, vn));

	/* We first multiply by the low order limb (or depending on optional function
	availability, limbs).  This result can be stored, not added, to rp.  We
	also avoid a loop for zeroing this way.  */

#if HAVE_NATIVE_mpn_mul_2
	if (vn >= 2)
	{
		rp[un + 1] = mpn_mul_2(rp, up, un, vp);
		rp += 2, vp += 2, vn -= 2;
	}
	else
	{
		rp[un] = mpn_mul_1(rp, up, un, vp[0]);
		return;
	}
#else
	rp[un] = mpn_mul_1(rp, up, un, vp[0]);
	rp += 1, vp += 1, vn -= 1;
#endif

	/* Now accumulate the product of up[] and the next higher limb (or depending
	on optional function availability, limbs) from vp[].  */

#define MAX_LEFT MP_SIZE_T_MAX	/* Used to simplify loops into if statements */


#if HAVE_NATIVE_mpn_addmul_6
	while (vn >= 6)
	{
		rp[un + 6 - 1] = mpn_addmul_6(rp, up, un, vp);
		if (MAX_LEFT == 6)
			return;
		rp += 6, vp += 6, vn -= 6;
		if (MAX_LEFT < 2 * 6)
			break;
	}
#undef MAX_LEFT
#define MAX_LEFT (6 - 1)
#endif

#if HAVE_NATIVE_mpn_addmul_5
	while (vn >= 5)
	{
		rp[un + 5 - 1] = mpn_addmul_5(rp, up, un, vp);
		if (MAX_LEFT == 5)
			return;
		rp += 5, vp += 5, vn -= 5;
		if (MAX_LEFT < 2 * 5)
			break;
	}
#undef MAX_LEFT
#define MAX_LEFT (5 - 1)
#endif

#if HAVE_NATIVE_mpn_addmul_4
	while (vn >= 4)
	{
		rp[un + 4 - 1] = mpn_addmul_4(rp, up, un, vp);
		if (MAX_LEFT == 4)
			return;
		rp += 4, vp += 4, vn -= 4;
		if (MAX_LEFT < 2 * 4)
			break;
	}
#undef MAX_LEFT
#define MAX_LEFT (4 - 1)
#endif

#if HAVE_NATIVE_mpn_addmul_3
	while (vn >= 3)
	{
		rp[un + 3 - 1] = mpn_addmul_3(rp, up, un, vp);
		if (MAX_LEFT == 3)
			return;
		rp += 3, vp += 3, vn -= 3;
		if (MAX_LEFT < 2 * 3)
			break;
	}
#undef MAX_LEFT
#define MAX_LEFT (3 - 1)
#endif

#if HAVE_NATIVE_mpn_addmul_2
	while (vn >= 2)
	{
		rp[un + 2 - 1] = mpn_addmul_2(rp, up, un, vp);
		if (MAX_LEFT == 2)
			return;
		rp += 2, vp += 2, vn -= 2;
		if (MAX_LEFT < 2 * 2)
			break;
	}
#undef MAX_LEFT
#define MAX_LEFT (2 - 1)
#endif

	while (vn >= 1)
	{
		rp[un] = mpn_addmul_1(rp, up, un, vp[0]);
		if (MAX_LEFT == 1)
			return;
		rp += 1, vp += 1, vn -= 1;
	}
}
