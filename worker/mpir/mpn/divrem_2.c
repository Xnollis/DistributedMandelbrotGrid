/* mpn_divrem_2 -- Divide natural numbers, producing both remainder and
   quotient.  The divisor is two limbs.

   THIS FILE CONTAINS INTERNAL FUNCTIONS WITH MUTABLE INTERFACES.  IT IS
   ONLY SAFE TO REACH THEM THROUGH DOCUMENTED INTERFACES.  IN FACT, IT IS
   ALMOST GUARANTEED THAT THEY'LL CHANGE OR DISAPPEAR IN A FUTURE GNU MP
   RELEASE.


Copyright 1993, 1994, 1995, 1996, 1999, 2000, 2001, 2002 Free Software
Foundation, Inc.

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

#ifndef DIVREM_2_THRESHOLD
#define DIVREM_2_THRESHOLD  0
#endif


/* Divide num (NP/NSIZE) by den (DP/2) and write
   the NSIZE-2 least significant quotient limbs at QP
   and the 2 long remainder at NP.  If QEXTRA_LIMBS is
   non-zero, generate that many fraction bits and append them after the
   other quotient limbs.
   Return the most significant limb of the quotient, this is always 0 or 1.

   Preconditions:
   0. NSIZE >= 2.
   1. The most significant bit of the divisor must be set.
   2. QP must either not overlap with the input operands at all, or
      QP + 2 >= NP must hold true.  (This means that it's
      possible to put the quotient in the high part of NUM, right after the
      remainder in NUM.
   3. NSIZE >= 2, even if QEXTRA_LIMBS is non-zero.  */

__GMP_DECLSPEC
mp_limb_t
mpn_divrem_2 (mp_ptr qp, mp_size_t qxn,
	      mp_ptr np, mp_size_t nn,
	      mp_srcptr dp)
{
  mp_limb_t most_significant_q_limb = 0;
  mp_size_t i;
  mp_limb_t n1, n0, n2;
  mp_limb_t d1, d0;
  mp_limb_t d1inv;
  int use_preinv;

  ASSERT (nn >= 2);
  ASSERT (qxn >= 0);
  ASSERT (dp[1] & GMP_NUMB_HIGHBIT);
  ASSERT (! MPN_OVERLAP_P (qp, nn-2+qxn, np, nn) || qp+2 >= np);
  ASSERT_MPN (np, nn);
  ASSERT_MPN (dp, 2);

#if HAVE_NATIVE_mpn_divrem_euclidean_qr_2
if (qxn==0) return mpn_divrem_euclidean_qr_2(qp,np,nn,dp);
#endif

  np += nn - 2;
  d1 = dp[1];
  d0 = dp[0];
  n1 = np[1];
  n0 = np[0];

  if (n1 >= d1 && (n1 > d1 || n0 >= d0))
    {
#if GMP_NAIL_BITS == 0
      sub_ddmmss (n1, n0, n1, n0, d1, d0);
#else
      n0 = n0 - d0;
      n1 = n1 - d1 - (n0 >> GMP_LIMB_BITS - 1);
      n0 &= GMP_NUMB_MASK;
#endif
      most_significant_q_limb = 1;
    }

  use_preinv = ABOVE_THRESHOLD (qxn + nn - 2, DIVREM_2_THRESHOLD);
  if (use_preinv)
    invert_limb (d1inv, d1);

  for (i = qxn + nn - 2 - 1; i >= 0; i--)
    {
      mp_limb_t q;
      mp_limb_t r;

      if (i >= qxn)
	np--;
      else
	np[0] = 0;

      if (n1 == d1)
	{
	  /* Q should be either 111..111 or 111..110.  Need special handling
	     of this rare case as normal division would give overflow.  */
	  q = GMP_NUMB_MASK;

	  r = (n0 + d1) & GMP_NUMB_MASK;
	  if (r < d1)	/* Carry in the addition? */
	    {
#if GMP_NAIL_BITS == 0
	      add_ssaaaa (n1, n0, r - d0, np[0], 0, d0);
#else
	      n0 = np[0] + d0;
	      n1 = (r - d0 + (n0 >> GMP_NUMB_BITS)) & GMP_NUMB_MASK;
	      n0 &= GMP_NUMB_MASK;
#endif
	      qp[i] = q;
	      continue;
	    }
	  n1 = d0 - (d0 != 0);
	  n0 = -d0 & GMP_NUMB_MASK;
	}
      else
	{
	  if (use_preinv)
	    udiv_qrnnd_preinv (q, r, n1, n0, d1, d1inv);
	  else
	    udiv_qrnnd (q, r, n1, n0 << GMP_NAIL_BITS, d1 << GMP_NAIL_BITS);
	  r >>= GMP_NAIL_BITS;
	  umul_ppmm (n1, n0, d0, q << GMP_NAIL_BITS);
	  n0 >>= GMP_NAIL_BITS;
	}

      n2 = np[0];

    q_test:
      if (n1 > r || (n1 == r && n0 > n2))
	{
	  /* The estimated Q was too large.  */
	  q--;

#if GMP_NAIL_BITS == 0
	  sub_ddmmss (n1, n0, n1, n0, 0, d0);
#else
	  n0 = n0 - d0;
	  n1 = n1 - (n0 >> GMP_LIMB_BITS - 1);
	  n0 &= GMP_NUMB_MASK;
#endif
	  r += d1;
	  if (r >= d1)	/* If not carry, test Q again.  */
	    goto q_test;
	}

      qp[i] = q;
#if GMP_NAIL_BITS == 0
      sub_ddmmss (n1, n0, r, n2, n1, n0);
#else
      n0 = n2 - n0;
      n1 = r - n1 - (n0 >> GMP_LIMB_BITS - 1);
      n0 &= GMP_NUMB_MASK;
#endif
    }
  np[1] = n1;
  np[0] = n0;

  return most_significant_q_limb;
}

__GMP_DECLSPEC
mp_limb_t
mpn_sb_div_qr(mp_ptr qp,
mp_ptr np, mp_size_t nn,
mp_srcptr dp, mp_size_t dn,
mp_limb_t dinv)
{
	mp_limb_t qh;
	mp_size_t i;
	mp_limb_t n1, n0;
	mp_limb_t d1, d0;
	mp_limb_t cy, cy2;
	mp_limb_t q;

	ASSERT(dn > 2);
	ASSERT(nn >= dn);
	ASSERT((dp[dn - 1] & GMP_NUMB_HIGHBIT) != 0);

	np += nn;

	qh = mpn_cmp(np - dn, dp, dn) >= 0;
	if (qh != 0)
		mpn_sub_n(np - dn, np - dn, dp, dn);

	d1 = dp[dn - 1];

	qp += nn - dn;

	dn -= 2;			/* offset dn by 2 for main division loops,
						saving two iterations in mpn_submul_1.  */
	d0 = dp[dn];

	np -= 2;

	n1 = np[1];

	for (i = nn - (dn + 2); i > 0; i--)
	{
		np--;
		if (UNLIKELY(n1 == d1) && np[1] == d0)
		{
			q = GMP_NUMB_MASK;
			mpn_submul_1(np - dn, dp, dn + 2, q);
			n1 = np[1];		/* update n1, last loop's value will now be invalid */
		}
		else
		{
			udiv_qr_3by2(q, n1, n0, n1, np[1], np[0], d1, d0, dinv);

			cy2 = mpn_submul_1(np - dn, dp, dn, q);

			sub_333(cy, n1, n0, 0, n1, n0, 0, 0, cy2);

			np[0] = n0;

			if (UNLIKELY(cy != 0))
			{
				n1 += d1 + mpn_add_n(np - dn, np - dn, dp, dn + 1);
				q--;
			}
		}

		*--qp = q;
	}
	np[1] = n1;

	return qh;
}
__GMP_DECLSPEC
mp_limb_t
mpn_dc_div_qr(mp_ptr qp,
mp_ptr np, mp_size_t nn,
mp_srcptr dp, mp_size_t dn,
mp_limb_t dinv)
{
	mp_size_t qn;
	mp_limb_t qh, cy;
	mp_ptr tp;
	TMP_DECL;

	TMP_MARK;

	ASSERT(dn >= 6);		/* to adhere to mpn_sb_div_qr's limits */
	ASSERT(nn - dn >= 3);	/* to adhere to mpn_sb_div_qr's limits */
	ASSERT(dp[dn - 1] & GMP_NUMB_HIGHBIT);

	tp = TMP_ALLOC_LIMBS(DC_DIVAPPR_Q_N_ITCH(dn));

	qn = nn - dn;
	qp += qn;
	np += nn;
	dp += dn;

	if (qn > dn)
	{
		/* Reduce qn mod dn without division, optimizing small operations.  */
		do
			qn -= dn;
		while (qn > dn);

		qp -= qn;			/* point at low limb of next quotient block */
		np -= qn;			/* point in the middle of partial remainder */

		/* Perform the typically smaller block first.  */
		if (qn == 1)
		{
			mp_limb_t q, n2, n1, n0, d1, d0;

			/* Handle qh up front, for simplicity. */
			qh = mpn_cmp(np - dn + 1, dp - dn, dn) >= 0;
			if (qh)
				ASSERT_NOCARRY(mpn_sub_n(np - dn + 1, np - dn + 1, dp - dn, dn));

			/* A single iteration of schoolbook: One 3/2 division,
			followed by the bignum update and adjustment. */
			n2 = np[0];
			n1 = np[-1];
			n0 = np[-2];
			d1 = dp[-1];
			d0 = dp[-2];

			ASSERT(n2 < d1 || (n2 == d1 && n1 <= d0));

			if (UNLIKELY(n2 == d1) && n1 == d0)
			{
				q = GMP_NUMB_MASK;
				cy = mpn_submul_1(np - dn, dp - dn, dn, q);
				ASSERT(cy == n2);
			}
			else
			{
				udiv_qr_3by2(q, n1, n0, n2, n1, n0, d1, d0, dinv);

				if (dn > 2)
				{
					mp_limb_t cy, cy1;
					cy = mpn_submul_1(np - dn, dp - dn, dn - 2, q);

					cy1 = n0 < cy;
					n0 = (n0 - cy) & GMP_NUMB_MASK;
					cy = n1 < cy1;
					n1 = (n1 - cy1) & GMP_NUMB_MASK;
					np[-2] = n0;

					if (UNLIKELY(cy != 0))
					{
						n1 += d1 + mpn_add_n(np - dn, np - dn, dp - dn, dn - 1);
						qh -= (q == 0);
						q = (q - 1) & GMP_NUMB_MASK;
					}
				}
				else
					np[-2] = n0;

				np[-1] = n1;
			}
			qp[0] = q;
		}
		else
		{
			/* Do a 2qn / qn division */
			if (qn == 2)
				qh = mpn_divrem_2(qp, 0L, np - 2, 4, dp - 2); /* FIXME: obsolete function. Use 5/3 division? */
			else if (BELOW_THRESHOLD(qn, DC_DIV_QR_THRESHOLD))
				qh = mpn_sb_div_qr(qp, np - qn, 2 * qn, dp - qn, qn, dinv);
			else
				qh = mpn_dc_div_qr_n(qp, np - qn, dp - qn, qn, dinv, tp);

			if (qn != dn)
			{
				if (qn > dn - qn)
					mpn_mul(tp, qp, qn, dp - dn, dn - qn);
				else
					mpn_mul(tp, dp - dn, dn - qn, qp, qn);

				cy = mpn_sub_n(np - dn, np - dn, tp, dn);
				if (qh != 0)
					cy += mpn_sub_n(np - dn + qn, np - dn + qn, dp - dn, dn - qn);

				while (cy != 0)
				{
					qh -= mpn_sub_1(qp, qp, qn, 1);
					cy -= mpn_add_n(np - dn, np - dn, dp - dn, dn);
				}
			}
		}

		qn = nn - dn - qn;
		do
		{
			qp -= dn;
			np -= dn;
			ASSERT_NOCARRY(mpn_dc_div_qr_n(qp, np - dn, dp - dn, dn, dinv, tp));
			qn -= dn;
		} while (qn > 0);
	}
	else
	{
		qp -= qn;			/* point at low limb of next quotient block */
		np -= qn;			/* point in the middle of partial remainder */

		if (BELOW_THRESHOLD(qn, DC_DIV_QR_THRESHOLD))
			qh = mpn_sb_div_qr(qp, np - qn, 2 * qn, dp - qn, qn, dinv);
		else
			qh = mpn_dc_div_qr_n(qp, np - qn, dp - qn, qn, dinv, tp);

		if (qn != dn)
		{
			if (qn > dn - qn)
				mpn_mul(tp, qp, qn, dp - dn, dn - qn);
			else
				mpn_mul(tp, dp - dn, dn - qn, qp, qn);

			cy = mpn_sub_n(np - dn, np - dn, tp, dn);
			if (qh != 0)
				cy += mpn_sub_n(np - dn + qn, np - dn + qn, dp - dn, dn - qn);

			while (cy != 0)
			{
				qh -= mpn_sub_1(qp, qp, qn, 1);
				cy -= mpn_add_n(np - dn, np - dn, dp - dn, dn);
			}
		}
	}

	TMP_FREE;
	return qh;
}
__GMP_DECLSPEC
mp_limb_t
mpn_dc_div_qr_n(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n,
mp_limb_t dinv, mp_ptr tp)
{
	mp_size_t lo, hi;
	mp_limb_t cy, qh, ql;

	lo = n >> 1; /* floor(n/2) */
	hi = n - lo;	/* ceil(n/2) */


	if (BELOW_THRESHOLD(hi, DC_DIV_QR_THRESHOLD))
		qh = mpn_sb_div_qr(qp + lo, np + 2 * lo, 2 * hi, dp + lo, hi, dinv);
	else
		qh = mpn_dc_div_qr_n(qp + lo, np + 2 * lo, dp + lo, hi, dinv, tp);

	mpn_mul(tp, qp + lo, hi, dp, lo);

	cy = mpn_sub_n(np + lo, np + lo, tp, n);
	if (qh != 0)
		cy += mpn_sub_n(np + n, np + n, dp, lo);

	while (cy != 0)
	{
		qh -= mpn_sub_1(qp + lo, qp + lo, hi, 1);
		cy -= mpn_add_n(np + lo, np + lo, dp, n);
	}

	if (BELOW_THRESHOLD(lo, DC_DIV_QR_THRESHOLD))
		ql = mpn_sb_div_qr(qp, np + hi, 2 * lo, dp + hi, lo, dinv);
	else
		ql = mpn_dc_div_qr_n(qp, np + hi, dp + hi, lo, dinv, tp);

	mpn_mul(tp, dp, hi, qp, lo);

	cy = mpn_sub_n(np, np, tp, n);
	if (ql != 0)
		cy += mpn_sub_n(np + lo, np + lo, dp, hi);

	while (cy != 0)
	{
		mpn_sub_1(qp, qp, lo, 1);
		cy -= mpn_add_n(np, np, dp, n);
	}

	return qh;
}
__GMP_DECLSPEC
mp_limb_t
mpn_inv_div_qr(mp_ptr qp,
mp_ptr np, mp_size_t nn,
mp_srcptr dp, mp_size_t dn,
mp_srcptr dinv)
{
	mp_size_t qn;
	mp_limb_t qh, cy, dinv2;
	mp_ptr tp;
	TMP_DECL;

	TMP_MARK;

	ASSERT(dn >= 6);		/* to adhere to mpn_sb_div_qr's limits */
	ASSERT(nn - dn >= 3);	/* to adhere to mpn_sb_div_qr's limits */
	ASSERT(dp[dn - 1] & GMP_NUMB_HIGHBIT);

	mpir_invert_pi1(dinv2, dp[dn - 1], dp[dn - 2]);

	tp = TMP_ALLOC_LIMBS(DC_DIVAPPR_Q_N_ITCH(dn));

	qn = nn - dn;
	qp += qn;
	np += nn;
	dp += dn;

	if (qn > dn)
	{
		/* Reduce qn mod dn without division, optimizing small operations.  */
		do
			qn -= dn;
		while (qn > dn);

		qp -= qn;			/* point at low limb of next quotient block */
		np -= qn;			/* point in the middle of partial remainder */

		/* Perform the typically smaller block first.  */
		if (qn == 1)
		{
			mp_limb_t q, n2, n1, n0, d1, d0, d11, d01;

			/* Handle qh up front, for simplicity. */
			qh = mpn_cmp(np - dn + 1, dp - dn, dn) >= 0;
			if (qh)
				ASSERT_NOCARRY(mpn_sub_n(np - dn + 1, np - dn + 1, dp - dn, dn));

			/* A single iteration of schoolbook: One 3/2 division,
			followed by the bignum update and adjustment. */
			n2 = np[0];
			n1 = np[-1];
			n0 = np[-2];
			d1 = dp[-1];
			d0 = dp[-2];
			d01 = d0 + 1;
			d11 = d1 + (d01 < d0);

			ASSERT(n2 < d1 || (n2 == d1 && n1 <= d0));

			if (UNLIKELY(n2 == d1) && n1 == d0)
			{
				q = GMP_NUMB_MASK;
				cy = mpn_submul_1(np - dn, dp - dn, dn, q);
				ASSERT(cy == n2);
			}
			else
			{
				udiv_qr_3by2(q, n1, n0, n2, n1, n0, d1, d0, dinv2);

				if (dn > 2)
				{
					mp_limb_t cy, cy1;
					cy = mpn_submul_1(np - dn, dp - dn, dn - 2, q);

					cy1 = n0 < cy;
					n0 = (n0 - cy) & GMP_NUMB_MASK;
					cy = n1 < cy1;
					n1 = (n1 - cy1) & GMP_NUMB_MASK;
					np[-2] = n0;

					if (UNLIKELY(cy != 0))
					{
						n1 += d1 + mpn_add_n(np - dn, np - dn, dp - dn, dn - 1);
						qh -= (q == 0);
						q = (q - 1) & GMP_NUMB_MASK;
					}
				}
				else
					np[-2] = n0;

				np[-1] = n1;
			}
			qp[0] = q;
			REMOVE_WARNINGS_OF_LOCAL_VAR(d11);
		}
		else
		{
			/* Do a 2qn / qn division */
			if (qn == 2)
				qh = mpn_divrem_2(qp, 0L, np - 2, 4, dp - 2); /* FIXME: obsolete function. Use 5/3 division? */
			else if (BELOW_THRESHOLD(qn, DC_DIV_QR_THRESHOLD))
				qh = mpn_sb_div_qr(qp, np - qn, 2 * qn, dp - qn, qn, dinv2);
			else if (BELOW_THRESHOLD(qn, INV_DIV_QR_THRESHOLD))
				qh = mpn_dc_div_qr_n(qp, np - qn, dp - qn, qn, dinv2, tp);
			else
			{
				mpn_invert_trunc(tp, qn, dinv, dn, dp - dn);
				qh = mpn_inv_div_qr_n(qp, np - qn, dp - qn, qn, tp);
			}

			if (qn != dn)
			{
				if (qn > dn - qn)
					mpn_mul(tp, qp, qn, dp - dn, dn - qn);
				else
					mpn_mul(tp, dp - dn, dn - qn, qp, qn);

				cy = mpn_sub_n(np - dn, np - dn, tp, dn);
				if (qh != 0)
					cy += mpn_sub_n(np - dn + qn, np - dn + qn, dp - dn, dn - qn);

				while (cy != 0)
				{
					qh -= mpn_sub_1(qp, qp, qn, 1);
					cy -= mpn_add_n(np - dn, np - dn, dp - dn, dn);
				}
			}
		}

		qn = nn - dn - qn;
		do
		{
			qp -= dn;
			np -= dn;
			ASSERT_NOCARRY(mpn_inv_div_qr_n(qp, np - dn, dp - dn, dn, dinv));
			qn -= dn;
		} while (qn > 0);
	}
	else
	{
		qp -= qn;			/* point at low limb of next quotient block */
		np -= qn;			/* point in the middle of partial remainder */

		if (BELOW_THRESHOLD(qn, DC_DIV_QR_THRESHOLD))
			qh = mpn_sb_div_qr(qp, np - qn, 2 * qn, dp - qn, qn, dinv2);
		else if (BELOW_THRESHOLD(qn, INV_DIV_QR_THRESHOLD))
			qh = mpn_dc_div_qr_n(qp, np - qn, dp - qn, qn, dinv2, tp);
		else
		{
			mpn_invert_trunc(tp, qn, dinv, dn, dp - dn);
			qh = mpn_inv_div_qr_n(qp, np - qn, dp - qn, qn, tp);
		}

		if (qn != dn)
		{
			if (qn > dn - qn)
				mpn_mul(tp, qp, qn, dp - dn, dn - qn);
			else
				mpn_mul(tp, dp - dn, dn - qn, qp, qn);

			cy = mpn_sub_n(np - dn, np - dn, tp, dn);
			if (qh != 0)
				cy += mpn_sub_n(np - dn + qn, np - dn + qn, dp - dn, dn - qn);

			while (cy != 0)
			{
				qh -= mpn_sub_1(qp, qp, qn, 1);
				cy -= mpn_add_n(np - dn, np - dn, dp - dn, dn);
			}
		}
	}

	TMP_FREE;
	return qh;
}
__GMP_DECLSPEC
mp_limb_t
mpn_inv_div_qr_n(mp_ptr qp, mp_ptr np,
mp_srcptr dp, mp_size_t dn, mp_srcptr inv)
{
	mp_limb_t cy, lo, ret = 0, ret2 = 0;
	mp_size_t m;
	mp_ptr tp;
	TMP_DECL;

	TMP_MARK;

	ASSERT(mpn_is_invert(inv, dp, dn));

	if (mpn_cmp(np + dn, dp, dn) >= 0)
	{
		ret2 = 1;
		mpn_sub_n(np + dn, np + dn, dp, dn);
	}

	tp = TMP_ALLOC_LIMBS(2 * dn + 1);
	mpn_mul(tp, np + dn - 1, dn + 1, inv, dn);
	add_ssaaaa(cy, lo, 0, np[dn - 1], 0, tp[dn]);
	ret += mpn_add_n(qp, tp + dn + 1, np + dn, dn);
	ret += mpn_add_1(qp, qp, dn, cy);

	/*
	Let X = B^dn + inv, D = { dp, dn }, N = { np, 2*dn }, then
	DX < B^{2*dn} <= D(X+1), thus
	Let N' = { np + n - 1, n + 1 }
	N'X/B^{dn+1} < B^{dn-1}N'/D <= N'X/B^{dn+1} + N'/B^{dn+1} < N'X/B^{dn+1} + 1
	N'X/B^{dn+1} < N/D <= N'X/B^{dn+1} + 1 + 2/B
	There is either one integer in this range, or two. However, in the latter case
	the left hand bound is either an integer or < 2/B below one.
	*/

	if (UNLIKELY(ret == 1))
	{
		ret -= mpn_sub_1(qp, qp, dn, 1);
		ASSERT(ret == 0);
	}

	ret -= mpn_sub_1(qp, qp, dn, 1);
	if (UNLIKELY(ret == ~CNST_LIMB(0)))
		ret += mpn_add_1(qp, qp, dn, 1);
	/* ret is now guaranteed to be 0 or 1*/
	ASSERT(ret == 0);

	m = dn + 1;
	if ((dn <= MPN_FFT_MUL_N_MINSIZE) || (ret))
	{
		mpn_mul_n(tp, qp, dp, dn);
	}
#if WANT_FFT
	else
	{
		mp_limb_t cy, cy2;

		if (m >= FFT_MULMOD_2EXPP1_CUTOFF)
			m = mpir_fft_adjust_limbs(m);
		cy = mpn_mulmod_Bexpp1_fft(tp, m, qp, dn, dp, dn);

		/* cy, {tp, m} = qp * dp mod (B^m+1) */
		cy2 = mpn_add_n(tp, tp, np + m, 2 * dn - m);
		mpn_add_1(tp + 2 * dn - m, tp + 2 * dn - m, 2 * m - 2 * dn, cy2);

		/* Make correction */
		mpn_sub_1(tp, tp, m, tp[0] - dp[0] * qp[0]);
	}
#endif
	mpn_sub_n(np, np, tp, m);
	MPN_ZERO(np + m, 2 * dn - m);
	while (np[dn] || mpn_cmp(np, dp, dn) >= 0)
	{
		ret += mpn_add_1(qp, qp, dn, 1);
		np[dn] -= mpn_sub_n(np, np, dp, dn);
	}

	/* Not possible for ret == 2 as we have qp*dp <= np */
	ASSERT(ret + ret2 < 2);

	TMP_FREE;

    REMOVE_WARNINGS_OF_LOCAL_VAR(lo);
	return ret + ret2;
}
