/* Implementation of the multiplication algorithm for Toom-Cook 8.5-way.

   Contributed to the GNU project by Marco Bodrato.

   THE FUNCTION IN THIS FILE IS INTERNAL WITH A MUTABLE INTERFACE.  IT IS ONLY
   SAFE TO REACH IT THROUGH DOCUMENTED INTERFACES.  IN FACT, IT IS ALMOST
   GUARANTEED THAT IT WILL CHANGE OR DISAPPEAR IN A FUTURE GNU MP RELEASE.

Copyright 2009, 2010 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library.  If not, see http://www.gnu.org/licenses/.  */


#include "mpir_inter_decl.h"

#if HAVE_NATIVE_mpn_addlsh_n
#define DO_mpn_addlsh_n(dst,src,n,s,ws) mpn_addlsh_n(dst,dst,src,n,s)
#else
static
__GMP_DECLSPEC
mp_limb_t
DO_mpn_addlsh_n(mp_ptr dst, mp_srcptr src, mp_size_t n, unsigned int s, mp_ptr ws)
{
#if USE_MUL_1 && 0
	return mpn_addmul_1(dst, src, n, CNST_LIMB(1) << (s));
#else
	mp_limb_t __cy;
	__cy = mpn_lshift(ws, src, n, s);
	return    __cy + mpn_add_n(dst, dst, ws, n);
#endif
}
#endif

/* Evaluates a polynomial of degree k >= 3. */
__GMP_DECLSPEC
int
mpn_toom_eval_pm2rexp(mp_ptr rp, mp_ptr rm,
unsigned int q, mp_srcptr ap, mp_size_t n, mp_size_t t,
unsigned int s, mp_ptr ws)
{
	unsigned int i;
	int neg;
	/* {ap,q*n+t} -> {rp,n+1} {rm,n+1} , with {ws, n+1}*/
	ASSERT(n >= t);
	ASSERT(s != 0); /* or _eval_pm1 should be used */
	ASSERT(q > 1);
	ASSERT(s*q < GMP_NUMB_BITS);
	rp[n] = mpn_lshift(rp, ap, n, s*q);
	ws[n] = mpn_lshift(ws, ap + n, n, s*(q - 1));
	if ((q & 1) != 0) {
		ASSERT_NOCARRY(mpn_add(ws, ws, n + 1, ap + n*q, t));
		rp[n] += DO_mpn_addlsh_n(rp, ap + n*(q - 1), n, s, rm);
	}
	else {
		ASSERT_NOCARRY(mpn_add(rp, rp, n + 1, ap + n*q, t));
	}
	for (i = 2; i < q - 1; i++)
	{
		rp[n] += DO_mpn_addlsh_n(rp, ap + n*i, n, s*(q - i), rm);
		i++;
		ws[n] += DO_mpn_addlsh_n(ws, ap + n*i, n, s*(q - i), rm);
	};

	neg = (mpn_cmp(rp, ws, n + 1) < 0) ? ~0 : 0;

#if HAVE_NATIVE_mpn_sumdiff_n
	if (neg)
		mpn_sumdiff_n(rp, rm, ws, rp, n + 1);
	else
		mpn_sumdiff_n(rp, rm, rp, ws, n + 1);
#else 
	if (neg)
		mpn_sub_n(rm, ws, rp, n + 1);
	else
		mpn_sub_n(rm, rp, ws, n + 1);

	ASSERT_NOCARRY(mpn_add_n(rp, rp, ws, n + 1));
#endif

	return neg;
}

/////////////////////////////////////////////////////////

/* Gets {pp,n} and (sign?-1:1)*{np,n}. Computes at once:
{pp,n} <- ({pp,n}+{np,n})/2^{ps+1}
{pn,n} <- ({pp,n}-{np,n})/2^{ns+1}
Finally recompose them obtaining:
{pp,n+off} <- {pp,n}+{np,n}*2^{off*GMP_NUMB_BITS}
*/
__GMP_DECLSPEC
void
mpn_toom_couple_handling(mp_ptr pp, mp_size_t n, mp_ptr np,
int nsign, mp_size_t off, int ps, int ns)
{
	if (nsign) {
#ifdef HAVE_NATIVE_mpn_rsh1sub_n
		mpn_rsh1sub_n(np, pp, np, n);
#else
		mpn_sub_n(np, pp, np, n);
		mpn_rshift(np, np, n, 1);
#endif
	}
	else {
#ifdef HAVE_NATIVE_mpn_rsh1add_n
		mpn_rsh1add_n(np, pp, np, n);
#else
		mpn_add_n(np, pp, np, n);
		mpn_rshift(np, np, n, 1);
#endif
	}

#ifdef HAVE_NATIVE_mpn_rsh1sub_n
	if (ps == 1)
		mpn_rsh1sub_n(pp, pp, np, n);
	else
#endif
	{
		mpn_sub_n(pp, pp, np, n);
		if (ps > 0)
			mpn_rshift(pp, pp, n, ps);
	}
	if (ns > 0)
		mpn_rshift(np, np, n, ns);
	pp[n] = mpn_add_n(pp + off, pp + off, np, n - off);
	ASSERT_NOCARRY(mpn_add_1(pp + n, np + n - off, off, pp[n]));
}

/////////////////////////////////////////////////////////

#if GMP_NUMB_BITS < 29
#error Not implemented.
#endif

#if GMP_NUMB_BITS < 43
#define BIT_CORRECTION 1
#define CORRECTION_BITS GMP_NUMB_BITS
#else
#define BIT_CORRECTION 0
#define CORRECTION_BITS 0
#endif

#define TOOM8H_MUL_N_REC(p, a, b, n)				\
  do {									\
    if (BELOW_THRESHOLD (n, MUL_TOOM8H_THRESHOLD))		\
      mpn_mul_n (p, a, b, n);				\
    else								\
      mpn_toom8h_mul (p, a, n, b, n);				\
  } while (0)

#define TOOM8H_MUL_REC(p, a, na, b, nb)		\
  do {	mpn_mul (p, a, na, b, nb);			\
  } while (0)

/* Toom-8.5 , compute the product {pp,an+bn} <- {ap,an} * {bp,bn}
   With: an >= bn >= 86, an*5 <  bn * 11.
   It _may_ work with bn<=?? and bn*?? < an*? < bn*??

   Evaluate in: infinity, +8,-8,+4,-4,+2,-2,+1,-1,+1/2,-1/2,+1/4,-1/4,+1/8,-1/8,0.
*/
/* Estimate on needed scratch:
   S(n) <= (n+7)\8*13+5+MAX(S((n+7)\8),1+2*(n+7)\8),
   since n>80; S(n) <= ceil(log(n/10)/log(8))*(13+5)+n*15\8 < n*15\8 + lg2(n)*6
 */

__GMP_DECLSPEC
void
mpn_toom8h_mul   (mp_ptr pp,
		  mp_srcptr ap, mp_size_t an,
		  mp_srcptr bp, mp_size_t bn)
{
  mp_size_t n, s, t;
  int p, q, half;
  int sign;
  mp_ptr scratch;

  TMP_DECL;

  TMP_MARK;

  /***************************** decomposition *******************************/

  ASSERT (an >= bn);
  /* Can not handle too small operands */
  ASSERT (bn >= 86);
  /* Can not handle too much unbalancement */
  ASSERT (an*4 <= bn*13);
  ASSERT (GMP_NUMB_BITS > 12*3 || an*4 <= bn*12);
  ASSERT (GMP_NUMB_BITS > 11*3 || an*5 <= bn*11);
  ASSERT (GMP_NUMB_BITS > 10*3 || an*6 <= bn*10);
  ASSERT (GMP_NUMB_BITS >  9*3 || an*7 <= bn* 9);

  /* Limit num/den is a rational number between
     (16/15)^(log(6)/log(2*6-1)) and (16/15)^(log(8)/log(2*8-1))             */
#define LIMIT_numerator (21)
#define LIMIT_denominat (20)

  if (LIKELY (an == bn) || an * (LIMIT_denominat>>1) < LIMIT_numerator * (bn>>1) ) /* is 8*... < 8*... */
    {
      half = 0;
      n = 1 + ((an - 1)>>3);
      p = q = 7;
      s = an - p * n;
      t = bn - q * n;
    }
  else
    {
      if (an * 13 < 16 * bn) /* (an*7*LIMIT_numerator<LIMIT_denominat*9*bn) */
	{ p = 9; q = 8; }
      else if (GMP_NUMB_BITS <= 9*3 ||
	       an *(LIMIT_denominat>>1) < (LIMIT_numerator/7*9) * (bn>>1))
	{ p = 9; q = 7; }
      else if (an * 10 < 33 * (bn>>1)) /* (an*3*LIMIT_numerator<LIMIT_denominat*5*bn) */
	{ p =10; q = 7; }
      else if (GMP_NUMB_BITS <= 10*3 ||
	       an * (LIMIT_denominat/5) < (LIMIT_numerator/3) * bn)
	{ p =10; q = 6; }
      else if (an * 6 < 13 * bn) /*(an * 5 * LIMIT_numerator < LIMIT_denominat *11 * bn)*/
	{ p =11; q = 6; }
      else if (GMP_NUMB_BITS <= 11*3 ||
	       an * 4 < 9 * bn)
	{ p =11; q = 5; }
      else if (an *(LIMIT_numerator/3) < LIMIT_denominat * bn )  /* is 4*... <12*... */
	{ p =12; q = 5; }
      else if (GMP_NUMB_BITS <= 12*3 ||
	       an * 9 < 28 * bn )  /* is 4*... <12*... */
	{ p =12; q = 4; }
      else
	{ p =13; q = 4; }

      half = (p+q)&1;
      n = 1 + (q * an >= p * bn ? (an - 1) / (size_t) p : (bn - 1) / (size_t) q);
      p--; q--;

      s = an - p * n;
      t = bn - q * n;

      if(half) { /* Recover from badly chosen splitting */
	if (s<1) {p--; s+=n; half=0;}
	else if (t<1) {q--; t+=n; half=0;}
      }
    }
#undef LIMIT_numerator
#undef LIMIT_denominat

  ASSERT (0 < s && s <= n);
  ASSERT (0 < t && t <= n);
  ASSERT (half || s + t > 3);
  ASSERT (n > 2);

  scratch = TMP_ALLOC_LIMBS(n*15 + GMP_LIMB_BITS*6);
  
#define   r6    (pp + 3 * n)			/* 3n+1 */
#define   r4    (pp + 7 * n)			/* 3n+1 */
#define   r2    (pp +11 * n)			/* 3n+1 */
#define   r0    (pp +15 * n)			/* s+t <= 2*n */
#define   r7    (scratch)			/* 3n+1 */
#define   r5    (scratch + 3 * n + 1)		/* 3n+1 */
#define   r3    (scratch + 6 * n + 2)		/* 3n+1 */
#define   r1    (scratch + 9 * n + 3)		/* 3n+1 */
#define   v0    (pp +11 * n)			/* n+1 */
#define   v1    (pp +12 * n+1)			/* n+1 */
#define   v2    (pp +13 * n+2)			/* n+1 */
#define   v3    (scratch +12 * n + 4)		/* n+1 */
#define   wsi   (scratch +12 * n + 4)		/* 3n+1 */

  /********************** evaluation and recursive calls *********************/

  /* $\pm1/8$ */
  sign = mpn_toom_eval_pm2rexp (v2, v0, p, ap, n, s, 3, pp) ^
	 mpn_toom_eval_pm2rexp (v3, v1, q, bp, n, t, 3, pp);
  TOOM8H_MUL_N_REC(pp, v0, v1, n + 1); /* A(-1/8)*B(-1/8)*8^. */
  TOOM8H_MUL_N_REC(r7, v2, v3, n + 1); /* A(+1/8)*B(+1/8)*8^. */
  mpn_toom_couple_handling (r7, 2 * n + 1 + BIT_CORRECTION, pp, sign, n, 3*(1+half), 3*(half));

  /* $\pm1/4$ */
  sign = mpn_toom_eval_pm2rexp (v2, v0, p, ap, n, s, 2, pp) ^
	 mpn_toom_eval_pm2rexp (v3, v1, q, bp, n, t, 2, pp);
  TOOM8H_MUL_N_REC(pp, v0, v1, n + 1); /* A(-1/4)*B(-1/4)*4^. */
  TOOM8H_MUL_N_REC(r5, v2, v3, n + 1); /* A(+1/4)*B(+1/4)*4^. */
  mpn_toom_couple_handling (r5, 2 * n + 1, pp, sign, n, 2*(1+half), 2*(half));

  /* $\pm2$ */
  sign = mpn_toom_eval_pm2 (v2, v0, p, ap, n, s, pp) ^
	 mpn_toom_eval_pm2 (v3, v1, q, bp, n, t, pp);
  TOOM8H_MUL_N_REC(pp, v0, v1, n + 1); /* A(-2)*B(-2) */
  TOOM8H_MUL_N_REC(r3, v2, v3, n + 1); /* A(+2)*B(+2) */
  mpn_toom_couple_handling (r3, 2 * n + 1, pp, sign, n, 1, 2);

  /* $\pm8$ */
  sign = mpn_toom_eval_pm2exp (v2, v0, p, ap, n, s, 3, pp) ^
	 mpn_toom_eval_pm2exp (v3, v1, q, bp, n, t, 3, pp);
  TOOM8H_MUL_N_REC(pp, v0, v1, n + 1); /* A(-8)*B(-8) */
  TOOM8H_MUL_N_REC(r1, v2, v3, n + 1); /* A(+8)*B(+8) */
  mpn_toom_couple_handling (r1, 2 * n + 1 + BIT_CORRECTION, pp, sign, n, 3, 6);

  /* $\pm1/2$ */
  sign = mpn_toom_eval_pm2rexp (v2, v0, p, ap, n, s, 1, pp) ^
	 mpn_toom_eval_pm2rexp (v3, v1, q, bp, n, t, 1, pp);
  TOOM8H_MUL_N_REC(pp, v0, v1, n + 1); /* A(-1/2)*B(-1/2)*2^. */
  TOOM8H_MUL_N_REC(r6, v2, v3, n + 1); /* A(+1/2)*B(+1/2)*2^. */
  mpn_toom_couple_handling (r6, 2 * n + 1, pp, sign, n, 1+half, half);

  /* $\pm1$ */
  sign = mpn_toom_eval_pm1 (v2, v0, p, ap, n, s, pp);
  if (q == 3)
    sign ^= mpn_toom_eval_dgr3_pm1 (v3, v1, bp, n, t, pp);
  else
    sign ^= mpn_toom_eval_pm1 (v3, v1, q, bp, n, t,    pp);
  TOOM8H_MUL_N_REC(pp, v0, v1, n + 1); /* A(-1)*B(-1) */
  TOOM8H_MUL_N_REC(r4, v2, v3, n + 1); /* A(1)*B(1) */
  mpn_toom_couple_handling (r4, 2 * n + 1, pp, sign, n, 0, 0);

  /* $\pm4$ */
  sign = mpn_toom_eval_pm2exp (v2, v0, p, ap, n, s, 2, pp) ^
	 mpn_toom_eval_pm2exp (v3, v1, q, bp, n, t, 2, pp);
  TOOM8H_MUL_N_REC(pp, v0, v1, n + 1); /* A(-4)*B(-4) */
  TOOM8H_MUL_N_REC(r2, v2, v3, n + 1); /* A(+4)*B(+4) */
  mpn_toom_couple_handling (r2, 2 * n + 1, pp, sign, n, 2, 4);

#undef v0
#undef v1
#undef v2
#undef v3

  /* A(0)*B(0) */
  TOOM8H_MUL_N_REC(pp, ap, bp, n);

  /* Infinity */
  if( half != 0) {
    if(s>t) {
      TOOM8H_MUL_REC(r0, ap + p * n, s, bp + q * n, t);
    } else {
      TOOM8H_MUL_REC(r0, bp + q * n, t, ap + p * n, s);
    };
  };

  mpn_toom_interpolate_16pts (pp, r1, r3, r5, r7, n, s+t, half, wsi);

  TMP_FREE;

#undef r0
#undef r1
#undef r2
#undef r3
#undef r4
#undef r5
#undef r6
#undef wsi
}
//////////////////////////////////////////////////////////////////////////
/* DO_addlsh2(d,a,b,n,cy) computes cy,{d,n} <- {a,n} + 4*(cy,{b,n}), it
can be used as DO_addlsh2(d,a,d,n,d[n]), for accumulation on {d,n+1}. */
#if HAVE_NATIVE_mpn_addlsh2_n
#define DO_addlsh2(d, a, b, n, cy)	\
do {					\
  (cy) <<= 2;				\
  (cy) += mpn_addlsh2_n(d, a, b, n);	\
} while (0)
#else
#if HAVE_NATIVE_mpn_addlsh_n
#define DO_addlsh2(d, a, b, n, cy)	\
do {					\
  (cy) <<= 2;				\
  (cy) += mpn_addlsh_n(d, a, b, n, 2);	\
} while (0)
#else
/* The following is not a general substitute for addlsh2.
It is correct if d == b, but it is not if d == a.	*/
#define DO_addlsh2(d, a, b, n, cy)	\
do {					\
  (cy) <<= 2;				\
  (cy) += mpn_lshift(d, b, n, 2);	\
  (cy) += mpn_add_n(d, d, a, n);	\
} while (0)
#endif
#endif

/* Evaluates a polynomial of degree 2 < k < GMP_NUMB_BITS, in the
points +2 and -2. */
__GMP_DECLSPEC
int
mpn_toom_eval_pm2(mp_ptr xp2, mp_ptr xm2, unsigned k,
mp_srcptr xp, mp_size_t n, mp_size_t hn, mp_ptr tp)
{
	int i;
	int neg;
	mp_limb_t cy;

	ASSERT(k >= 3);
	ASSERT(k < GMP_NUMB_BITS);

	ASSERT(hn > 0);
	ASSERT(hn <= n);

	/* The degree k is also the number of full-size coefficients, so
	* that last coefficient, of size hn, starts at xp + k*n. */

	cy = 0;
	DO_addlsh2(xp2, xp + (k - 2) * n, xp + k * n, hn, cy);
	if (hn != n)
		cy = mpn_add_1(xp2 + hn, xp + (k - 2) * n + hn, n - hn, cy);
	for (i = k - 4; i >= 0; i -= 2)
		DO_addlsh2(xp2, xp + i * n, xp2, n, cy);
	xp2[n] = cy;

	k--;

	cy = 0;
	DO_addlsh2(tp, xp + (k - 2) * n, xp + k * n, n, cy);
	for (i = k - 4; i >= 0; i -= 2)
		DO_addlsh2(tp, xp + i * n, tp, n, cy);
	tp[n] = cy;

	if (k & 1)
		ASSERT_NOCARRY(mpn_lshift(tp, tp, n + 1, 1));
	else
		ASSERT_NOCARRY(mpn_lshift(xp2, xp2, n + 1, 1));

	neg = (mpn_cmp(xp2, tp, n + 1) < 0) ? ~0 : 0;

#if HAVE_NATIVE_mpn_sumdiff_n
	if (neg)
		mpn_sumdiff_n(xp2, xm2, tp, xp2, n + 1);
	else
		mpn_sumdiff_n(xp2, xm2, xp2, tp, n + 1);
#else 
	if (neg)
		mpn_sub_n(xm2, tp, xp2, n + 1);
	else
		mpn_sub_n(xm2, xp2, tp, n + 1);

	mpn_add_n(xp2, xp2, tp, n + 1);
#endif

	ASSERT(xp2[n] < (1 << (k + 2)) - 1);
	ASSERT(xm2[n] < ((1 << (k + 3)) - 1 - (1 ^ k & 1)) / 3);

	neg ^= ((k & 1) - 1);

	return neg;
}

#undef DO_addlsh2

//////////////////////////////////////////////////////////////////////////

// k degree poly so have k+1 coeffs and first k are size n
// k>3 so we can do the first add unconditionally 
__GMP_DECLSPEC int	mpn_toom_eval_pm1(mp_ptr pp, mp_ptr mp, unsigned int k, mp_srcptr xp, mp_size_t n, mp_size_t m, mp_ptr tp)
{
	int isneg = 0; unsigned int i;

	ASSERT(k > 3); ASSERT(n >= m); ASSERT(m > 0); ASSERT_MPN(xp, n*k + m);
	//ASSERT_SPACE(pp,n+1);ASSERT_SPACE(mp,n+1);ASSERT_SPACE(tp,n+1);
	ASSERT(!MPN_OVERLAP_P(pp, n + 1, mp, n + 1)); ASSERT(!MPN_OVERLAP_P(pp, n + 1, xp, n*k + m)); ASSERT(!MPN_OVERLAP_P(pp, n + 1, tp, n + 1));
	ASSERT(!MPN_OVERLAP_P(mp, n + 1, xp, n*k + m)); ASSERT(!MPN_OVERLAP_P(xp, n*k + m, tp, n + 1));
#if ! HAVE_NATIVE_mpn_sumdiff_n
	ASSERT(!MPN_OVERLAP_P(mp, n + 1, tp, n + 1));
#endif
#if HAVE_NATIVE_mpn_addadd_n
	if (k == 4){ pp[n] = mpn_add_n(pp, xp, xp + 2 * n, n); tp[n] = mpn_add_n(tp, xp + n, xp + 3 * n, n); }
	else
		if (k == 5){ pp[n] = mpn_addadd_n(pp, xp, xp + 2 * n, xp + 4 * n, n); tp[n] = mpn_add_n(tp, xp + n, xp + 3 * n, n); }
		else
		{
			pp[n] = mpn_addadd_n(pp, xp, xp + 2 * n, xp + 4 * n, n); tp[n] = mpn_addadd_n(tp, xp + n, xp + 3 * n, xp + 5 * n, n);
			for (i = 7; i < k - 2; i += 4){ pp[n] += mpn_addadd_n(pp, pp, xp + (i - 1)*n, xp + (i + 1)*n, n); tp[n] += mpn_addadd_n(tp, tp, xp + i*n, xp + (i + 2)*n, n); }
			if (k % 4 == 3){ pp[n] += mpn_add_n(pp, pp, xp + (k - 1)*n, n); }
			if (k % 4 == 0){ pp[n] += mpn_add_n(pp, pp, xp + (k - 2)*n, n); tp[n] += mpn_add_n(tp, tp, xp + (k - 1)*n, n); }
			if (k % 4 == 1){ pp[n] += mpn_addadd_n(pp, pp, xp + (k - 3)*n, xp + (k - 1)*n, n); tp[n] += mpn_add_n(tp, tp, xp + (k - 2)*n, n); }
		}
	if (k % 2 == 0){ pp[n] += mpn_add(pp, pp, n, xp + k*n, m); }
	else{ tp[n] += mpn_add(tp, tp, n, xp + k*n, m); }
#else
	// pp is xp+0 xp+2n xp+4n xp+6n ... xp+jn where j<=k-1
	// mp is xp+1 xp+3n xp+5n xp+7n ... xp+jn where j<=k-1
	pp[n] = mpn_add_n(pp, xp, xp + 2 * n, n); tp[n] = mpn_add_n(tp, xp + n, xp + 3 * n, n);
	for (i = 5; i < k; i += 2)
	{
		pp[n] += mpn_add_n(pp, pp, xp + (i - 1)*n, n);
		tp[n] += mpn_add_n(tp, tp, xp + i*n, n);
	}
	if (k % 2 == 1)
	{
		pp[n] += mpn_add_n(pp, pp, xp + (k - 1)*n, n);
		tp[n] += mpn_add(tp, tp, n, xp + k*n, m);
	}
	else
	{
		pp[n] += mpn_add(pp, pp, n, xp + k*n, m);
	}
#endif
	if (mpn_cmp(tp, pp, n + 1)>0)
		isneg = -1;
#if HAVE_NATIVE_mpn_sumdiff_n
	if (isneg == 0)
	{
		mpn_sumdiff_n(pp, mp, pp, tp, n + 1);
	}
	else
	{
		mpn_sumdiff_n(pp, mp, tp, pp, n + 1);
	}
#else
	if (isneg == 0)
	{
		mpn_sub_n(mp, pp, tp, n + 1);
	}
	else
	{
		mpn_sub_n(mp, tp, pp, n + 1);
	}
	mpn_add_n(pp, pp, tp, n + 1);
#endif
	return isneg;
}

/* Evaluates a polynomial of degree k > 2, in the points +2^shift and -2^shift. */
__GMP_DECLSPEC
int
mpn_toom_eval_pm2exp(mp_ptr xp2, mp_ptr xm2, unsigned k,
mp_srcptr xp, mp_size_t n, mp_size_t hn, unsigned shift,
mp_ptr tp)
{
	unsigned i;
	int neg;
#if HAVE_NATIVE_mpn_addlsh_n
	mp_limb_t cy;
#endif

	ASSERT(k >= 3);
	ASSERT(shift*k < GMP_NUMB_BITS);

	ASSERT(hn > 0);
	ASSERT(hn <= n);

	/* The degree k is also the number of full-size coefficients, so
	* that last coefficient, of size hn, starts at xp + k*n. */

#if HAVE_NATIVE_mpn_addlsh_n
	xp2[n] = mpn_addlsh_n(xp2, xp, xp + 2 * n, n, 2 * shift);
	for (i = 4; i < k; i += 2)
		xp2[n] += mpn_addlsh_n(xp2, xp2, xp + i*n, n, i*shift);

	tp[n] = mpn_lshift(tp, xp + n, n, shift);
	for (i = 3; i < k; i += 2)
		tp[n] += mpn_addlsh_n(tp, tp, xp + i*n, n, i*shift);

	if (k & 1)
	{
		cy = mpn_addlsh_n(tp, tp, xp + k*n, hn, k*shift);
		MPN_INCR_U(tp + hn, n + 1 - hn, cy);
	}
	else
	{
		cy = mpn_addlsh_n(xp2, xp2, xp + k*n, hn, k*shift);
		MPN_INCR_U(xp2 + hn, n + 1 - hn, cy);
	}

#else /* !HAVE_NATIVE_mpn_addlsh_n */
	xp2[n] = mpn_lshift(tp, xp + 2 * n, n, 2 * shift);
	xp2[n] += mpn_add_n(xp2, xp, tp, n);
	for (i = 4; i < k; i += 2)
	{
		xp2[n] += mpn_lshift(tp, xp + i*n, n, i*shift);
		xp2[n] += mpn_add_n(xp2, xp2, tp, n);
	}

	tp[n] = mpn_lshift(tp, xp + n, n, shift);
	for (i = 3; i < k; i += 2)
	{
		tp[n] += mpn_lshift(xm2, xp + i*n, n, i*shift);
		tp[n] += mpn_add_n(tp, tp, xm2, n);
	}

	xm2[hn] = mpn_lshift(xm2, xp + k*n, hn, k*shift);
	if (k & 1)
		mpn_add(tp, tp, n + 1, xm2, hn + 1);
	else
		mpn_add(xp2, xp2, n + 1, xm2, hn + 1);
#endif /* !HAVE_NATIVE_mpn_addlsh_n */

	neg = (mpn_cmp(xp2, tp, n + 1) < 0) ? ~0 : 0;

#if HAVE_NATIVE_mpn_sumdiff_n
	if (neg)
		mpn_sumdiff_n(xp2, xm2, tp, xp2, n + 1);
	else
		mpn_sumdiff_n(xp2, xm2, xp2, tp, n + 1);
#else 
	if (neg)
		mpn_sub_n(xm2, tp, xp2, n + 1);
	else
		mpn_sub_n(xm2, xp2, tp, n + 1);

	mpn_add_n(xp2, xp2, tp, n + 1);
#endif

	/* FIXME: the following asserts are useless if (k+1)*shift >= GMP_LIMB_BITS */
	ASSERT((k + 1)*shift >= GMP_LIMB_BITS ||
		xp2[n] < ((CNST_LIMB(1) << ((k + 1)*shift)) - 1) / ((CNST_LIMB(1) << shift) - 1));
	ASSERT((k + 2)*shift >= GMP_LIMB_BITS ||
		xm2[n] < ((CNST_LIMB(1) << ((k + 2)*shift)) - ((k & 1) ? (CNST_LIMB(1) << shift) : 1)) / ((CNST_LIMB(1) << (2 * shift)) - 1));

	return neg;
}

__GMP_DECLSPEC
int
mpn_toom_eval_dgr3_pm1(mp_ptr xp1, mp_ptr xm1,
mp_srcptr xp, mp_size_t n, mp_size_t x3n, mp_ptr tp)
{
	int neg;

	ASSERT(x3n > 0);
	ASSERT(x3n <= n);

	xp1[n] = mpn_add_n(xp1, xp, xp + 2 * n, n);
	tp[n] = mpn_add(tp, xp + n, n, xp + 3 * n, x3n);

	neg = (mpn_cmp(xp1, tp, n + 1) < 0) ? ~0 : 0;

#if HAVE_NATIVE_mpn_sumdiff_n
	if (neg)
		mpn_sumdiff_n(xp1, xm1, tp, xp1, n + 1);
	else
		mpn_sumdiff_n(xp1, xm1, xp1, tp, n + 1);
#else
	if (neg)
		mpn_sub_n(xm1, tp, xp1, n + 1);
	else
		mpn_sub_n(xm1, xp1, tp, n + 1);

	mpn_add_n(xp1, xp1, tp, n + 1);
#endif

	ASSERT(xp1[n] <= 3);
	ASSERT(xm1[n] <= 1);

	return neg;
}
//////////////////////////////////////////////////////////////////////////
#if GMP_NUMB_BITS < 29
#error Not implemented: Both sublsh_n(,,,28) should be corrected; r2 and r5 need one more LIMB.
#endif

#if GMP_NUMB_BITS < 28
#error Not implemented: divexact_by188513325 and _by182712915 will not work.
#endif


#if HAVE_NATIVE_mpn_sublsh_n
#define DO_mpn_sublsh_n(dst,src,n,s,ws) mpn_sublsh_n(dst,dst,src,n,s)
#else
static
__GMP_DECLSPEC
mp_limb_t
DO_mpn_sublsh_n(mp_ptr dst, mp_srcptr src, mp_size_t n, unsigned int s, mp_ptr ws)
{
#if USE_MUL_1 && 0
	return mpn_submul_1(dst, src, n, CNST_LIMB(1) << (s));
#else
	mp_limb_t __cy;
	__cy = mpn_lshift(ws, src, n, s);
	return    __cy + mpn_sub_n(dst, dst, ws, n);
#endif
}
#endif


#if HAVE_NATIVE_mpn_subrsh
#define DO_mpn_subrsh(dst,nd,src,ns,s,ws) mpn_subrsh(dst,nd,src,ns,s)
#else
/* FIXME: This is not a correct definition, it assumes no carry */
#define DO_mpn_subrsh(dst,nd,src,ns,s,ws)				\
do {									\
  mp_limb_t __cy;							\
  MPN_DECR_U (dst, nd, src[0] >> s);					\
  __cy = DO_mpn_sublsh_n (dst, src + 1, ns - 1, GMP_NUMB_BITS - s, ws);	\
  MPN_DECR_U (dst + ns - 1, nd - ns + 1, __cy);				\
} while (0)
#endif


/* FIXME: tuneup should decide the best variant */
#ifndef AORSMUL_FASTER_AORS_AORSLSH
#define AORSMUL_FASTER_AORS_AORSLSH 1
#endif
#ifndef AORSMUL_FASTER_AORS_2AORSLSH
#define AORSMUL_FASTER_AORS_2AORSLSH 1
#endif
#ifndef AORSMUL_FASTER_2AORSLSH
#define AORSMUL_FASTER_2AORSLSH 1
#endif
#ifndef AORSMUL_FASTER_3AORSLSH
#define AORSMUL_FASTER_3AORSLSH 1
#endif

#if GMP_NUMB_BITS < 43
#define BIT_CORRECTION 1
#define CORRECTION_BITS GMP_NUMB_BITS
#else
#define BIT_CORRECTION 0
#define CORRECTION_BITS 0
#endif

#define BINVERT_9 \
  ((((GMP_NUMB_MAX / 9) << (6 - GMP_NUMB_BITS % 6)) * 8 & GMP_NUMB_MAX) | 0x39)

#define BINVERT_255 \
  (GMP_NUMB_MAX - ((GMP_NUMB_MAX / 255) << (8 - GMP_NUMB_BITS % 8)))

/* FIXME: find some more general expressions for inverses */
#if GMP_LIMB_BITS == 32
#define BINVERT_2835  (GMP_NUMB_MASK &		CNST_LIMB(0x53E3771B))
#define BINVERT_42525 (GMP_NUMB_MASK &		CNST_LIMB(0x9F314C35))
#define BINVERT_182712915 (GMP_NUMB_MASK &	CNST_LIMB(0x550659DB))
#define BINVERT_188513325 (GMP_NUMB_MASK &	CNST_LIMB(0xFBC333A5))
#define BINVERT_255x182712915L (GMP_NUMB_MASK &	CNST_LIMB(0x6FC4CB25))
#define BINVERT_255x188513325L (GMP_NUMB_MASK &	CNST_LIMB(0x6864275B))
#if GMP_NAIL_BITS == 0
#define BINVERT_255x182712915H CNST_LIMB(0x1B649A07)
#define BINVERT_255x188513325H CNST_LIMB(0x06DB993A)
#else /* GMP_NAIL_BITS != 0 */
#define BINVERT_255x182712915H \
  (GMP_NUMB_MASK & CNST_LIMB((0x1B649A07<<GMP_NAIL_BITS) | (0x6FC4CB25>>GMP_NUMB_BITS)))
#define BINVERT_255x188513325H \
  (GMP_NUMB_MASK & CNST_LIMB((0x06DB993A<<GMP_NAIL_BITS) | (0x6864275B>>GMP_NUMB_BITS)))
#endif
#else
#if GMP_LIMB_BITS == 64
#define BINVERT_2835  (GMP_NUMB_MASK &	CNST_LIMB(0x938CC70553E3771B))
#define BINVERT_42525 (GMP_NUMB_MASK &	CNST_LIMB(0xE7B40D449F314C35))
#define BINVERT_255x182712915  (GMP_NUMB_MASK &	CNST_LIMB(0x1B649A076FC4CB25))
#define BINVERT_255x188513325  (GMP_NUMB_MASK &	CNST_LIMB(0x06DB993A6864275B))
#endif
#endif

#ifndef mpn_divexact_by255
#if GMP_NUMB_BITS % 8 == 0
#define mpn_divexact_by255(dst,src,size) \
  (255 & 1 * mpn_divexact_byfobm1(dst, src, size, 255, __GMP_CAST (mp_limb_t, GMP_NUMB_MASK / 255)))
#else
#if HAVE_NATIVE_mpn_pi1_bdiv_q_1
#define mpn_divexact_by255(dst,src,size) mpn_pi1_bdiv_q_1(dst,src,size,CNST_LIMB(255),BINVERT_255,0)
#else
#define mpn_divexact_by255(dst,src,size) mpn_divexact_1(dst,src,size,CNST_LIMB(255))
#endif
#endif
#endif

#ifndef mpn_divexact_by255x4
#if HAVE_NATIVE_mpn_pi1_bdiv_q_1
#define mpn_divexact_by255x4(dst,src,size) mpn_pi1_bdiv_q_1(dst,src,size,CNST_LIMB(255),BINVERT_255,2)
#else
#define mpn_divexact_by255x4(dst,src,size) mpn_divexact_1(dst,src,size,CNST_LIMB(255)<<2)
#endif
#endif

#ifndef mpn_divexact_by9x16
#if HAVE_NATIVE_mpn_pi1_bdiv_q_1
#define mpn_divexact_by9x16(dst,src,size) mpn_pi1_bdiv_q_1(dst,src,size,CNST_LIMB(9),BINVERT_9,4)
#else
#define mpn_divexact_by9x16(dst,src,size) mpn_divexact_1(dst,src,size,CNST_LIMB(9)<<4)
#endif
#endif

#ifndef mpn_divexact_by42525x16
#if HAVE_NATIVE_mpn_pi1_bdiv_q_1 && defined(BINVERT_42525)
#define mpn_divexact_by42525x16(dst,src,size) mpn_pi1_bdiv_q_1(dst,src,size,CNST_LIMB(42525),BINVERT_42525,4)
#else
#define mpn_divexact_by42525x16(dst,src,size) mpn_divexact_1(dst,src,size,CNST_LIMB(42525)<<4)
#endif
#endif

#ifndef mpn_divexact_by2835x64
#if HAVE_NATIVE_mpn_pi1_bdiv_q_1 && defined(BINVERT_2835)
#define mpn_divexact_by2835x64(dst,src,size) mpn_pi1_bdiv_q_1(dst,src,size,CNST_LIMB(2835),BINVERT_2835,6)
#else
#define mpn_divexact_by2835x64(dst,src,size) mpn_divexact_1(dst,src,size,CNST_LIMB(2835)<<6)
#endif
#endif

#ifndef  mpn_divexact_by255x182712915
#if GMP_NUMB_BITS < 36
#if HAVE_NATIVE_mpn_bdiv_q_2_pi2 && defined(BINVERT_255x182712915H)
/* FIXME: use mpn_bdiv_q_2_pi2 */
#endif
#if HAVE_NATIVE_mpn_pi1_bdiv_q_1 && defined(BINVERT_182712915)
#define mpn_divexact_by255x182712915(dst,src,size)				\
  do {										\
    mpn_pi1_bdiv_q_1(dst,src,size,CNST_LIMB(182712915),BINVERT_182712915,0);	\
    mpn_divexact_by255(dst,dst,size);						\
    } while(0)
#else
#define mpn_divexact_by255x182712915(dst,src,size)	\
  do {							\
    mpn_divexact_1(dst,src,size,CNST_LIMB(182712915));	\
    mpn_divexact_by255(dst,dst,size);			\
    } while(0)
#endif
#else /* GMP_NUMB_BITS > 35 */
#if HAVE_NATIVE_mpn_pi1_bdiv_q_1 && defined(BINVERT_255x182712915)
#define mpn_divexact_by255x182712915(dst,src,size) \
  mpn_pi1_bdiv_q_1(dst,src,size,255*CNST_LIMB(182712915),BINVERT_255x182712915,0)
#else
#define mpn_divexact_by255x182712915(dst,src,size) mpn_divexact_1(dst,src,size,255*CNST_LIMB(182712915))
#endif
#endif /* GMP_NUMB_BITS >?< 36 */
#endif

#ifndef  mpn_divexact_by255x188513325
#if GMP_NUMB_BITS < 36
#if HAVE_NATIVE_mpn_bdiv_q_1_pi2 && defined(BINVERT_255x188513325H)
/* FIXME: use mpn_bdiv_q_1_pi2 */
#endif
#if HAVE_NATIVE_mpn_pi1_bdiv_q_1 && defined(BINVERT_188513325)
#define mpn_divexact_by255x188513325(dst,src,size)			\
  do {									\
    mpn_pi1_bdiv_q_1(dst,src,size,CNST_LIMB(188513325),BINVERT_188513325,0);	\
    mpn_divexact_by255(dst,dst,size);					\
    } while(0)
#else
#define mpn_divexact_by255x188513325(dst,src,size)	\
  do {							\
    mpn_divexact_1(dst,src,size,CNST_LIMB(188513325));	\
    mpn_divexact_by255(dst,dst,size);			\
    } while(0)
#endif
#else /* GMP_NUMB_BITS > 35 */
#if HAVE_NATIVE_mpn_pi1_bdiv_q_1 && defined(BINVERT_255x188513325)
#define mpn_divexact_by255x188513325(dst,src,size) \
  mpn_pi1_bdiv_q_1(dst,src,size,255*CNST_LIMB(188513325),BINVERT_255x188513325,0)
#else
#define mpn_divexact_by255x188513325(dst,src,size) mpn_divexact_1(dst,src,size,255*CNST_LIMB(188513325))
#endif
#endif /* GMP_NUMB_BITS >?< 36 */
#endif

/* Interpolation for Toom-8.5 (or Toom-8), using the evaluation
points: infinity(8.5 only), +-8, +-4, +-2, +-1, +-1/4, +-1/2,
+-1/8, 0. More precisely, we want to compute
f(2^(GMP_NUMB_BITS * n)) for a polynomial f of degree 15 (or
14), given the 16 (rsp. 15) values:

r0 = limit at infinity of f(x) / x^7,
r1 = f(8),f(-8),
r2 = f(4),f(-4),
r3 = f(2),f(-2),
r4 = f(1),f(-1),
r5 = f(1/4),f(-1/4),
r6 = f(1/2),f(-1/2),
r7 = f(1/8),f(-1/8),
r8 = f(0).

All couples of the form f(n),f(-n) must be already mixed with
toom_couple_handling(f(n),...,f(-n),...)

The result is stored in {pp, spt + 7*n (or 8*n)}.
At entry, r8 is stored at {pp, 2n},
r6 is stored at {pp + 3n, 3n + 1}.
r4 is stored at {pp + 7n, 3n + 1}.
r2 is stored at {pp +11n, 3n + 1}.
r0 is stored at {pp +15n, spt}.

The other values are 3n+1 limbs each (with most significant limbs small).

Negative intermediate results are stored two-complemented.
Inputs are destroyed.
*/

__GMP_DECLSPEC
void
mpn_toom_interpolate_16pts(mp_ptr pp, mp_ptr r1, mp_ptr r3, mp_ptr r5, mp_ptr r7,
mp_size_t n, mp_size_t spt, int half, mp_ptr wsi)
{
	mp_limb_t cy;
	mp_size_t n3;
	mp_size_t n3p1;
	n3 = 3 * n;
	n3p1 = n3 + 1;

#define   r6    (pp + n3)			/* 3n+1 */
#define   r4    (pp + 7 * n)			/* 3n+1 */
#define   r2    (pp +11 * n)			/* 3n+1 */
#define   r0    (pp +15 * n)			/* s+t <= 2*n */

	ASSERT(spt <= 2 * n);
	/******************************* interpolation *****************************/
	if (half != 0) {
		cy = mpn_sub_n(r4, r4, r0, spt);
		MPN_DECR_U(r4 + spt, n3p1 - spt, cy);

		cy = DO_mpn_sublsh_n(r3, r0, spt, 14, wsi);
		MPN_DECR_U(r3 + spt, n3p1 - spt, cy);
		DO_mpn_subrsh(r6, n3p1, r0, spt, 2, wsi);

		cy = DO_mpn_sublsh_n(r2, r0, spt, 28, wsi);
		MPN_DECR_U(r2 + spt, n3p1 - spt, cy);
		DO_mpn_subrsh(r5, n3p1, r0, spt, 4, wsi);

		cy = DO_mpn_sublsh_n(r1 + BIT_CORRECTION, r0, spt, 42 - CORRECTION_BITS, wsi);
		//MPN_DECR_U (r1 + spt + BIT_CORRECTION, n3p1 - spt - BIT_CORRECTION, cy);
#if BIT_CORRECTION
		cy = mpn_sub_1(r1 + spt + BIT_CORRECTION, r1 + spt + BIT_CORRECTION,
			n3p1 - spt - BIT_CORRECTION, cy);
		ASSERT(BIT_CORRECTION > 0 || cy == 0);
		/* FIXME: assumes r7[n3p1] is writable (it is if r5 follows). */
		cy = r7[n3p1];
		r7[n3p1] = 0x80;
#else
		MPN_DECR_U(r1 + spt + BIT_CORRECTION, n3p1 - spt - BIT_CORRECTION, cy);
#endif
		DO_mpn_subrsh(r7, n3p1 + BIT_CORRECTION, r0, spt, 6, wsi);
#if BIT_CORRECTION
		/* FIXME: assumes r7[n3p1] is writable. */
		ASSERT(BIT_CORRECTION > 0 || r7[n3p1] == 0x80);
		r7[n3p1] = cy;
#endif
	};

	r5[n3] -= DO_mpn_sublsh_n(r5 + n, pp, 2 * n, 28, wsi);
	DO_mpn_subrsh(r2 + n, 2 * n + 1, pp, 2 * n, 4, wsi);

#if HAVE_NATIVE_mpn_sumdiff_n
	mpn_sumdiff_n(r2, wsi, r5, r2, n3p1);
	MP_PTR_SWAP(r5, wsi);
#else
	mpn_sub_n(wsi, r5, r2, n3p1); /* can be negative */
	ASSERT_NOCARRY(mpn_add_n(r2, r2, r5, n3p1));
	MP_PTR_SWAP(r5, wsi);
#endif

	r6[n3] -= DO_mpn_sublsh_n(r6 + n, pp, 2 * n, 14, wsi);
	DO_mpn_subrsh(r3 + n, 2 * n + 1, pp, 2 * n, 2, wsi);

#if HAVE_NATIVE_mpn_sumdiff_n
	mpn_sumdiff_n(wsi, r6, r6, r3, n3p1);
	MP_PTR_SWAP(r3, wsi);
#else
	ASSERT_NOCARRY(mpn_add_n(wsi, r3, r6, n3p1));
	mpn_sub_n(r6, r6, r3, n3p1); /* can be negative */
	MP_PTR_SWAP(r3, wsi);
#endif

	r7[n3] -= DO_mpn_sublsh_n(r7 + n + BIT_CORRECTION, pp, 2 * n, 42 - CORRECTION_BITS, wsi)
		* (1 - BIT_CORRECTION); /* if BIT_CORRECTION != 0, discard the carry. */
#if BIT_CORRECTION
	MPN_DECR_U(r1 + n, 2 * n + 1, pp[0] >> 6);
	cy = DO_mpn_sublsh_n(r1 + n, pp + 1, 2 * n - 1, GMP_NUMB_BITS - 6, wsi);
	cy = mpn_sub_1(r1 + 3 * n - 1, r1 + 3 * n - 1, 2, cy);
	ASSERT(BIT_CORRECTION > 0 || cy != 0);
#else
	DO_mpn_subrsh(r1 + n, 2 * n + 1, pp, 2 * n, 6, wsi);
#endif

#if HAVE_NATIVE_mpn_sumdiff_n
	mpn_sumdiff_n(r1, wsi, r7, r1, n3p1);
	MP_PTR_SWAP(r7, wsi);
#else
	mpn_sub_n(wsi, r7, r1, n3p1); /* can be negative */
	mpn_add_n(r1, r1, r7, n3p1);  /* if BIT_CORRECTION != 0, can give a carry. */
	MP_PTR_SWAP(r7, wsi);
#endif

	r4[n3] -= mpn_sub_n(r4 + n, r4 + n, pp, 2 * n);

#if AORSMUL_FASTER_2AORSLSH
	mpn_submul_1(r5, r6, n3p1, 1028); /* can be negative */
#else
	DO_mpn_sublsh_n(r5, r6, n3p1, 2, wsi); /* can be negative */
	DO_mpn_sublsh_n(r5, r6, n3p1, 10, wsi); /* can be negative */
#endif

	mpn_submul_1(r7, r5, n3p1, 1300); /* can be negative */
#if AORSMUL_FASTER_3AORSLSH
	mpn_submul_1(r7, r6, n3p1, 1052688); /* can be negative */
#else
	DO_mpn_sublsh_n(r7, r6, n3p1, 4, wsi); /* can be negative */
	DO_mpn_sublsh_n(r7, r6, n3p1, 12, wsi); /* can be negative */
	DO_mpn_sublsh_n(r7, r6, n3p1, 20, wsi); /* can be negative */
#endif
	mpn_divexact_by255x188513325(r7, r7, n3p1);

	mpn_submul_1(r5, r7, n3p1, 12567555); /* can be negative */
	/* A division by 2835x64 followsi. Warning: the operand can be negative! */
	mpn_divexact_by2835x64(r5, r5, n3p1);
	if ((r5[n3] & (GMP_NUMB_MAX << (GMP_NUMB_BITS - 7))) != 0)
		r5[n3] |= (GMP_NUMB_MAX << (GMP_NUMB_BITS - 6));

#if AORSMUL_FASTER_AORS_AORSLSH
	mpn_submul_1(r6, r7, n3p1, 4095); /* can be negative */
#else
	mpn_add_n(r6, r6, r7, n3p1); /* can give a carry */
	DO_mpn_sublsh_n(r6, r7, n3p1, 12, wsi); /* can be negative */
#endif
#if AORSMUL_FASTER_2AORSLSH
	mpn_addmul_1(r6, r5, n3p1, 240); /* can be negative */
#else
	DO_mpn_addlsh_n(r6, r5, n3p1, 8, wsi); /* can give a carry */
	DO_mpn_sublsh_n(r6, r5, n3p1, 4, wsi); /* can be negative */
#endif
	/* A division by 255x4 followsi. Warning: the operand can be negative! */
	mpn_divexact_by255x4(r6, r6, n3p1);
	if ((r6[n3] & (GMP_NUMB_MAX << (GMP_NUMB_BITS - 3))) != 0)
		r6[n3] |= (GMP_NUMB_MAX << (GMP_NUMB_BITS - 2));

	ASSERT_NOCARRY(DO_mpn_sublsh_n(r3, r4, n3p1, 7, wsi));

	ASSERT_NOCARRY(DO_mpn_sublsh_n(r2, r4, n3p1, 13, wsi));
	ASSERT_NOCARRY(mpn_submul_1(r2, r3, n3p1, 400));

	/* If GMP_NUMB_BITS < 42 next operations on r1 can give a carry!*/
	DO_mpn_sublsh_n(r1, r4, n3p1, 19, wsi);
	mpn_submul_1(r1, r2, n3p1, 1428);
	mpn_submul_1(r1, r3, n3p1, 112896);
	mpn_divexact_by255x182712915(r1, r1, n3p1);

	ASSERT_NOCARRY(mpn_submul_1(r2, r1, n3p1, 15181425));
	mpn_divexact_by42525x16(r2, r2, n3p1);

#if AORSMUL_FASTER_AORS_2AORSLSH
	ASSERT_NOCARRY(mpn_submul_1(r3, r1, n3p1, 3969));
#else
	ASSERT_NOCARRY(mpn_sub_n(r3, r3, r1, n3p1));
	ASSERT_NOCARRY(DO_mpn_addlsh_n(r3, r1, n3p1, 7, wsi));
	ASSERT_NOCARRY(DO_mpn_sublsh_n(r3, r1, n3p1, 12, wsi));
#endif
	ASSERT_NOCARRY(mpn_submul_1(r3, r2, n3p1, 900));
	mpn_divexact_by9x16(r3, r3, n3p1);

	ASSERT_NOCARRY(mpn_sub_n(r4, r4, r1, n3p1));
	ASSERT_NOCARRY(mpn_sub_n(r4, r4, r3, n3p1));
	ASSERT_NOCARRY(mpn_sub_n(r4, r4, r2, n3p1));

	mpn_add_n(r6, r2, r6, n3p1);
	ASSERT_NOCARRY(mpn_rshift(r6, r6, n3p1, 1));
	ASSERT_NOCARRY(mpn_sub_n(r2, r2, r6, n3p1));

	mpn_sub_n(r5, r3, r5, n3p1);
	ASSERT_NOCARRY(mpn_rshift(r5, r5, n3p1, 1));
	ASSERT_NOCARRY(mpn_sub_n(r3, r3, r5, n3p1));

	mpn_add_n(r7, r1, r7, n3p1);
	ASSERT_NOCARRY(mpn_rshift(r7, r7, n3p1, 1));
	ASSERT_NOCARRY(mpn_sub_n(r1, r1, r7, n3p1));

	/* last interpolation steps... */
	/* ... could be mixed with recomposition
	||H-r7|M-r7|L-r7|   ||H-r5|M-r5|L-r5|
	*/

	/***************************** recomposition *******************************/
	/*
	pp[] prior to operations:
	|M r0|L r0|___||H r2|M r2|L r2|___||H r4|M r4|L r4|___||H r6|M r6|L r6|____|H_r8|L r8|pp

	summation scheme for remaining operations:
	|__16|n_15|n_14|n_13|n_12|n_11|n_10|n__9|n__8|n__7|n__6|n__5|n__4|n__3|n__2|n___|n___|pp
	|M r0|L r0|___||H r2|M r2|L r2|___||H r4|M r4|L r4|___||H r6|M r6|L r6|____|H_r8|L r8|pp
	||H r1|M r1|L r1|   ||H r3|M r3|L r3|   ||H_r5|M_r5|L_r5|   ||H r7|M r7|L r7|
	*/

	cy = mpn_add_n(pp + n, pp + n, r7, n);
	cy = mpn_add_1(pp + 2 * n, r7 + n, n, cy);
#if HAVE_NATIVE_mpn_add_nc
	cy = r7[n3] + mpn_add_nc(pp + n3, pp + n3, r7 + 2 * n, n, cy);
#else
	MPN_INCR_U(r7 + 2 * n, n + 1, cy);
	cy = r7[n3] + mpn_add_n(pp + n3, pp + n3, r7 + 2 * n, n);
#endif
	MPN_INCR_U(pp + 4 * n, 2 * n + 1, cy);

	pp[2 * n3] += mpn_add_n(pp + 5 * n, pp + 5 * n, r5, n);
	cy = mpn_add_1(pp + 2 * n3, r5 + n, n, pp[2 * n3]);
#if HAVE_NATIVE_mpn_add_nc
	cy = r5[n3] + mpn_add_nc(pp + 7 * n, pp + 7 * n, r5 + 2 * n, n, cy);
#else
	MPN_INCR_U(r5 + 2 * n, n + 1, cy);
	cy = r5[n3] + mpn_add_n(pp + 7 * n, pp + 7 * n, r5 + 2 * n, n);
#endif
	MPN_INCR_U(pp + 8 * n, 2 * n + 1, cy);

	pp[10 * n] += mpn_add_n(pp + 9 * n, pp + 9 * n, r3, n);
	cy = mpn_add_1(pp + 10 * n, r3 + n, n, pp[10 * n]);
#if HAVE_NATIVE_mpn_add_nc
	cy = r3[n3] + mpn_add_nc(pp + 11 * n, pp + 11 * n, r3 + 2 * n, n, cy);
#else
	MPN_INCR_U(r3 + 2 * n, n + 1, cy);
	cy = r3[n3] + mpn_add_n(pp + 11 * n, pp + 11 * n, r3 + 2 * n, n);
#endif
	MPN_INCR_U(pp + 12 * n, 2 * n + 1, cy);

	pp[14 * n] += mpn_add_n(pp + 13 * n, pp + 13 * n, r1, n);
	if (half) {
		cy = mpn_add_1(pp + 14 * n, r1 + n, n, pp[14 * n]);
#if HAVE_NATIVE_mpn_add_nc
		if (LIKELY(spt > n)) {
			cy = r1[n3] + mpn_add_nc(pp + 15 * n, pp + 15 * n, r1 + 2 * n, n, cy);
			MPN_INCR_U(pp + 16 * n, spt - n, cy);
		}
		else {
			ASSERT_NOCARRY(mpn_add_nc(pp + 15 * n, pp + 15 * n, r1 + 2 * n, spt, cy));
		}
#else
		MPN_INCR_U(r1 + 2 * n, n + 1, cy);
		if (LIKELY(spt > n)) {
			cy = r1[n3] + mpn_add_n(pp + 15 * n, pp + 15 * n, r1 + 2 * n, n);
			MPN_INCR_U(pp + 16 * n, spt - n, cy);
		}
		else {
			ASSERT_NOCARRY(mpn_add_n(pp + 15 * n, pp + 15 * n, r1 + 2 * n, spt));
		}
#endif
	}
	else {
		ASSERT_NOCARRY(mpn_add_1(pp + 14 * n, r1 + n, spt, pp[14 * n]));
	}

#undef   r0
#undef   r2
#undef   r4
#undef   r6
}
//////////////////////////////////////////////////////////////////////////
__GMP_DECLSPEC
void
mpn_divexact_1(mp_ptr dst, mp_srcptr src, mp_size_t size, mp_limb_t divisor)
{
	mp_size_t  i;
	mp_limb_t  c, h, l, ls, s, s_next, inverse, dummy;
	unsigned   shift;

	ASSERT(size >= 1);
	ASSERT(divisor != 0);
	ASSERT(MPN_SAME_OR_SEPARATE_P(dst, src, size));
	ASSERT_MPN(src, size);
	ASSERT_LIMB(divisor);

	s = src[0];

	if (size == 1)
	{
		dst[0] = s / divisor;
		return;
	}

	if ((divisor & 1) == 0)
	{
		count_trailing_zeros(shift, divisor);
		divisor >>= shift;
	}
	else
		shift = 0;

	modlimb_invert(inverse, divisor);
	divisor <<= GMP_NAIL_BITS;

	if (shift != 0)
	{
		c = 0;
		i = 0;
		size--;

		do
		{
			s_next = src[i + 1];
			ls = ((s >> shift) | (s_next << (GMP_NUMB_BITS - shift))) & GMP_NUMB_MASK;
			s = s_next;

			SUBC_LIMB(c, l, ls, c);

			l = (l * inverse) & GMP_NUMB_MASK;
			dst[i] = l;

			umul_ppmm(h, dummy, l, divisor);
			c += h;

			i++;
		} while (i < size);

		ls = s >> shift;
		l = ls - c;
		l = (l * inverse) & GMP_NUMB_MASK;
		dst[i] = l;
	}
	else
	{
		l = (s * inverse) & GMP_NUMB_MASK;
		dst[0] = l;
		i = 1;
		c = 0;

		do
		{
			umul_ppmm(h, dummy, l, divisor);
			c += h;

			s = src[i];
			SUBC_LIMB(c, l, s, c);

			l = (l * inverse) & GMP_NUMB_MASK;
			dst[i] = l;
			i++;
		} while (i < size);
	}
}
