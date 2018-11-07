/* mpn_divrem_1 -- mpn by limb division.

Copyright 1991, 1993, 1994, 1996, 1998, 1999, 2000, 2002, 2003 Free Software
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

#ifndef DIVREM_1_NORM_THRESHOLD
#define DIVREM_1_NORM_THRESHOLD  0
#endif
#ifndef DIVREM_1_UNNORM_THRESHOLD
#define DIVREM_1_UNNORM_THRESHOLD  0
#endif


/* modlimb_invert_table[i] is the multiplicative inverse of 2*i+1 mod 256,
 ie. (modlimb_invert_table[i] * (2*i+1)) % 256 == 1 */


__GMP_DECLSPEC mp_limb_t getValFrom_invert_table(mp_limb_t idx)
{
	if (idx >= 128)return 0;
	const unsigned char  modlimb_invert_table[128] = {
		0x01, 0xAB, 0xCD, 0xB7, 0x39, 0xA3, 0xC5, 0xEF,
		0xF1, 0x1B, 0x3D, 0xA7, 0x29, 0x13, 0x35, 0xDF,
		0xE1, 0x8B, 0xAD, 0x97, 0x19, 0x83, 0xA5, 0xCF,
		0xD1, 0xFB, 0x1D, 0x87, 0x09, 0xF3, 0x15, 0xBF,
		0xC1, 0x6B, 0x8D, 0x77, 0xF9, 0x63, 0x85, 0xAF,
		0xB1, 0xDB, 0xFD, 0x67, 0xE9, 0xD3, 0xF5, 0x9F,
		0xA1, 0x4B, 0x6D, 0x57, 0xD9, 0x43, 0x65, 0x8F,
		0x91, 0xBB, 0xDD, 0x47, 0xC9, 0xB3, 0xD5, 0x7F,
		0x81, 0x2B, 0x4D, 0x37, 0xB9, 0x23, 0x45, 0x6F,
		0x71, 0x9B, 0xBD, 0x27, 0xA9, 0x93, 0xB5, 0x5F,
		0x61, 0x0B, 0x2D, 0x17, 0x99, 0x03, 0x25, 0x4F,
		0x51, 0x7B, 0x9D, 0x07, 0x89, 0x73, 0x95, 0x3F,
		0x41, 0xEB, 0x0D, 0xF7, 0x79, 0xE3, 0x05, 0x2F,
		0x31, 0x5B, 0x7D, 0xE7, 0x69, 0x53, 0x75, 0x1F,
		0x21, 0xCB, 0xED, 0xD7, 0x59, 0xC3, 0xE5, 0x0F,
		0x11, 0x3B, 0x5D, 0xC7, 0x49, 0x33, 0x55, 0xFF
	};
	return (mp_limb_t)(modlimb_invert_table[idx]);
}


/***************************************************************/

/* set to 1 = store or 0 = not store */
#define STORE_QUOTIENT 1
/* set to 0 = udiv  1 = gmp-preinv   2 = barrett */
#define UDIV_METHOD 1

#if UDIV_NEEDS_NORMALIZATION == 1 || UDIV_METHOD == 1
#define NORMALIZE 1
#else
#define NORMALIZE 0
#endif

#if UDIV_METHOD == 0
#define UDIV(q, r, h, l, d, i) udiv_qrnnd(q, r, h, l, d)
#endif

#if UDIV_METHOD == 1
#define UDIV udiv_qrnnd_preinv
#endif

#if UDIV_METHOD == 2
#define UDIV udiv_qrnnd_barrett
#endif

/***************************************************************/

__GMP_DECLSPEC
mp_limb_t mpn_rsh_divrem_hensel_qr_1_1(mp_ptr qp, mp_srcptr xp,
                                       mp_size_t n, mp_limb_t d, int s, mp_limb_t cin)
{
    mp_size_t j;
    mp_limb_t c, h, q, dummy, h1, t, m, qo;
    
    ASSERT(n > 0);
    ASSERT(d%2 == 1);
    ASSERT_MPN(xp, n);
    ASSERT(MPN_SAME_OR_SEPARATE_P(qp, xp, n));
    ASSERT(s >= 0);
    
    modlimb_invert(m, d); /*should we allow s = 0 ?? */
    h1 = xp[0];
    c = 0;
    h = cin;
    t = h + c;
    
    if (t > h1)
    {
        h1 = h1 - t;
        c = 1;
    }
    else
    {
        h1 = h1 - t;
        c = 0;
    }
    
    q = h1*m;
    qo = q>>s;
    umul_ppmm(h, dummy, q, d);
    
    for (j = 1; j <= n - 1; j++)
    {
        h1 = xp[j];
        t = h + c;
        if (t > h1)
        {
            h1 = h1 - t;
            c = 1;
        }
        else
        {
            h1 = h1 - t;
            c = 0;
        }
        
        q = h1*m;
        qo = qo | (q<<(GMP_LIMB_BITS - 1 - s)<<1);
        qp[j - 1] = qo;
        qo = q>>s;
        umul_ppmm(h, dummy, q, d);
        
        ASSERT(dummy == h1);
    }
    
    qp[n - 1] = qo;
    return h + c;
}

/*
 using a two limb inverse of a one limb divisor
 (xp,n) = (qp,n)*d - ret*B^n and 0 <= ret < d
 */
__GMP_DECLSPEC
mp_limb_t mpn_rsh_divrem_hensel_qr_1_2(mp_ptr qp, mp_srcptr xp,
                                       mp_size_t n, mp_limb_t d, int s, mp_limb_t cin)
{
    mp_size_t j;
    mp_limb_t c, h, q, dummy, h1, t, ml, mh, xl, xh, ql, qh, qo;
    
    ASSERT(n >= 2);
    ASSERT_MPN(xp, n);
    ASSERT(MPN_SAME_OR_SEPARATE_P(qp, xp, n));
    ASSERT(d%2 == 1);
    ASSERT(s >= 0);
    
    modlimb_invert(ml, d);
    umul_ppmm(h, dummy, d, ml);
    
    ASSERT(dummy == 1);
    
    h = -h;
    mh = ml*h; /* (mh, ml) is our two limb inverse */
    h1 = xp[0];
    h = cin;
    c = 0;
    t = h + c;
    
    if (t > h1)
    {
        h1 = h1 - t;
        c = 1;
    } else
    {
        h1 = h1 - t;
        c = 0;
    }
    
    q = h1*ml;
    qo = q>>s;
    umul_ppmm(h, dummy, q, d);
    
    for (j = 1; j + 1 <= n - 1; j += 2)
    {
        xl = xp[j];
        xh = xp[j + 1];
        t = h + c;
        
        if (xh == 0 && t > xl)
            c = 1;
        else
            c = 0;
        
        sub_ddmmss(xh, xl, xh, xl, 0, t);
        umul_ppmm(qh, ql, xl, ml);
        qh = qh + xh*ml + xl*mh;
        
        qo = qo|(ql<<(GMP_LIMB_BITS - 1 - s)<<1);
        qp[j - 1] = qo;
        qo = ql>>s;
        
        qo = qo | (qh<<(GMP_LIMB_BITS - 1 - s)<<1);
        qp[j + 1 - 1] = qo;
        qo = qh>>s;
        
        umul_ppmm(h, h1, qh, d);
        
        if (h1 > xh)
            h++;
    }
    
    if (j <= n-1)
    {
        h1 = xp[j];
        t = h + c;
        
        if (t > h1)
        {
            h1 = h1 - t;
            c = 1;
        }
        else
        {
            h1 = h1 - t;
            c = 0;
        }
        
        q = h1*ml;
        qo = qo | (q<<(GMP_LIMB_BITS - 1 - s)<<1);
        qp[j - 1] = qo;
        qo = q>>s;
        umul_ppmm(h, dummy, q, d);
        
        ASSERT(dummy == h1);
    }
    
    qp[n - 1] = qo;
    
    return h+c;
}
__GMP_DECLSPEC
mp_limb_t mpn_rsh_divrem_hensel_qr_1(mp_ptr qp, mp_srcptr xp,
                                     mp_size_t n, mp_limb_t d, int s, mp_limb_t cin)
{
    ASSERT(n > 0);
    ASSERT(s >= 0);
    ASSERT_MPN(xp, n);
    ASSERT(MPN_SAME_OR_SEPARATE_P(qp, xp, n));
    ASSERT(d%2 == 1);
    
    if (BELOW_THRESHOLD(n, RSH_DIVREM_HENSEL_QR_1_THRESHOLD))
        return mpn_rsh_divrem_hensel_qr_1_1(qp, xp, n, d, s, cin);
    
    return mpn_rsh_divrem_hensel_qr_1_2(qp, xp, n, d, s, cin);
}
/***************************************************************/
__GMP_DECLSPEC
mp_limb_t
mpn_divrem_1 (mp_ptr qp, mp_size_t qxn,
	      mp_srcptr up, mp_size_t un, mp_limb_t d)
{
  mp_size_t  n;
  mp_size_t  i;
  mp_limb_t  n1, n0;
  mp_limb_t  r = 0;

  ASSERT (qxn >= 0);
  ASSERT (un >= 0);
  ASSERT (d != 0);
  /* FIXME: What's the correct overlap rule when qxn!=0? */
  ASSERT (MPN_SAME_OR_SEPARATE_P (qp+qxn, up, un));

  n = un + qxn;
  if (n == 0)
    return 0;

  d <<= GMP_NAIL_BITS;

  if(qxn==0)
    {
     if(d<=GMP_LIMB_HIGHBIT/2+1 && ABOVE_THRESHOLD(un,DIVREM_EUCLID_HENSEL_THRESHOLD))
       {r=mpn_divrem_euclidean_r_1(up,un,d);
        count_trailing_zeros(i,d);
        mpn_rsh_divrem_hensel_qr_1(qp,up,un,d>>i,i,r);
        return r;}
  #if HAVE_NATIVE_mpn_divrem_euclidean_qr_1
     return mpn_divrem_euclidean_qr_1(qp,0,up,un,d);
  #endif
    }
  qp += (n - 1);   /* Make qp point at most significant quotient limb */

  if ((d & GMP_LIMB_HIGHBIT) != 0)
    {
      if (un != 0)
	{
	  /* High quotient limb is 0 or 1, skip a divide step. */
	  mp_limb_t q;
	  r = up[un - 1] << GMP_NAIL_BITS;
	  q = (r >= d);
	  *qp-- = q;
	  r -= (d & -q);
	  r >>= GMP_NAIL_BITS;
	  n--;
	  un--;
	}

      if (BELOW_THRESHOLD (n, DIVREM_1_NORM_THRESHOLD))
	{
	plain:
	  for (i = un - 1; i >= 0; i--)
	    {
	      n0 = up[i] << GMP_NAIL_BITS;
	      udiv_qrnnd (*qp, r, r, n0, d);
	      r >>= GMP_NAIL_BITS;
	      qp--;
	    }
	  for (i = qxn - 1; i >= 0; i--)
	    {
	      udiv_qrnnd (*qp, r, r, CNST_LIMB(0), d);
	      r >>= GMP_NAIL_BITS;
	      qp--;
	    }
	  return r;
	}
      else
	{
	  /* Multiply-by-inverse, divisor already normalized. */
	  mp_limb_t dinv;
	  invert_limb (dinv, d);

	  for (i = un - 1; i >= 0; i--)
	    {
	      n0 = up[i] << GMP_NAIL_BITS;
	      udiv_qrnnd_preinv (*qp, r, r, n0, d, dinv);
	      r >>= GMP_NAIL_BITS;
	      qp--;
	    }
	  for (i = qxn - 1; i >= 0; i--)
	    {
	      udiv_qrnnd_preinv (*qp, r, r, CNST_LIMB(0), d, dinv);
	      r >>= GMP_NAIL_BITS;
	      qp--;
	    }
	  return r;
	}
    }
  else
    {
      /* Most significant bit of divisor == 0.  */
      int norm;

      /* Skip a division if high < divisor (high quotient 0).  Testing here
	 before normalizing will still skip as often as possible.  */
      if (un != 0)
	{
	  n1 = up[un - 1] << GMP_NAIL_BITS;
	  if (n1 < d)
	    {
	      r = n1 >> GMP_NAIL_BITS;
	      *qp-- = 0;
	      n--;
	      if (n == 0)
		return r;
	      un--;
	    }
	}

      if (! UDIV_NEEDS_NORMALIZATION
	  && BELOW_THRESHOLD (n, DIVREM_1_UNNORM_THRESHOLD))
	goto plain;

      count_leading_zeros (norm, d);
      d <<= norm;
      r <<= norm;

      if (UDIV_NEEDS_NORMALIZATION
	  && BELOW_THRESHOLD (n, DIVREM_1_UNNORM_THRESHOLD))
	{
	  if (un != 0)
	    {
	      n1 = up[un - 1] << GMP_NAIL_BITS;
	      r |= (n1 >> (GMP_LIMB_BITS - norm));
	      for (i = un - 2; i >= 0; i--)
		{
		  n0 = up[i] << GMP_NAIL_BITS;
		  udiv_qrnnd (*qp, r, r,
			      (n1 << norm) | (n0 >> (GMP_NUMB_BITS - norm)),
			      d);
		  r >>= GMP_NAIL_BITS;
		  qp--;
		  n1 = n0;
		}
	      udiv_qrnnd (*qp, r, r, n1 << norm, d);
	      r >>= GMP_NAIL_BITS;
	      qp--;
	    }
	  for (i = qxn - 1; i >= 0; i--)
	    {
	      udiv_qrnnd (*qp, r, r, CNST_LIMB(0), d);
	      r >>= GMP_NAIL_BITS;
	      qp--;
	    }
	  return r >> norm;
	}
      else
	{
	  mp_limb_t  dinv;
	  invert_limb (dinv, d);
	  if (un != 0)
	    {
	      n1 = up[un - 1] << GMP_NAIL_BITS;
	      r |= (n1 >> (GMP_LIMB_BITS - norm));
	      for (i = un - 2; i >= 0; i--)
		{
		  n0 = up[i] << GMP_NAIL_BITS;
		  udiv_qrnnd_preinv (*qp, r, r, 
				     ((n1 << norm) | (n0 >> (GMP_NUMB_BITS - norm))),
				     d, dinv);
		  r >>= GMP_NAIL_BITS;
		  qp--;
		  n1 = n0;
		}
	      udiv_qrnnd_preinv (*qp, r, r, n1 << norm, d, dinv);
	      r >>= GMP_NAIL_BITS;
	      qp--;
	    }
	  for (i = qxn - 1; i >= 0; i--)
	    {
	      udiv_qrnnd_preinv (*qp, r, r, CNST_LIMB(0), d, dinv);
	      r >>= GMP_NAIL_BITS;
	      qp--;
	    }
	  return r >> norm;
	}
    }
}

#undef STORE_QUOTIENT
#undef UDIV_METHOD