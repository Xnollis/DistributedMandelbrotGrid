//#include "cuda_runtime.h"
#include<stdio.h>
#include<stdlib.h>
#include<memory.h>
#include "mpir_inter_decl.h"
/////////////////////////////////////////
const static int __gmp_default_fp_limb_precision=2;
/////////////////////////////////////////
void mpf_init (mpf_ptr r)
{
  mp_size_t prec = __gmp_default_fp_limb_precision;
  r->_mp_size = 0;
  r->_mp_exp = 0;
  r->_mp_prec = prec;
  memset(r->_mp_d, 0, sizeof(r->_mp_d));
}
void mpf_clear (mpf_ptr m)
{
	//Do nothing here
}
void mpf_set_d (mpf_ptr r, double d)
{
  int negative;
/* check if is a valid double value
  DOUBLE_NAN_INF_ACTION (d,
                         __gmp_invalid_operation (),
                         __gmp_invalid_operation ());
*/
  if (UNLIKELY (d == 0))
    {
      SIZ(r) = 0;
      EXP(r) = 0;
      return;
    }
  negative = d < 0;
  d = ABS (d);

  SIZ(r) = negative ? -LIMBS_PER_DOUBLE : LIMBS_PER_DOUBLE;
  EXP(r) = __gmp_extract_double (PTR(r), d);
}
double mpf_get_d(mpf_srcptr src)
{
	mp_size_t  size, abs_size;
	mp_exp_t       exp;

	size = SIZ(src);
	if (UNLIKELY(size == 0))
		return 0.0;

	abs_size = ABS(size);
	exp = (EXP(src) - abs_size) * GMP_NUMB_BITS;
	return mpn_get_d(PTR(src), abs_size, size, exp);//32-bit, 参数总大小应该是10h字节，但是现在由于有64位的参数，总大小变为18h了 《---由于注释这里有中文，会导致无法编译通过
}
double mpf_get_d_2exp (mp_exp_t *exp2, mpf_srcptr src)
{
  mp_size_t size, abs_size;
  mp_srcptr ptr;
  int cnt;
  mp_exp_t exp;

  size = SIZ(src);
  if (UNLIKELY (size == 0))
    {
      *exp2 = 0;
      return 0.0;
    }

  ptr = PTR(src);
  abs_size = ABS (size);
  count_leading_zeros (cnt, ptr[abs_size - 1]);
  cnt -= GMP_NAIL_BITS;

  exp = EXP(src) * GMP_NUMB_BITS - cnt;
  *exp2 = exp;

  return mpn_get_d (ptr, abs_size, size,
                    (long) - (abs_size * GMP_NUMB_BITS - cnt));
}
void mpf_mul_ui(mpf_ptr r, mpf_srcptr u, mpir_ui v)
{
	mp_srcptr up;
	mp_size_t usize;
	mp_size_t size;
	mp_size_t prec, excess;
	mp_limb_t cy_limb, vl, cbit, cin;
	mp_ptr rp;

	usize = u->_mp_size;
	if (UNLIKELY (v == 0) || UNLIKELY (usize == 0))
	{
		r->_mp_size = 0;
		r->_mp_exp = 0;
		return;
	}

#if BITS_PER_UI > GMP_NUMB_BITS  /* avoid warnings about shift amount */
	if (v > GMP_NUMB_MAX)
	{
		mpf_t     vf;
		mp_limb_t vp[2];
		vp[0] = v & GMP_NUMB_MASK;
		vp[1] = v >> GMP_NUMB_BITS;
		PTR(vf) = vp;
		SIZ(vf) = 2;
		ASSERT_CODE (PREC(vf) = 2);
		EXP(vf) = 2;
		mpf_mul (r, u, vf);
		return;
	}
#endif

	size = ABS (usize);
	prec = r->_mp_prec;
	up = u->_mp_d;
	vl = v;
	excess = size - prec;
	cin = 0;

	if (excess > 0)
	{
		/* up is bigger than desired rp, shorten it to prec limbs and
		determine a carry-in */

		mp_limb_t  vl_shifted = vl << GMP_NAIL_BITS;
		mp_limb_t  hi, lo, next_lo, sum;
		mp_size_t  i;

		/* high limb of top product */
		i = excess - 1;
		umul_ppmm(cin, lo, up[i], vl_shifted);

		/* and carry bit out of products below that, if any */
		for (;;)
		{
			i--;
			if (i < 0)
				break;

			umul_ppmm(hi, next_lo, up[i], vl_shifted);
			lo >>= GMP_NAIL_BITS;
			ADDC_LIMB(cbit, sum, hi, lo);
			cin += cbit;
			lo = next_lo;

			/* Continue only if the sum is GMP_NUMB_MAX.  GMP_NUMB_MAX is the
			only value a carry from below can propagate across.  If we've
			just seen the carry out (ie. cbit!=0) then sum!=GMP_NUMB_MAX,
			so this test stops us for that case too.  */
			if (LIKELY(sum != GMP_NUMB_MAX))
				break;
		}

		up += excess;
		size = prec;
	}

	rp = r->_mp_d;
#if HAVE_NATIVE_mpn_mul_1c
	cy_limb = mpn_mul_1c(rp, up, size, vl, cin);
#else
	cy_limb = mpn_mul_1(rp, up, size, vl);
	__GMPN_ADD_1(cbit, rp, rp, size, cin);
	cy_limb += cbit;
#endif
	rp[size] = cy_limb;
	cy_limb = cy_limb != 0;
	r->_mp_exp = u->_mp_exp + cy_limb;
	size += cy_limb;
	r->_mp_size = usize >= 0 ? size : -size;
}
void mpf_div_ui (mpf_ptr r, mpf_srcptr u, mpir_ui v)
{
    mp_srcptr up;
    mp_ptr rp, tp, rtp;
    mp_size_t usize;
    mp_size_t rsize, tsize;
    mp_size_t sign_quotient;
    mp_size_t prec;
    mp_limb_t q_limb;
    mp_exp_t rexp;
    TMP_DECL;
    
#if BITS_PER_UI > GMP_NUMB_BITS  /* avoid warnings about shift amount */
    if (v > GMP_NUMB_MAX)
    {
        mpf_t vf;
        mp_limb_t vl[2];
        SIZ(vf) = 2;
        EXP(vf) = 2;
        PTR(vf) = vl;
        vl[0] = v & GMP_NUMB_MASK;
        vl[1] = v >> GMP_NUMB_BITS;
        mpf_div (r, u, vf);
        return;
    }
#endif
    
    usize = u->_mp_size;
    sign_quotient = usize;
    usize = ABS (usize);
    prec = r->_mp_prec;
    
    if (v == 0)
    {
        DIVIDE_BY_ZERO;
        r->_mp_size=0;
        r->_mp_exp=0;
        return;
    }
    
    if (usize == 0)
    {
        r->_mp_size = 0;
        r->_mp_exp = 0;
        return;
    }
    
    TMP_MARK;
    
    rp = r->_mp_d;
    up = u->_mp_d;
    
    tsize = 1 + prec;
    tp = (mp_ptr) TMP_ALLOC ((tsize + 1) * BYTES_PER_MP_LIMB);
    
    if (usize > tsize)
    {
        up += usize - tsize;
        usize = tsize;
        rtp = tp;
    }
    else
    {
        MPN_ZERO (tp, tsize - usize);
        rtp = tp + (tsize - usize);
    }
    
    /* Move the dividend to the remainder.  */
    MPN_COPY (rtp, up, usize);
    
    mpn_divmod_1 (rp, tp, tsize, (mp_limb_t) v);
    q_limb = rp[tsize - 1];
    
    rsize = tsize - (q_limb == 0);
    rexp = u->_mp_exp - (q_limb == 0);
    r->_mp_size = sign_quotient >= 0 ? rsize : -rsize;
    r->_mp_exp = rexp;
    TMP_FREE;
}

void mpf_add (mpf_ptr r, mpf_srcptr u, mpf_srcptr v)
{
    mp_srcptr up, vp;
    mp_ptr rp, tp;
    mp_size_t usize, vsize, rsize;
    mp_size_t prec;
    mp_exp_t uexp;
    mp_size_t ediff;
    mp_limb_t cy;
    int negate;
    TMP_DECL;
    
    usize = u->_mp_size;
    vsize = v->_mp_size;
    
    /* Handle special cases that don't work in generic code below.  */
    if (usize == 0)
    {
    set_r_v_maybe:
        if (r != v)
            mpf_set (r, v);
        return;
    }
    if (vsize == 0)
    {
        v = u;
        goto set_r_v_maybe;
    }
    
    /* If signs of U and V are different, perform subtraction.  */
    if ((usize ^ vsize) < 0)
    {
        __mpf_struct v_negated;
        v_negated._mp_size = -vsize;
        v_negated._mp_exp = v->_mp_exp;
//        v_negated._mp_d = v->_mp_d;
        memcpy(v_negated._mp_d,v->_mp_d,v->_mp_size*sizeof(mp_limb_t));// assign ptr CHANGE to memcpy
        mpf_sub (r, u, &v_negated);
        return;
    }
    
    TMP_MARK;
    
    /* Signs are now known to be the same.  */
    negate = usize < 0;
    
    /* Make U be the operand with the largest exponent.  */
    if (u->_mp_exp < v->_mp_exp)
    {
        mpf_srcptr t;
        t = u; u = v; v = t;
        usize = u->_mp_size;
        vsize = v->_mp_size;
    }
    
    usize = ABS (usize);
    vsize = ABS (vsize);
    up = u->_mp_d;
    vp = v->_mp_d;
    rp = r->_mp_d;
    prec = r->_mp_prec;
    uexp = u->_mp_exp;
    ediff = u->_mp_exp - v->_mp_exp;
    
    /* If U extends beyond PREC, ignore the part that does.  */
    if (usize > prec)
    {
        up += usize - prec;
        usize = prec;
    }
    
    /* If V extends beyond PREC, ignore the part that does.
     Note that this may make vsize negative.  */
    if (vsize + ediff > prec)
    {
        vp += vsize + ediff - prec;
        vsize = prec - ediff;
    }
    
#if 0
    /* Locate the least significant non-zero limb in (the needed parts
     of) U and V, to simplify the code below.  */
    while (up[0] == 0)
        up++, usize--;
    while (vp[0] == 0)
        vp++, vsize--;
#endif
    
    /* Allocate temp space for the result.  Allocate
     just vsize + ediff later???  */
    tp = (mp_ptr) TMP_ALLOC (prec * BYTES_PER_MP_LIMB);
    
    if (ediff >= prec)
    {
        /* V completely cancelled.  */
        if (rp != up)
            MPN_COPY_INCR (rp, up, usize);
        rsize = usize;
    }
    else
    {
        /* uuuu     |  uuuu     |  uuuu     |  uuuu     |  uuuu    */
        /* vvvvvvv  |  vv       |    vvvvv  |    v      |       vv */
        
        if (usize > ediff)
        {
            /* U and V partially overlaps.  */
            if (vsize + ediff <= usize)
            {
                /* uuuu     */
                /*   v      */
                mp_size_t size;
                size = usize - ediff - vsize;
                MPN_COPY (tp, up, size);
                cy = mpn_add (tp + size, up + size, usize - size, vp, vsize);
                rsize = usize;
            }
            else
            {
                /* uuuu     */
                /*   vvvvv  */
                mp_size_t size;
                size = vsize + ediff - usize;
                MPN_COPY (tp, vp, size);
                cy = mpn_add (tp + size, up, usize, vp + size, usize - ediff);
                rsize = vsize + ediff;
            }
        }
        else
        {
            /* uuuu     */
            /*      vv  */
            mp_size_t size;
            size = vsize + ediff - usize;
            MPN_COPY (tp, vp, vsize);
            MPN_ZERO (tp + vsize, ediff - usize);
            MPN_COPY (tp + size, up, usize);
            cy = 0;
            rsize = size + usize;
        }
        
        MPN_COPY (rp, tp, rsize);
        rp[rsize] = cy;
        rsize += cy;
        uexp += cy;
    }
    
    r->_mp_size = negate ? -rsize : rsize;
    r->_mp_exp = uexp;
    TMP_FREE;
}

void mpf_sub (mpf_ptr r, mpf_srcptr u, mpf_srcptr v)
{
    mp_srcptr up, vp;
    mp_ptr rp, tp;
    mp_size_t usize, vsize, rsize;
    mp_size_t prec;
    mp_exp_t exp;
    mp_size_t ediff;
    int negate;
    TMP_DECL;
    
    usize = u->_mp_size;
    vsize = v->_mp_size;
    
    /* Handle special cases that don't work in generic code below.  */
    if (usize == 0)
    {
        mpf_neg (r, v);
        return;
    }
    if (vsize == 0)
    {
        if (r != u)
            mpf_set (r, u);
        return;
    }
    
    /* If signs of U and V are different, perform addition.  */
    if ((usize ^ vsize) < 0)
    {
        __mpf_struct v_negated;
        v_negated._mp_size = -vsize;
        v_negated._mp_exp = v->_mp_exp;
        //v_negated._mp_d = v->_mp_d;
        memcpy(v_negated._mp_d,v->_mp_d,v->_mp_size*sizeof(mp_limb_t));// assign ptr CHANGE to memcpy
        mpf_add (r, u, &v_negated);
        return;
    }
    
    TMP_MARK;
    
    /* Signs are now known to be the same.  */
    negate = usize < 0;
    
    /* Make U be the operand with the largest exponent.  */
    if (u->_mp_exp < v->_mp_exp)
    {
        mpf_srcptr t;
        t = u; u = v; v = t;
        negate ^= 1;
        usize = u->_mp_size;
        vsize = v->_mp_size;
    }
    
    usize = ABS (usize);
    vsize = ABS (vsize);
    up = u->_mp_d;
    vp = v->_mp_d;
    rp = r->_mp_d;
    prec = r->_mp_prec + 1;
    exp = u->_mp_exp;
    ediff = u->_mp_exp - v->_mp_exp;
    
    /* If ediff is 0 or 1, we might have a situation where the operands are
     extremely close.  We need to scan the operands from the most significant
     end ignore the initial parts that are equal.  */
    if (ediff <= 1)
    {
        if (ediff == 0)
        {
            /* Skip leading limbs in U and V that are equal.  */
            if (up[usize - 1] == vp[vsize - 1])
            {
                /* This loop normally exits immediately.  Optimize for that.  */
                do
                {
                    usize--;
                    vsize--;
                    exp--;
                    
                    if (usize == 0)
                    {
                        /* u cancels high limbs of v, result is rest of v */
                        negate ^= 1;
                    cancellation:
                        /* strip high zeros before truncating to prec */
                        while (vsize != 0 && vp[vsize - 1] == 0)
                        {
                            vsize--;
                            exp--;
                        }
                        if (vsize > prec)
                        {
                            vp += vsize - prec;
                            vsize = prec;
                        }
                        MPN_COPY_INCR (rp, vp, vsize);
                        rsize = vsize;
                        goto done;
                    }
                    if (vsize == 0)
                    {
                        vp = up;
                        vsize = usize;
                        goto cancellation;
                    }
                }
                while (up[usize - 1] == vp[vsize - 1]);
            }
            
            if (up[usize - 1] < vp[vsize - 1])
            {
                /* For simplicity, swap U and V.  Note that since the loop above
                 wouldn't have exited unless up[usize - 1] and vp[vsize - 1]
                 were non-equal, this if-statement catches all cases where U
                 is smaller than V.  */
                MPN_SRCPTR_SWAP (up,usize, vp,vsize);
                negate ^= 1;
                /* negating ediff not necessary since it is 0.  */
            }
            
            /* Check for
             x+1 00000000 ...
             x  ffffffff ... */
            if (up[usize - 1] != vp[vsize - 1] + 1)
                goto general_case;
            usize--;
            vsize--;
            exp--;
        }
        else /* ediff == 1 */
        {
            /* Check for
             1 00000000 ...
             0 ffffffff ... */
            
            if (up[usize - 1] != 1 || vp[vsize - 1] != GMP_NUMB_MAX
                || (usize >= 2 && up[usize - 2] != 0))
                goto general_case;
            
            usize--;
            exp--;
        }
        
        /* Skip sequences of 00000000/ffffffff */
        while (vsize != 0 && usize != 0 && up[usize - 1] == 0
               && vp[vsize - 1] == GMP_NUMB_MAX)
        {
            usize--;
            vsize--;
            exp--;
        }
        
        if (usize == 0)
        {
            while (vsize != 0 && vp[vsize - 1] == GMP_NUMB_MAX)
            {
                vsize--;
                exp--;
            }
        }
        
        if (usize > prec - 1)
        {
            up += usize - (prec - 1);
            usize = prec - 1;
        }
        if (vsize > prec - 1)
        {
            vp += vsize - (prec - 1);
            vsize = prec - 1;
        }
        
        tp = (mp_ptr) TMP_ALLOC (prec * BYTES_PER_MP_LIMB);
        {
            mp_limb_t cy_limb;
            if (vsize == 0)
            {
                mp_size_t size, i;
                size = usize;
                for (i = 0; i < size; i++)
                    tp[i] = up[i];
                tp[size] = 1;
                rsize = size + 1;
                exp++;
                goto normalize;
            }
            if (usize == 0)
            {
                mp_size_t size, i;
                size = vsize;
                for (i = 0; i < size; i++)
                    tp[i] = ~vp[i] & GMP_NUMB_MASK;
                cy_limb = 1 - mpn_add_1 (tp, tp, vsize, (mp_limb_t) 1);
                rsize = vsize;
                if (cy_limb == 0)
                {
                    tp[rsize] = 1;
                    rsize++;
                    exp++;
                }
                goto normalize;
            }
            if (usize >= vsize)
            {
                /* uuuu     */
                /* vv       */
                mp_size_t size;
                size = usize - vsize;
                MPN_COPY (tp, up, size);
                cy_limb = mpn_sub_n (tp + size, up + size, vp, vsize);
                rsize = usize;
            }
            else /* (usize < vsize) */
            {
                /* uuuu     */
                /* vvvvvvv  */
                mp_size_t size, i;
                size = vsize - usize;
                for (i = 0; i < size; i++)
                    tp[i] = ~vp[i] & GMP_NUMB_MASK;
                cy_limb = mpn_sub_n (tp + size, up, vp + size, usize);
                cy_limb+= mpn_sub_1 (tp + size, tp + size, usize, (mp_limb_t) 1);
                cy_limb-= mpn_add_1 (tp, tp, vsize, (mp_limb_t) 1);
                rsize = vsize;
            }
            if (cy_limb == 0)
            {
                tp[rsize] = 1;
                rsize++;
                exp++;
            }
            goto normalize;
        }
    }
    
general_case:
    /* If U extends beyond PREC, ignore the part that does.  */
    if (usize > prec)
    {
        up += usize - prec;
        usize = prec;
    }
    
    /* If V extends beyond PREC, ignore the part that does.
     Note that this may make vsize negative.  */
    if (vsize + ediff > prec)
    {
        vp += vsize + ediff - prec;
        vsize = prec - ediff;
    }
    
    /* Allocate temp space for the result.  Allocate
     just vsize + ediff later???  */
    tp = (mp_ptr) TMP_ALLOC (prec * BYTES_PER_MP_LIMB);
    
    if (ediff >= prec)
    {
        /* V completely cancelled.  */
        if (tp != up)
            MPN_COPY (rp, up, usize);
        rsize = usize;
    }
    else
    {
        /* Locate the least significant non-zero limb in (the needed
         parts of) U and V, to simplify the code below.  */
        for (;;)
        {
            if (vsize == 0)
            {
                MPN_COPY (rp, up, usize);
                rsize = usize;
                goto done;
            }
            if (vp[0] != 0)
                break;
            vp++, vsize--;
        }
        for (;;)
        {
            if (usize == 0)
            {
                MPN_COPY (rp, vp, vsize);
                rsize = vsize;
                negate ^= 1;
                goto done;
            }
            if (up[0] != 0)
                break;
            up++, usize--;
        }
        
        /* uuuu     |  uuuu     |  uuuu     |  uuuu     |  uuuu    */
        /* vvvvvvv  |  vv       |    vvvvv  |    v      |       vv */
        
        if (usize > ediff)
        {
            /* U and V partially overlaps.  */
            if (ediff == 0)
            {
                /* Have to compare the leading limbs of u and v
                 to determine whether to compute u - v or v - u.  */
                if (usize >= vsize)
                {
                    /* uuuu     */
                    /* vv       */
                    mp_size_t size;
                    size = usize - vsize;
                    MPN_COPY (tp, up, size);
                    mpn_sub_n (tp + size, up + size, vp, vsize);
                    rsize = usize;
                }
                else /* (usize < vsize) */
                {
                    /* uuuu     */
                    /* vvvvvvv  */
                    mp_size_t size, i;
                    size = vsize - usize;
                    tp[0] = -vp[0] & GMP_NUMB_MASK;
                    for (i = 1; i < size; i++)
                        tp[i] = ~vp[i] & GMP_NUMB_MASK;
                    mpn_sub_n (tp + size, up, vp + size, usize);
                    mpn_sub_1 (tp + size, tp + size, usize, (mp_limb_t) 1);
                    rsize = vsize;
                }
            }
            else
            {
                if (vsize + ediff <= usize)
                {
                    /* uuuu     */
                    /*   v      */
                    mp_size_t size;
                    size = usize - ediff - vsize;
                    MPN_COPY (tp, up, size);
                    mpn_sub (tp + size, up + size, usize - size, vp, vsize);
                    rsize = usize;
                }
                else
                {
                    /* uuuu     */
                    /*   vvvvv  */
                    mp_size_t size, i;
                    size = vsize + ediff - usize;
                    tp[0] = -vp[0] & GMP_NUMB_MASK;
                    for (i = 1; i < size; i++)
                        tp[i] = ~vp[i] & GMP_NUMB_MASK;
                    mpn_sub (tp + size, up, usize, vp + size, usize - ediff);
                    mpn_sub_1 (tp + size, tp + size, usize, (mp_limb_t) 1);
                    rsize = vsize + ediff;
                }
            }
        }
        else
        {
            /* uuuu     */
            /*      vv  */
            mp_size_t size, i;
            size = vsize + ediff - usize;
            tp[0] = -vp[0] & GMP_NUMB_MASK;
            for (i = 1; i < vsize; i++)
                tp[i] = ~vp[i] & GMP_NUMB_MASK;
            for (i = vsize; i < size; i++)
                tp[i] = GMP_NUMB_MAX;
            mpn_sub_1 (tp + size, up, usize, (mp_limb_t) 1);
            rsize = size + usize;
        }
        
    normalize:
        /* Full normalize.  Optimize later.  */
        while (rsize != 0 && tp[rsize - 1] == 0)
        {
            rsize--;
            exp--;
        }
        MPN_COPY (rp, tp, rsize);
    }
    
done:
    r->_mp_size = negate ? -rsize : rsize;
    if (rsize == 0)
        exp = 0;
    r->_mp_exp = exp;
    TMP_FREE;
}
void mpf_neg (mpf_ptr r, mpf_srcptr u)
{
    mp_size_t size;
    
    size = -u->_mp_size;
    if (r != u)
    {
        mp_size_t prec;
        mp_size_t asize;
        mp_ptr rp;
        mp_srcptr up;
        
        prec = r->_mp_prec + 1;    /* lie not to lose precision in assignment */
        asize = ABS (size);
        rp = r->_mp_d;
        up = u->_mp_d;
        
        if (asize > prec)
        {
            up += asize - prec;
            asize = prec;
        }
        
        MPN_COPY (rp, up, asize);
        r->_mp_exp = u->_mp_exp;
        size = size >= 0 ? asize : -asize;
    }
    r->_mp_size = size;
}
void mpf_set (mpf_ptr r, mpf_srcptr u)
{
    mp_ptr rp;
    mp_srcptr up;
    mp_size_t size, asize;
    mp_size_t prec;
    
    prec = r->_mp_prec + 1;        /* lie not to lose precision in assignment */
    size = u->_mp_size;
    asize = ABS (size);
    rp = r->_mp_d;
    up = u->_mp_d;
    
    if (asize > prec)
    {
        up += asize - prec;
        asize = prec;
    }
    
    r->_mp_exp = u->_mp_exp;
    r->_mp_size = size >= 0 ? asize : -asize;
    MPN_COPY_INCR (rp, up, asize);
}

void *TMP_ALLOC_FUNC(size_t n,char pLocalBuf[TMP_ALLOC_LOCAL_BUF_SIZE],size_t *pCnt)
{
    if(!pCnt||!pLocalBuf||n<=0)return NULL;
    size_t nRemain=TMP_ALLOC_LOCAL_BUF_SIZE-*pCnt;
    if(n>nRemain)return NULL;
    char *p=pLocalBuf+*pCnt;
    (*pCnt)+=n;
    return p;
}
/////////////////////////////////////////
#if 0
int IsCUDA_Supported(int bPrintInfoToConsole)
{
	cudaError_t err;
	int deviceCount;
	int device = 0,iGoodDeviceCount=0;

	err = cudaGetDeviceCount(&deviceCount);
	if (err != cudaSuccess)
		return 0;
	if (deviceCount <= 0)
		return 0;
	printf("找到%d个支持CUDA的显卡。以下是其详细信息：\n",deviceCount);
	struct cudaDeviceProp deviceProp;
	for (device = 0; device < deviceCount; ++device)
	{
		err = cudaGetDeviceProperties(&deviceProp, device);
		if (err == cudaSuccess)
		{
			if (bPrintInfoToConsole)
			{
				printf("显卡设备 %d 详情:\n%s 的计算能力(Compute Capability)为： %d.%d\n", device, deviceProp.name, deviceProp.major, deviceProp.minor);
				//枚举详细信息
				//printf("Identify: %s\n", deviceProp.name);
				//printf("Host Memory: %d\n", deviceProp.canMapHostMemory);
				printf("其他信息：\n显存大小： %zu 字节\n", deviceProp.totalGlobalMem);
				printf("核心处理器个数: %d\n", deviceProp.multiProcessorCount);
				printf("时钟频率 %d khz\n", deviceProp.clockRate);
#ifdef DEBUG
				printf("Compute Mode: %d\n", deviceProp.computeMode);
				printf("Device Overlap: %d\n", deviceProp.deviceOverlap);
				printf("Integrated: %d\n", deviceProp.integrated);
				printf("Kernel Exec Timeout Enabled: %d\n", deviceProp.kernelExecTimeoutEnabled);
				printf("Max Grid Size: %d * %d * %d\n", deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]);
				printf("Max Threads Dim: %d * %d * %d\n", deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]);
				printf("Max Threads per Block: %d\n", deviceProp.maxThreadsPerBlock);
				printf("Maximum Pitch: %d bytes\n", deviceProp.memPitch);
				printf("32bit Registers Availble per Block: %d\n", deviceProp.regsPerBlock);
				printf("Shared Memory Available per Block: %d bytes\n", deviceProp.sharedMemPerBlock);
				printf("Alignment Requirement for Textures: %d\n", deviceProp.textureAlignment);
				printf("Constant Memory Available: %d bytes\n", deviceProp.totalConstMem);
				printf("Warp Size: %d threads\n", deviceProp.warpSize);
#endif
				printf("=========================================\n");
			}
			if (deviceProp.major >= 2)//, deviceProp.minor);
				++iGoodDeviceCount;
		}
	}
	return (iGoodDeviceCount>0);
}
#endif
