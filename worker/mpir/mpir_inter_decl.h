#ifndef __MPIR_INTER_DECL_H__
#define __MPIR_INTER_DECL_H__
#include"mpir.h"
#include "longlong.h"
mp_limb_t mpn_add_n (mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n);
mp_limb_t mpn_sub_n (mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n);
mp_limb_t mpn_add (mp_ptr __gmp_wp, mp_srcptr __gmp_xp, mp_size_t __gmp_xsize, mp_srcptr __gmp_yp, mp_size_t __gmp_ysize);
mp_limb_t mpn_add_1 (mp_ptr __gmp_dst, mp_srcptr __gmp_src, mp_size_t __gmp_size, mp_limb_t __gmp_n);
int mpn_cmp (mp_srcptr __gmp_xp, mp_srcptr __gmp_yp, mp_size_t __gmp_size);
mp_limb_t mpn_sub (mp_ptr __gmp_wp, mp_srcptr __gmp_xp, mp_size_t __gmp_xsize, mp_srcptr __gmp_yp, mp_size_t __gmp_ysize);
mp_limb_t mpn_sub_1 (mp_ptr __gmp_dst, mp_srcptr __gmp_src, mp_size_t __gmp_size, mp_limb_t __gmp_n);
mp_limb_t mpn_invert_limb (mp_limb_t);
#include "gmpn_add.h"
int __gmp_extract_double(mp_ptr rp, double d);
double mpn_get_d(mp_srcptr ptr, mp_size_t size, mp_size_t sign, long exp);
mp_limb_t mpn_mul_1 (mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl);
mp_limb_t mpn_divrem_1 (mp_ptr qp, mp_size_t qxn,mp_srcptr up, mp_size_t un, mp_limb_t d);
mp_limb_t mpn_mod_1 (mp_srcptr up, mp_size_t un, mp_limb_t d);
void mpn_mod_1_1(mp_ptr rem, mp_srcptr xp, mp_size_t xn, mp_srcptr db);
void mpn_mod_1_2(mp_ptr rem, mp_srcptr xp, mp_size_t xn, mp_srcptr db);
void mpn_mod_1_3(mp_ptr rem, mp_srcptr xp, mp_size_t xn, mp_srcptr db);
void mpn_mul_n(mp_ptr p, mp_srcptr a, mp_srcptr b, mp_size_t n);
mp_limb_t mpn_mul(mp_ptr prodp, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn);
void mpn_sqr(mp_ptr p, mp_srcptr a, mp_size_t n);

void mpn_mul_basecase(mp_ptr rp,mp_srcptr up, mp_size_t un,mp_srcptr vp, mp_size_t vn);
void mpn_toom4_interpolate(mp_ptr rp, mp_size_t * rpn, mp_size_t sn,
	mp_ptr tp, mp_size_t s4, mp_size_t n4, mp_size_t n6, mp_limb_t r30);
void mpn_toom53_mul(mp_ptr rp, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn);
void mpn_toom8h_mul(mp_ptr pp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn);
void mpn_toom4_mul(mp_ptr rp, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn);
void mpn_toom3_mul(mp_ptr c, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_ptr t);
void mpn_toom3_mul_n(mp_ptr c, mp_srcptr a, mp_srcptr b, mp_size_t n, mp_ptr t);
void mpn_toom4_mul_n(mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n);
mp_limb_t mpn_lshift(mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt);
mp_limb_t mpn_rshift(mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt);
mp_limb_t mpn_addmul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl);
mp_limb_t mpn_addmul_2(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_srcptr vp);
mp_limb_t mpn_submul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl);
void mpn_kara_mul_n(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n, mp_ptr tp);
#define mpn_half(__xp,__n) mpn_rshift1((__xp),(__xp),(__n))
#define mpn_addlsh1_n(__xp,__yp,__zp,__n) mpn_addlsh_n((__xp),(__yp),(__zp),(__n),1)
#define mpn_rshift1(__xp,__yp,__n) mpn_rshift((__xp),(__yp),(__n),1)
#define mpn_lshift1(__xp,__yp,__n) mpn_lshift((__xp),(__yp),(__n),1)
#define mpn_double(__xp,__n) mpn_lshift1((__xp),(__xp),(__n))
#define mpn_not(__xp,__n) mpn_com_n((__xp),(__xp),(__n))

#define MPN_DECR_U(ptr, size, n)   mpn_decr_u (ptr, n)
#define MPN_INCR_U(ptr, size, n)   mpn_incr_u (ptr, n)
#define mpn_divexact_by3(dst,src,size)  mpn_divexact_by3c (dst, src, size, __GMP_CAST (mp_limb_t, 0))

int mpn_toom_eval_pm2(mp_ptr xp2, mp_ptr xm2, unsigned k, mp_srcptr xp, mp_size_t n, mp_size_t hn, mp_ptr tp);
int mpn_toom_eval_pm2exp(mp_ptr xp2, mp_ptr xm2, unsigned k,mp_srcptr xp, mp_size_t n, mp_size_t hn, unsigned shift,mp_ptr tp);
int mpn_toom_eval_pm2rexp(mp_ptr rp, mp_ptr rm, unsigned int q, mp_srcptr ap, mp_size_t n, mp_size_t t, unsigned int s, mp_ptr ws);
int mpn_toom_eval_dgr3_pm1(mp_ptr xp1, mp_ptr xm1, mp_srcptr xp, mp_size_t n, mp_size_t x3n, mp_ptr tp);
void mpn_toom_interpolate_16pts(mp_ptr pp, mp_ptr r1, mp_ptr r3, mp_ptr r5, mp_ptr r7, mp_size_t n, mp_size_t spt, int half, mp_ptr wsi);
void mpn_divexact_1(mp_ptr dst, mp_srcptr src, mp_size_t size, mp_limb_t divisor);
mp_limb_t mpn_divexact_by3c(mp_ptr qp, mp_srcptr xp, mp_size_t n, mp_limb_t ci);
mp_limb_t mpn_divexact_byfobm1(mp_ptr qp, mp_srcptr xp, mp_size_t n,mp_limb_t f, mp_limb_t Bm1of);

mp_limb_t mpn_divrem_euclidean_qr_1(mp_ptr qp, mp_size_t qxn, mp_srcptr xp, mp_size_t n, mp_limb_t d);
mp_limb_t mpn_divrem_euclidean_r_1(mp_srcptr xp, mp_size_t n, mp_limb_t d);
mp_limb_t mpn_rsh_divrem_hensel_qr_1(mp_ptr qp, mp_srcptr xp, mp_size_t n, mp_limb_t d, int s, mp_limb_t cin);
mp_limb_t mpn_divrem_1(mp_ptr qp, mp_size_t qxn,mp_srcptr up, mp_size_t un, mp_limb_t d);
mp_limb_t mpn_divrem_2(mp_ptr qp, mp_size_t qxn,mp_ptr np, mp_size_t nn,mp_srcptr dp);
mp_limb_t mpn_sb_div_qr(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv);
mp_limb_t mpn_dc_div_qr_n(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n,mp_limb_t dinv, mp_ptr tp);
mp_limb_t mpn_inv_div_qr(mp_ptr qp,mp_ptr np, mp_size_t nn,mp_srcptr dp, mp_size_t dn,mp_srcptr dinv);
mp_limb_t mpn_inv_div_qr_n(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t dn, mp_srcptr inv);
#define DC_DIVAPPR_Q_N_ITCH(n) ((n)*4 + 64)
#define DC_BDIV_Q_N_ITCH(n) ((n)/2 + 2)
#define DC_BDIV_QR_N_ITCH(n) (n)

#endif
