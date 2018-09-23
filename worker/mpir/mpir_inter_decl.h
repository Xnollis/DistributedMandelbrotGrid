#ifndef __MPIR_INTER_DECL_H__
#define __MPIR_INTER_DECL_H__
#include"mpir.h"
#include "longlong.h"
#ifdef MPIR_CUDA_ACC
extern const unsigned char __GMP_DECLSPEC_G_VALUE  modlimb_invert_tabled[128];
#ifdef COUNT_LEADING_ZEROS_NEED_CLZ_TAB
extern const unsigned char __GMP_DECLSPEC_G_VALUE __clz_tabd[129];
#define THE_CLZ_TAB __clz_tabd
#define THE_INVERT_TABLE modlimb_invert_tabled
#endif
#else
extern const unsigned char __GMP_DECLSPEC_G_VALUE  modlimb_invert_table[128];
#ifdef COUNT_LEADING_ZEROS_NEED_CLZ_TAB
extern const unsigned char __GMP_DECLSPEC_G_VALUE __clz_tab[129];
#define THE_CLZ_TAB __clz_tab
#define THE_INVERT_TABLE modlimb_invert_table
#endif
#endif // MPIR_CUDA_ACC
__GMP_DECLSPEC mp_limb_t mpn_add_n(mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n);
__GMP_DECLSPEC mp_limb_t mpn_sub_n(mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n);
__GMP_DECLSPEC mp_limb_t mpn_add(mp_ptr __gmp_wp, mp_srcptr __gmp_xp, mp_size_t __gmp_xsize, mp_srcptr __gmp_yp, mp_size_t __gmp_ysize);
__GMP_DECLSPEC mp_limb_t mpn_add_1(mp_ptr __gmp_dst, mp_srcptr __gmp_src, mp_size_t __gmp_size, mp_limb_t __gmp_n);
__GMP_DECLSPEC int mpn_cmp(mp_srcptr __gmp_xp, mp_srcptr __gmp_yp, mp_size_t __gmp_size);
__GMP_DECLSPEC mp_limb_t mpn_sub(mp_ptr __gmp_wp, mp_srcptr __gmp_xp, mp_size_t __gmp_xsize, mp_srcptr __gmp_yp, mp_size_t __gmp_ysize);
__GMP_DECLSPEC mp_limb_t mpn_sub_1(mp_ptr __gmp_dst, mp_srcptr __gmp_src, mp_size_t __gmp_size, mp_limb_t __gmp_n);
__GMP_DECLSPEC mp_limb_t mpn_invert_limb(mp_limb_t);
#include "gmpn_add.h"
__GMP_DECLSPEC int __gmp_extract_double(mp_ptr rp, double d);
__GMP_DECLSPEC double mpn_get_d(mp_srcptr ptr, mp_size_t size, mp_size_t sign, long exp);
__GMP_DECLSPEC mp_limb_t mpn_mul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl);
__GMP_DECLSPEC mp_limb_t mpn_divrem_1(mp_ptr qp, mp_size_t qxn, mp_srcptr up, mp_size_t un, mp_limb_t d);
__GMP_DECLSPEC mp_limb_t mpn_mod_1(mp_srcptr up, mp_size_t un, mp_limb_t d);
__GMP_DECLSPEC void mpn_mod_1_1(mp_ptr rem, mp_srcptr xp, mp_size_t xn, mp_srcptr db);
__GMP_DECLSPEC void mpn_mod_1_2(mp_ptr rem, mp_srcptr xp, mp_size_t xn, mp_srcptr db);
__GMP_DECLSPEC void mpn_mod_1_3(mp_ptr rem, mp_srcptr xp, mp_size_t xn, mp_srcptr db);
__GMP_DECLSPEC void mpn_mul_n(mp_ptr p, mp_srcptr a, mp_srcptr b, mp_size_t n);
__GMP_DECLSPEC mp_limb_t mpn_mul(mp_ptr prodp, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn);
__GMP_DECLSPEC void mpn_sqr(mp_ptr p, mp_srcptr a, mp_size_t n);

__GMP_DECLSPEC void mpn_mul_basecase(mp_ptr rp, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn);
__GMP_DECLSPEC void mpn_toom4_interpolate(mp_ptr rp, mp_size_t * rpn, mp_size_t sn, mp_ptr tp, mp_size_t s4, mp_size_t n4, mp_size_t n6, mp_limb_t r30);
__GMP_DECLSPEC void mpn_toom53_mul(mp_ptr rp, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn);
__GMP_DECLSPEC void mpn_toom8h_mul(mp_ptr pp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn);
__GMP_DECLSPEC void mpn_toom4_mul(mp_ptr rp, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn);
__GMP_DECLSPEC void mpn_toom3_mul(mp_ptr c, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_ptr t);
__GMP_DECLSPEC void mpn_toom3_mul_n(mp_ptr c, mp_srcptr a, mp_srcptr b, mp_size_t n, mp_ptr t);
__GMP_DECLSPEC void mpn_toom4_mul_n(mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n);
__GMP_DECLSPEC void mpn_toom42_mul(mp_ptr c, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_ptr t);
__GMP_DECLSPEC void mpn_toom32_mul(mp_ptr c, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_ptr t);
__GMP_DECLSPEC void mpn_toom3_sqr_n(mp_ptr c, mp_srcptr a, mp_size_t n, mp_ptr t);
__GMP_DECLSPEC void mpn_toom4_sqr_n(mp_ptr rp, mp_srcptr up, mp_size_t n);
__GMP_DECLSPEC void mpn_toom8_sqr_n(mp_ptr pp, mp_srcptr ap, mp_size_t an);
__GMP_DECLSPEC void mpn_toom3_interpolate(mp_ptr c, mp_ptr v1, mp_ptr v2, mp_ptr vm1, mp_ptr vinf, mp_size_t k, mp_size_t rr2, int sa, mp_limb_t vinf0, mp_ptr ws);
__GMP_DECLSPEC void tc4_addlsh1_unsigned(mp_ptr rp, mp_size_t * rn, mp_srcptr xp, mp_size_t xn);
__GMP_DECLSPEC void tc4_addmul_1(mp_ptr wp, mp_size_t * wn, mp_srcptr xp, mp_size_t xn, mp_limb_t y);
__GMP_DECLSPEC void tc4_add(mp_ptr rp, mp_size_t * rn, mp_srcptr r1, mp_size_t r1n, mp_srcptr r2, mp_size_t r2n);
__GMP_DECLSPEC void mpn_toom_couple_handling(mp_ptr pp, mp_size_t n, mp_ptr np, int nsign, mp_size_t off, int ps, int ns);
__GMP_DECLSPEC void mpn_kara_sqr_n(mp_ptr rp, mp_srcptr xp, mp_size_t n, mp_ptr tp);
__GMP_DECLSPEC void mpn_invert(mp_ptr xp, mp_srcptr ap, mp_size_t n);
__GMP_DECLSPEC mp_limb_t mpn_lshift(mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt);
__GMP_DECLSPEC mp_limb_t mpn_rshift(mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt);
__GMP_DECLSPEC mp_limb_t mpn_addmul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl);
__GMP_DECLSPEC mp_limb_t mpn_addmul_2(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_srcptr vp);
__GMP_DECLSPEC mp_limb_t mpn_submul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl);
__GMP_DECLSPEC void mpn_kara_mul_n(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n, mp_ptr tp);
__GMP_DECLSPEC void mpn_com_n(mp_ptr rp, mp_srcptr up, mp_size_t n);
#define mpn_half(__xp,__n) mpn_rshift1((__xp),(__xp),(__n))
#define mpn_addlsh1_n(__xp,__yp,__zp,__n) mpn_addlsh_n((__xp),(__yp),(__zp),(__n),1)
#define mpn_rshift1(__xp,__yp,__n) mpn_rshift((__xp),(__yp),(__n),1)
#define mpn_lshift1(__xp,__yp,__n) mpn_lshift((__xp),(__yp),(__n),1)
#define mpn_double(__xp,__n) mpn_lshift1((__xp),(__xp),(__n))
#define mpn_not(__xp,__n) mpn_com_n((__xp),(__xp),(__n))

#define MPN_DECR_U(ptr, size, n)   mpn_decr_u (ptr, n)
#define MPN_INCR_U(ptr, size, n)   mpn_incr_u (ptr, n)
#define mpn_divexact_by3(dst,src,size)  mpn_divexact_by3c (dst, src, size, __GMP_CAST (mp_limb_t, 0))

__GMP_DECLSPEC int	mpn_toom_eval_pm1(mp_ptr pp, mp_ptr mp, unsigned int k, mp_srcptr xp, mp_size_t n, mp_size_t m, mp_ptr tp);
__GMP_DECLSPEC int mpn_toom_eval_pm2(mp_ptr xp2, mp_ptr xm2, unsigned k, mp_srcptr xp, mp_size_t n, mp_size_t hn, mp_ptr tp);
__GMP_DECLSPEC int mpn_toom_eval_pm2exp(mp_ptr xp2, mp_ptr xm2, unsigned k, mp_srcptr xp, mp_size_t n, mp_size_t hn, unsigned shift, mp_ptr tp);
__GMP_DECLSPEC int mpn_toom_eval_pm2rexp(mp_ptr rp, mp_ptr rm, unsigned int q, mp_srcptr ap, mp_size_t n, mp_size_t t, unsigned int s, mp_ptr ws);
__GMP_DECLSPEC int mpn_toom_eval_dgr3_pm1(mp_ptr xp1, mp_ptr xm1, mp_srcptr xp, mp_size_t n, mp_size_t x3n, mp_ptr tp);
__GMP_DECLSPEC void mpn_toom_interpolate_16pts(mp_ptr pp, mp_ptr r1, mp_ptr r3, mp_ptr r5, mp_ptr r7, mp_size_t n, mp_size_t spt, int half, mp_ptr wsi);
__GMP_DECLSPEC void mpn_divexact_1(mp_ptr dst, mp_srcptr src, mp_size_t size, mp_limb_t divisor);
__GMP_DECLSPEC void mpn_sqr_basecase(mp_ptr rp, mp_srcptr up, mp_size_t n);
__GMP_DECLSPEC mp_limb_t mpn_divexact_by3c(mp_ptr qp, mp_srcptr xp, mp_size_t n, mp_limb_t ci);
__GMP_DECLSPEC mp_limb_t mpn_divexact_byfobm1(mp_ptr qp, mp_srcptr xp, mp_size_t n, mp_limb_t f, mp_limb_t Bm1of);

__GMP_DECLSPEC mp_limb_t mpn_divrem_euclidean_qr_1(mp_ptr qp, mp_size_t qxn, mp_srcptr xp, mp_size_t n, mp_limb_t d);
__GMP_DECLSPEC mp_limb_t mpn_divrem_euclidean_qr_2(mp_ptr qp, mp_ptr xp, mp_size_t xn, mp_srcptr dp);
__GMP_DECLSPEC mp_limb_t mpn_divrem_euclidean_r_1(mp_srcptr xp, mp_size_t n, mp_limb_t d);
__GMP_DECLSPEC mp_limb_t mpn_rsh_divrem_hensel_qr_1(mp_ptr qp, mp_srcptr xp, mp_size_t n, mp_limb_t d, int s, mp_limb_t cin);
__GMP_DECLSPEC mp_limb_t mpn_divrem_1(mp_ptr qp, mp_size_t qxn, mp_srcptr up, mp_size_t un, mp_limb_t d);
__GMP_DECLSPEC mp_limb_t mpn_divrem_2(mp_ptr qp, mp_size_t qxn, mp_ptr np, mp_size_t nn, mp_srcptr dp);
__GMP_DECLSPEC mp_limb_t mpn_sb_div_qr(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv);
__GMP_DECLSPEC mp_limb_t mpn_dc_div_qr_n(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n, mp_limb_t dinv, mp_ptr tp);
__GMP_DECLSPEC mp_limb_t mpn_dc_div_qr(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv);
__GMP_DECLSPEC mp_limb_t mpn_inv_div_qr(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_srcptr dinv);
__GMP_DECLSPEC mp_limb_t mpn_inv_div_qr_n(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t dn, mp_srcptr inv);
__GMP_DECLSPEC void mpn_invert_trunc(mp_ptr x_new, mp_size_t m, mp_srcptr xp, mp_size_t n, mp_srcptr ap);
__GMP_DECLSPEC void tc4_add_unsigned(mp_ptr rp, mp_size_t * rn, mp_srcptr r1, mp_size_t r1n, mp_srcptr r2, mp_size_t r2n);
__GMP_DECLSPEC void tc4_sub(mp_ptr rp, mp_size_t * rn, mp_srcptr r1, mp_size_t r1n, mp_srcptr r2, mp_size_t r2n);
__GMP_DECLSPEC void tc4_lshift(mp_ptr rp, mp_size_t * rn, mp_srcptr xp, mp_size_t xn, mp_size_t bits);
__GMP_DECLSPEC void mpn_tdiv_qr(mp_ptr qp, mp_ptr rp, mp_size_t qxn, mp_srcptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn);
#define DC_DIVAPPR_Q_N_ITCH(n) ((n)*4 + 64)
#define DC_BDIV_Q_N_ITCH(n) ((n)/2 + 2)
#define DC_BDIV_QR_N_ITCH(n) (n)

#endif
