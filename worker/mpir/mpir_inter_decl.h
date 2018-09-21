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
#endif
