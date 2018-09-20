#ifndef __MPIR_INTER_DECL_H__
#define __MPIR_INTER_DECL_H__
#include"mpir.h"
int __gmp_extract_double(mp_ptr rp, double d);
double mpn_get_d(mp_srcptr ptr, mp_size_t size, mp_size_t sign, long exp);
#endif