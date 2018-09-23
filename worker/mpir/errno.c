/* gmp_errno, __gmp_exception -- exception handling and reporting.

   THE FUNCTIONS IN THIS FILE, APART FROM gmp_errno, ARE FOR INTERNAL USE
   ONLY.  THEY'RE ALMOST CERTAIN TO BE SUBJECT TO INCOMPATIBLE CHANGES OR
   DISAPPEAR COMPLETELY IN FUTURE GNU MP RELEASES.

Copyright 2000, 2001, 2003 Free Software Foundation, Inc.

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
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdlib.h>
#include "mpir.h"

#ifdef MPIR_CUDA_ACC
__GMP_DECLSPEC int gmp_errnod = 0;
__GMP_DECLSPEC const int __gmpd_0 = 10;
__GMP_DECLSPEC int __gmp_junkd = 0;
__GMP_DECLSPEC
void
__gmp_exception(int error_bit)
{
	int a = 9;
	++a;
	__gmp_junkd = 10 / (__gmpd_0 - a);
}
#else
int gmp_errno = 0;
const int __gmp_0=10;
int __gmp_junk = 0;
/* The deliberate divide by zero triggers an exception on most systems.  On
   those where it doesn't, for example power and powerpc, use abort instead.

   Enhancement: Perhaps raise(SIGFPE) (or the same with kill()) would be
   better than abort.  Perhaps it'd be possible to get the BSD style
   FPE_INTDIV_TRAP parameter in there too.  */
__GMP_DECLSPEC
void
__gmp_exception (int error_bit)
{
    int a=9;
    ++a;
  __gmp_junk = 10 / (__gmp_0-a);
#ifndef MPIR_CUDA_ACC
  abort();
#endif // !MPIR_CUDA_ACC
}
#endif

/* These functions minimize the amount of code required in functions raising
   exceptions.  Since they're "noreturn" and don't take any parameters, a
   test and call might even come out as a simple conditional jump.  */
__GMP_DECLSPEC
void
__gmp_sqrt_of_negative (void)
{
  __gmp_exception (GMP_ERROR_SQRT_OF_NEGATIVE);
}
__GMP_DECLSPEC
void
__gmp_divide_by_zero (void)
{
  __gmp_exception (GMP_ERROR_DIVISION_BY_ZERO);
}
