#include "mpir/mpir.h"
#include "mpir_CUDA/mpir_CUDA.h"
#include <math.h>
#include <conio.h>
double MPIRWorking_HOST(double fStartValue);
int main()
{
#pragma message(VAR_NAME_VALUE(GMP_NUMB_BITS))
    int i=0;

	char sBuf[3600]; sBuf[0] = 0;
    double d1,d2,d3,d4,d5,d6,d7;
    mpf_t f1,f2,f3,f4,f5,f6,f7;
    mpf_init(f1);
    mpf_init(f2);
    mpf_init(f3);
    mpf_init(f4);
	mpf_init(f5);
	mpf_init(f6);
	mpf_init(f7);
    mpf_set_d(f1,123.456);

	d1 = mpf_get_d(f1); 
    mpf_mul_ui(f2, f1, 10);//f2=f1*10,1234.56
    d2 = mpf_get_d(f2);
    mpf_add(f3,f2,f1);//f3=f2+f1,1358.016
    d3 = mpf_get_d(f3);
    
    mpf_sub(f4,f3,f1);//f4=f3-f1
    d4 = mpf_get_d(f4);
    mpf_div_ui(f4,f3,10);//f4=f3/10;135.8016
    d4 = mpf_get_d(f4);
    mp_exp_t ee;
    mpf_get_d_2exp(&ee, f4);
    mpf_set_d(f1,10.0);

	mpf_mul(f5, f4, f1);// f2=f4 * f1;1358.016
	d5 = mpf_get_d(f5);

	mpf_div(f6, f4, f1);// f5=f4/f1;13.58016
	d6 = mpf_get_d(f6);
	mpf_div_ui(f7, f6, 10);// f5=f4/f1;13.58016
	d7 = mpf_get_d(f7);
	if (mpf_cmp_d(f7, 10) > 0)
	{
		printf("big");
	}
/*    gmp_sprintf(sBuf,"\n%Ff,%ld\n",f5,ee);*/
    printf("%ld %f %s\n",ee,d4,sBuf);
    mpf_clear(f1);
#define EPSILON 1e-12
    bool b;
    d1=123.4561;
    d2=MPIRWorking_HOST(d1);
    b=(d1==d2);
    d3=fabs(d1)-fabs(d2);
    b=(d3<EPSILON);
    d1=124;
    d2=MPIRWorking_HOST(d1);
    b=(d1==d2);
    d3=fabs(d1)-fabs(d2);
    b=(d3<EPSILON);
	//i=IsCUDA_Supported(1);
	printf("Host-->MPIR sizes:sizeof(mp_limb_t)=%d,sizeof(mp_limb_signed_t)=%d,sizeof(mp_size_t)=%d,sizeof(mp_exp_t)=%d\n", sizeof(mp_limb_t), sizeof(mp_limb_signed_t), sizeof(mp_size_t), sizeof(mp_exp_t));
	sayhello(); getch();
    return i;
}
double MPIRWorking_HOST(double fStartValue)
{
    mpf_t f1,f2,f3,f4,f5;
    mpf_init(f1);
    mpf_init(f2);
    mpf_init(f3);
    mpf_init(f4);
    mpf_init(f5);
    mpf_set_d(f1,fStartValue);
    mpf_mul(f2, f1, f1);// f2=f1 * f1;15241.383936
    mpf_div_ui(f3, f2, 10);//f3=f2/10;1524.1383936
    mpf_set_d(f2,1000L);
    mpf_mul(f3,f3,f2);//f3=f3*f2;1524138.3936
    mpf_div(f4,f3,f1);//f4=f3/f1;12345.6
    mpf_div_ui(f5,f4,100L);//f4=f3/f1;123.456
    return mpf_get_d(f5);
}
