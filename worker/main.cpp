#include "mpir/mpir.h"
#include "mpir_CUDA/mpir_CUDA.h"
int main()
{
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
	
	mpf_mul(f5, f4, f1);// f2=f4 * f1;
	d5 = mpf_get_d(f5);

	mpf_div(f6,f4,f1);// f5=f4/f1;13.58016
	d6 = mpf_get_d(f6);
/*    gmp_sprintf(sBuf,"\n%Ff,%ld\n",f5,ee);*/
    printf("%ld %f %s\n",ee,d4,sBuf);
    mpf_clear(f1);

    //i=IsCUDA_Supported(1);
	//sayhello();
    return i;
}
