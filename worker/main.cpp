#include "mpir/mpir.h"
#include "mpir_CUDA/mpir_CUDA.h"
int main()
{
    int i=0;

    char sBuf[3600];sBuf[0]=0;
    mpf_t f1,f2,f3,f4,f5;
    mpf_init(f1);
    mpf_init(f2);
    mpf_init(f3);
    mpf_init(f4);
    mpf_init(f5);
    mpf_set_d(f1,123.456);
    
    /*mpf_mul_ui(f2, f1, 10);//f2=f1*10,1234.56
    mpf_add(f3,f2,f1);//f3=f2+f1,1358.016
    mpf_div_ui(f4,f3,10);//f4=f3/10;135.8016
    mpir_si ee;
    mpf_get_d_2exp(&ee, f4);
    mpf_set_d(f1,10.0);
    mpf_div(f5,f4,f1);// f5=f4/f1;13.58016
    gmp_sprintf(sBuf,"\n%Ff,%ld\n",f5,ee);*/
    printf("%s",sBuf);
    mpf_clear(f1);

    //i=IsCUDA_Supported(1);
	//sayhello();
    return i;
}