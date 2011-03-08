#include <stdio.h>
#include <math.h>
#include "multiuser_load.h"
main()
{
	double SNR;
	double p_e;
	while(1) {
		printf("Enter SNR in dB >");
		scanf("%lf",&SNR);
		printf("Enter p_e in dB >");
		scanf("%lf",&p_e);
//		printf("Gamma = %lf\n",calc_gamma(p_e,SNR));	
	}
}
