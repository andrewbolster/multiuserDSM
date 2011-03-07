#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "multiuser_load.h"


int main(int argc, char **argv)
{

	double x;
	x = atof(argv[1]);
	printf("%4.2lf dbm/hz = %4.2e watts = %4.2lf mW\n",x,dbmhz_to_watts(x),1e3*dbmhz_to_watts(x));

}


double dbmhz_to_watts(double psd)
{

        double watts;

        if (psd == 0)                   // using 0dBm/Hz to represent off state as this psd level is unrealistic anyway
                return 0;

        watts = pow(10,psd/10) * 1e-3 * CHANNEL_BANDWIDTH;

        return watts;

}

