#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "multiuser_load.h"
void main(int argc, char **argv)
{
	double x;
	x = atof(argv[1]);
	printf("%4.2lf watts = %g dbm/hz \n",x,watts_to_dbmhz(x));
}
double watts_to_dbmhz(double e)
{
        double psd;
        psd = 10*log10((e*1e3)/CHANNEL_BANDWIDTH);
        return psd;
}
