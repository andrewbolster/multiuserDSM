#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "multiuser_load.h"
static double l_erfc(double);
double gammaln(double);
void gser(double,double,double *,double *);
void gcf(double,double,double *,double *);
double gammp(double,double);
double gammq(double,double);
/*
void main(int argc, char** argv)
{
        double x =0;
        printf("Enter double value >");
        scanf("%lf",&x);
        printf("erfc(%lf) = %lf\n",x,erfc(x));
}*/
double symerr(double SNR, int b)
{
        double param,param1,symerr;
        double M = pow(2,b);
        param = sqrt(3 * SNR/(2 * (M-1)));
        //printf("sqrt(3 * SNR/(2 * (M-1))) = %lf\n",param);
        if (param < 0)
                param = -param;
        //printf("1/sqrt(M) = %lf\n",1/sqrt(M));
        //printf("erfc(%lf) = %lf\n",param,erfc(param));
        param1 = 1 - ((1 - (1/sqrt(M))) * l_erfc(param));
        //printf("1- ((1 - (1/sqrt(M))) * erfc(param)) = %lf\n",param1);
        symerr = 1 - pow(param1,2);
        return symerr;
}
double symerr_inf(double SNR, double b)
{
        double param,param1,symerr;
        double M = pow(2,b);
        param = sqrt(3 * SNR/(2 * (M-1)));
        //printf("sqrt(3 * SNR/(2 * (M-1))) = %lf\n",param);
        if (param < 0)
                param = -param;
        //printf("1/sqrt(M) = %lf\n",1/sqrt(M));
        //printf("erfc(%lf) = %lf\n",param,erfc(param));
        param1 = 1 - ((1 - (1/sqrt(M))) * l_erfc(param));
        //printf("1- ((1 - (1/sqrt(M))) * erfc(param)) = %lf\n",param1);
        symerr = 1 - pow(param1,2);
        return symerr;
}
/* Interface to matlab added May, 2000 */
/* Can be compiled as a mex file */
/* Numerical recipies in Pascal conversion of gammln */
/* David I. Laurenson, July 31st, 1997 */
/* Department of Electrical Engineering, The University of Edinburgh */
/* gammln returns the value ln(Gamma(xx)) for xx>0.  Full accuracy is */
/* obtained for xx>1.  For 0<xx<1 the reflection formula can be used. */
/*                      Pi                   Pi z          */
/* Gamma(1-z) = ------------------ = --------------------  */
/*              Gamma(z) sin(Pi z)   Gamma(1+z) sin(Pi z)  */
#define stp 2.50662827465
double gammln(double xx)
{
  double x, tmp, ser;
  x = xx - 1.0;
  tmp = x + 5.5;
  tmp = (x + 0.5) * log(tmp) - tmp;
  ser = 1.0 + 76.18009173/(x+1.0) - 86.505323033/(x+2.0) +
    24.01409822/(x+3.0) - 1.231739516/(x+4.0) + 0.120858003e-2/(x+5.0)
      - 0.536382e-5/(x+6.0);
  return(tmp + log(stp*ser));
}
#undef stp
/* Numerical recipies in Pascal conversion of gser */
/* David I. Laurenson, June 24th, 1993 */
/* Department of Electrical Engineering, The University of Edinburgh */
/* Calculates the incomplete gamma function, P(a,x) for real a and x */
void gser(double a, double x, double *gamser, double *gln)
{
  int itmax = 100;
  double eps = 3.0e-7;
  int n;
  double sum, del, ap;
  *gln = gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) {
      fprintf(stderr, "GSER - x is less than zero\n");
      exit(-1);
    }
    *gamser = 0.0;
  }
  else {
    ap = a;
    sum = 1.0 / a;
    del = sum;
    for(n=1; n<=itmax; n++) {
      ap += 1.0;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*eps)
	break;
    }
    if (n>itmax)
      fprintf(stderr, "GSER - Too many iterations\n");
    *gamser = sum * exp(-x + a * log(x) - *gln);
  }
}
/* Numerical recipies in Pascal conversion of gcf */
/* David I. Laurenson, June 24th, 1993 */
/* Department of Electrical Engineering, The University of Edinburgh */
/* Calculates the incomplete gamma function, Q(a,x) for real a and x */
void gcf(double a, double x, double *gammcf, double *gln)
{
  int itmax = 100;
  double eps = 3.0e-7;
  int n;
  double gold, g, fac, b1, b0, anf, ana, an, a1, a0;
  *gln = gammln(a);
  gold = 0.0;
  g = 0.0; // compiler warning
  a0 = 1.0;
  a1 = x;
  b0 = 0.0;
  b1 = 1.0;
  fac = 1.0;
  for(n=1; n<=itmax; n++) {
    an = (double) n;
    ana = an - a;
    a0 = (a1 + a0 * ana) * fac;
    b0 = (b1 + b0 * ana) * fac;
    anf = an * fac;
    a1 = x * a0 + anf * a1;
    b1 = x * b0 + anf * b1;
    if (a1 != 0.0) {
      fac = 1.0 / a1;
      g = b1 * fac;
      if (fabs((g - gold) / g) < eps)
	break;
      gold = g;
    }
  }
  if (n>itmax) {
    fprintf(stderr, "GCF - too many iterations\n");
    getchar();
  }
  *gammcf = exp(-x + a * log(x) - *gln) * g;
}
/* Numerical recipies in Pascal conversion of gammp */
/* David I. Laurenson, June 24th, 1993 */
/* Department of Electrical Engineering, The University of Edinburgh */
/* Calculates the incomplete gamma function, P(a,x) for real a and x */
double gammp(double a,double x)
{
  double gamser, gammcf, gln;
  if ((x < 0.0) || (a <= 0.0)) {
    fprintf(stderr, "GAMMP - Invalid inputs: (a,x) = (%lg,%lg)\n", a,
	    x);
    exit(-1);
  }
  if (x < a+1.0) {
    gser(a, x, &gamser, &gln);
    return(gamser);
  }
  else {
    gcf(a, x, &gammcf, &gln);
    return(1.0 - gammcf);
  }
}
/* Numerical recipies in Pascal conversion of gammq */
/* David I. Laurenson, June 24th, 1993 */
/* Department of Electrical Engineering, The University of Edinburgh */
/* Calculates the incomplete gamma function, Q(a,x) for real a and x */
double gammq(double a,double x)
{
  double gamser, gammcf, gln;
  if ((x < 0.0) || (a <= 0.0)) {
    fprintf(stderr, "GAMMQ - Invalid inputs: (a,x) = (%lg,%lg)\n", a,
	    x);
    exit(-1);
  }
  if (x < a+1.0) {
    gser(a, x, &gamser, &gln);
    return(1.0 - gamser);
  }
  else {
    gcf(a, x, &gammcf, &gln);
    return(gammcf);
  }
}
static double l_erfc(double x)
{
  if (x<0.0) {
    return(1+gammp(0.5, x*x));
  }
  return(gammq(0.5, x*x));
}
double calc_gamma(double pe_target,double SNR)
{
        //double SNR = 2e9;       // doesnt matter what exact value is, sorry did i say it doesnt? it does!
        double gamma = 9.95 - (pe_target/1e-7)/2500;
        double b;
        //double M;
        double max_err = pe_target/1000;
        double max_iters = 3000;
        int i = 1;
	SNR = pow(10,SNR/10);
        do {
                printf("gamma = %lf\n",gamma);
                b = log(1 + SNR/pow(10,gamma/10))/log(2.0);
                //printf("b = log_2(1 + SNR/GAMMA) = %lf\n",b);
                //M = pow(2,b);
                //printf("M = 2^b = %lf\n",M);
                //printf("symerrprob(%lf,%lf) = %e\n",SNR,M,symerr(SNR,M));
                gamma = gamma + 0.5 * (symerr_inf(SNR,b)/pe_target - 1);
                if (i++ > max_iters) {
                        printf("max_iters exceeded in calc_gamma\n");
                        return -1;
                }
        } while(fabs(symerr_inf(SNR,b) - pe_target) > max_err);
        printf("Found gamma = %lf in %d iterations\n",gamma,i);
        return gamma;
}
