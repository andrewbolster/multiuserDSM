#include <math.h>
#include <stdio.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "multiuser_load.h"

double do_transfer_fn(int,double,double);
extern struct line* list_head;

//typedef struct cmplx {double real,imag;} fcomplex;

#ifdef STANDALONE

main()
{
	int type;
	double length;
	double freq;

	printf("Enter line type > ");
	scanf("%d",&type);
	printf("Enter line length in km > ");
	scanf("%lf",&length);
	printf("Enter frequency in Hz > ");
	scanf("%lf",&freq);

	printf("transfer_fn(%d,%lf,%lf) = %lf\n",type,length,freq,do_transfer_fn(type,length,freq));
}

#endif

#ifndef STANDALONE

double transfer_fn(int line_id,double freq) {

	struct line* current=list_head;
	int type = 2;
	double ret;
	
	while (current->line_id != line_id){
		current = current->next;
	}

	//printf("line_id = %d\tlt = %d nt = %d length = %d\n",current->line_id,current->lt,current->nt,current->length);

	ret = do_transfer_fn(type,current->length/1000,freq);	// change length to kms
	
	return ret;

}

#endif

double insertion_loss(int type, double length, double freq)
{
	double ret;
	//double factor = pow(10,1.2);		// insertion loss is 6dB less than transfer function, apparently!
	double factor = pow(10,1.2);		// insertion loss is 6dB less than transfer function, apparently! quick reduce fext!
	// to reproduce ICT2008 results, change factor to 20dB!! remove h3 from case1 xtalk also	

//	ret = do_transfer_fn(type,length,freq)/factor;
	ret = do_transfer_fn(type,length,freq);		// removved 6dB term to artifically increase xtalk for now!
	return ret;
}


double do_transfer_fn(int type, double length, double f)
{

        double r_oc,a_c,R;
        double l_0,l_inf,b,f_m,L;
        double C_inf,c_0,c_e,C;
        double g_0,g_e,G;

        gsl_complex Z,Y,Z_0,tmp;
        gsl_complex Z_s,Z_l;

        gsl_complex gammad;

        gsl_complex top,denom1,denom2,H;

        double w = 2 * M_PI * f;

        type=3;


        GSL_SET_COMPLEX(&Z_s,100,0);
        GSL_SET_COMPLEX(&Z_l,100,0);


        switch(type) {			// awg 26
                case 1:
			r_oc=286.17578;         // ohms/km
			a_c=0.14769620;

			l_0=675.36888e-6;       // H/km
			l_inf=488.95186e-6;     // H/km

			b=0.92930728;
			f_m=806.33863;          // kHz  

			C_inf=49e-9;            // F/km
			c_0=0;
			c_e=0;

			g_0=43e-9;              // S/km
			g_e=0.70;
                break;
                case 2:			// BT_DWUG
			r_oc=179;               // ohms/km
			a_c=35.89e-3;

			l_0=0.695e-3;   // H/km
			l_inf=585e-6;   // H/km

			b=1.2;
			f_m=1000;               // kHz  

			C_inf=55e-9;            // F/km
			c_0=1e-9;
			c_e=0.1;

			g_0=0.5e-9;             // S/km
			g_e=1.033;
                break;
                case 3:			// awg 24
			r_oc=174.55888;               // ohms/km
			a_c=0.053073;

			l_0=617.29e-6;   // H/km
			l_inf=478.97e-6;   // H/km

			b=1.1529;
			f_m=553.760;               // kHz  

			C_inf=50e-9;            // F/km
			c_0=0.0;
			c_e=0.0;

			g_0=234.87476e-15;             // S/km
			g_e=1.38;
                break;
	}


        R = pow((pow(r_oc,4.0)+a_c*f*f),0.25);
        L = (l_0+l_inf*pow(f*1e-3/f_m,b))/(1+pow(f*1e-3/f_m,b));
        C = C_inf + c_0*pow(f,-c_e);
        G = g_0 * pow(f,g_e);

	printf("f=%lf\tR = %lf\tL = %e\tC = %e G = %e\n",f,R,L,C,G);
	//getchar();

        GSL_SET_COMPLEX(&Z,R,w*L);      // set impedance per unit length
        GSL_SET_COMPLEX(&Y,G,w*C);      // set admittance per unit length

        tmp = gsl_complex_mul(Z,Y);	// gamma = sqrt(Z.Y)

        gammad = gsl_complex_mul_real(gsl_complex_sqrt(tmp),length);    // propogation constant*lenght=gamma*d

        tmp = gsl_complex_div(Z,Y);	// Z_0 = sqrt(Z/Y)

        Z_0 = gsl_complex_sqrt(tmp);    // characteristic impedance

        top = gsl_complex_mul(Z_0,gsl_complex_sech(gammad));    // top = Z_0.sech(gamma*d);

        denom1 = gsl_complex_mul(Z_s,gsl_complex_add(gsl_complex_div(Z_0,Z_l),gsl_complex_tanh(gammad))); // denom1=Z_s(Z_0/Z_l + tanh(gammad))

        denom2 = gsl_complex_mul(Z_0,gsl_complex_add_real(gsl_complex_mul(gsl_complex_div(Z_0,Z_l),gsl_complex_tanh(gammad)),1));

        H = gsl_complex_div(top,gsl_complex_add(denom1,denom2));

        //return gsl_complex_abs(H);
        return gsl_complex_abs(gsl_complex_mul(H,H));


/*
  double r_0c, a_c, l_0, l_inf, f_m, b, R, L, C;
  double Z_L, Z_S, omega, alpha, tmp;
  fcomplex Z_0, gammad, hypsin, hypcos, expgammad, expngammad;
  fcomplex Z_0sinh, Z_0cosh, denom, H;

  type=2;

  switch (type) {
  case 1:
    r_0c = 0.409e3;
    a_c = 0.3822;
    l_0 = 0.6075e-3;
    l_inf = 0.5e-3;
    f_m = 0.6090e3;
    b = 5.269;
    C = 40e-9;
    break;
  case 2:
    r_0c = 0.28e3;
    a_c = 0.0969;
    l_0 = 0.5873e-3;
    l_inf = 0.426e-3;
    f_m = 0.7459e3;
    b = 1.385;
    C = 50e-9; 
    break;
  case 3:
    r_0c = 0.1792e3;
    a_c = 0.0561;
    l_0 = 0.6746e-3;
    l_inf = 0.5327e-3;
    f_m = 0.6647e3;
    b = 1.195;
    C = 50e-9;
    break;
  case 4:
    r_0c = 0.113e3;
    a_c = 0.0257;
    l_0 = 0.6994e-3;
    l_inf = 0.4772e-3;
    f_m = 0.2658e3;
    b = 1.0956;
    C = 45e-9;
    break;
  case 5:
    r_0c = 0.0551e3;
    a_c = 0.0094;
    l_0 = 0.7509e-3;
    l_inf = 0.5205e-3;
    f_m = 0.1238e3;
    b = 0.9604;
    C = 40e-9;
    break;
  }

  R = pow(pow(r_0c, 4.0) + a_c*freq*freq, 0.25);
  L = (l_0 + l_inf * pow(freq*1e-3/f_m, b)) / (1 + pow(freq*1e-3/f_m, b));
  C = C;
  Z_L = 100;
  Z_S = 100;
  omega = 2 * M_PI * freq;

  alpha = atan2(-R, omega*L);
  tmp = sqrt(sqrt(omega*L*omega*L + R*R)/(omega*C));
  Z_0.real = tmp*cos(alpha/2.0);
  Z_0.imag = tmp*sin(alpha/2.0);

  alpha = atan2(R, -omega*L);
  tmp = sqrt(sqrt(omega*L*omega*L + R*R)*omega*C);
  gammad.real = length*tmp*cos(alpha/2.0);
  gammad.imag = length*tmp*sin(alpha/2.0);

  expgammad.real = exp(gammad.real)*cos(gammad.imag);
  expgammad.imag = exp(gammad.real)*sin(gammad.imag);
  expngammad.real = exp(-gammad.real)*cos(gammad.imag);
  expngammad.imag = exp(-gammad.real)*sin(-gammad.imag);

  hypsin.real = 0.5 * (expgammad.real - expngammad.real);
  hypsin.imag = 0.5 * (expgammad.imag - expngammad.imag);
  hypcos.real = 0.5 * (expgammad.real + expngammad.real);
  hypcos.imag = 0.5 * (expgammad.imag + expngammad.imag);

  Z_0sinh.real = Z_0.real*hypsin.real - Z_0.imag*hypsin.imag;
  Z_0sinh.imag = Z_0.real*hypsin.imag + Z_0.imag*hypsin.real;
  Z_0cosh.real = Z_0.real*hypcos.real - Z_0.imag*hypcos.imag;
  Z_0cosh.imag = Z_0.real*hypcos.imag + Z_0.imag*hypcos.real;

  denom.real = Z_0.real*(Z_L * hypcos.real + Z_0sinh.real) -
    Z_0.imag*(Z_L * hypcos.imag + Z_0sinh.imag);
  denom.imag = Z_0.real*(Z_L * hypcos.imag + Z_0sinh.imag) +
    Z_0.imag*(Z_L * hypcos.real + Z_0sinh.real);
  denom.real = (1/Z_L)*(Z_S*(Z_0cosh.real + Z_L*hypsin.real) + denom.real);
  denom.imag = (1/Z_L)*(Z_S*(Z_0cosh.imag + Z_L*hypsin.imag) + denom.imag);

  tmp = denom.real*denom.real + denom.imag*denom.imag;
  H.real = (Z_0.real * denom.real + Z_0.imag * denom.imag)/tmp;
  H.imag = (Z_0.imag * denom.real - Z_0.real * denom.imag)/tmp;

  
  return(sqrt(H.real*H.real + H.imag*H.imag)); 
*/
}

