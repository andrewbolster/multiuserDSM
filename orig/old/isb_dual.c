#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "multiuser_load.h"
#include <string.h>
#include <float.h>
#include <stdbool.h>

#define P_MAT p_mat[2][2][2][MAXBITSPERTONE+1][MAXBITSPERTONE+1][DMTCHANNELS]
#define B_MAT b_mat[2][DMTCHANNELS]
#define S_MAT s_mat[2][DMTCHANNELS]

#define RATE_GOOD(x,y,z) (abs(rate(x,y,p,z) - p->rate_target[x][y]) < e)
#define RATE_OK(x,y,z) (rate(x,y,p,z) > p->rate_target[x][y])
#define RATE_HIGH(x,y,z) ((rate(x,y,p,z) - p->rate_target[x][y]) > e)
#define RATE_LOW(x,y,z) ((rate(x,y,p,z) - p->rate_target[x][y]) < -e)
#define ALL_RATES_GOOD (RATE_GOOD(0,0,true) && RATE_GOOD(0,1,true) && RATE_GOOD(1,0,true) && RATE_GOOD(1,1,true))

#define GAMMA0 9.95
#define GAMMA1 6.02

struct isb_2_dual_params {

	double P_MAT;
	int B_MAT;
	int S_MAT;
	double w1;
	double w1_max;
	double w1_min;
	double w2;
	double w2_max;
	double w2_min;
	double w3;
	double w3_max;
	double w3_min;
	double w4;
	double w4_max;
	double w4_min;
	double l1;
	double l1_min;
	double l1_max;
	double l2;
	double l2_min;
	double l2_max;
	/*
	int rate0s0_target;
	int rate0s1_target;
	int rate1s0_target;
	int rate1s1_target;
	*/
	int rate_target[2][2];
	double p0_budget;
	double p1_budget;
	double p0_budget_min;
	double p1_budget_min;
	double p0_budget_max;
	double p1_budget_max;
};

static void optimise_w1(struct isb_2_dual_params *p);
static void optimise_w2(struct isb_2_dual_params *p);
static void optimise_w(struct isb_2_dual_params *p);
static void optimise_l1(struct isb_2_dual_params *p);
static void optimise_l2(struct isb_2_dual_params *p);
static void isb_init_params(struct isb_2_dual_params *p);
static void optimise_p(struct isb_2_dual_params *p);
static void optimise_p2(struct isb_2_dual_params *p);
static double l_k(int b1,int b2,int s0,int s1,int channel,struct isb_2_dual_params *p);
static void isb_init_p_matrix(double P_MAT);
static void psd(double P_MAT,int b1,int b2,int s0,int s1,int channel);
static void _psd(int line_id,int s0,int s1,int b1,int b2,double P_MAT,int channel,double gain);
static double _cnr(int line_id,double P_MAT,int channel,double gain,int b1,int b2,int s0,int s1);
static double osb_fext(int line_id,double P_MAT,int channel,int b1,int b2,int s0,int s1);
static int rate(int line_id,int service,struct isb_2_dual_params *p,bool print);
static double tot_pow(int line_id,struct isb_2_dual_params *p,bool print);
static void isb_dual_init_lines(struct isb_2_dual_params *p);

void isb_init_params(struct isb_2_dual_params *p)
{

        memset(p->b_mat,0,sizeof(int)*lines*DMTCHANNELS);
        memset(p->p_mat,0,sizeof(double)*lines*(MAXBITSPERTONE+1)*(MAXBITSPERTONE+1)*DMTCHANNELS);
        memset(p->s_mat,0,sizeof(int)*lines*DMTCHANNELS);

        p->w1_min=0.0;
        p->w1_max=1.0;

        p->w2_min=0.0;
        p->w2_max=1.0;
        
	p->w3_min=0.0;
        p->w3_max=1.0;

        p->w4_min=0.0;
        p->w4_max=1.0;

        p->l1_min=0.0;
        p->l1_max=1.0;

        p->l2_min=0.0;
        p->l2_max=1.0;
/*
        p->rate0s0_target=1092;
        p->rate0s1_target=504;

        p->rate1s0_target=844;
        p->rate1s1_target=272;
*/

	//p->rate_target[0][0] = 1075;
	p->rate_target[0][1] = 600;
	//p->rate_target[1][0] = 530;
	p->rate_target[1][1] = 600;
       
	p->p0_budget=0.1;
        p->p1_budget=0.1;
	p->p0_budget_max=0.1;
	p->p0_budget_min=0.0;
	p->p1_budget_max=0.1;
	p->p1_budget_min=0.0;


        isb_init_p_matrix(p->p_mat);

}


struct isb_2_dual_params *isb_2_dual()
{

	struct isb_2_dual_params *p;

	FILE *fp;

	fp=fopen("isb_dual_test.txt","w");
	if (fp==NULL)
		exit(1);

	p = (struct isb_2_dual_params *)malloc(sizeof(struct isb_2_dual_params));

	if (p==NULL) {
		printf("malloc sucks\n");
		exit(1);
	}

	
	isb_init_params(p);

	//printf("line0 service0 target rate = %d\n",p->rate0s0_target);
	//printf("line0 service1 target rate = %d\n",p->rate0s1_target);

	optimise_w(p);
/*
	p->w1=0.50781250;
	p->w2=0.35156250;
	p->w3=0.64062500;
	p->w4=0.45;
*/
	//optimise_l1(p);

/*
	fprintf(fp,"w1\tw2\tl1\tl2\t\tb0s0\tb0s1\tp0\tb1s0\tb1s1\tp1\n");
	for (p->w1=0.0;p->w1<=1.0;p->w1+=0.05) {
		for (p->w2=0.0;p->w2<=1.0;p->w2+=0.05) {
			optimise_l1(p);	
			fprintf(fp,"%4.2lf\t%4.2lf\t%6.4lf\t%6.4lf\t%d\t%d\t%4.2lf\t%d\t%d\t%4.2lf\n",p->w1,p->w2,p->l1,p->l2,
									    rate(0,0,p),rate(0,1,p),tot_pow(0,p,false),
								            rate(1,0,p),rate(1,1,p),tot_pow(1,p,false));
			printf("%4.2lf\t%4.2lf\t%6.4lf\t%6.4lf\t%d\t%d\t%4.2lf\t%d\t%d\t%4.2lf\n",p->w1,p->w2,p->l1,p->l2,
									    rate(0,0,p),rate(0,1,p),tot_pow(0,p,false),
								            rate(1,0,p),rate(1,1,p),tot_pow(1,p,false));
		}
	}
*/



	
	rate(0,0,p,true);
	rate(0,1,p,true);
	rate(1,0,p,true);
	rate(1,1,p,true);
	tot_pow(0,p,true);
	tot_pow(1,p,true);

	isb_dual_init_lines(p);

	return p;
}

/*
void optimise_w1(struct osb_2_dual_params *p)		// w1 = line0service0 weight
{

	int e = 5;
	while(abs(rate(0,0,p) - p->rate0s0_target) > e) {
		p->w1=(p->w1_max+p->w1_min)/2;
		printf("w1=%lf\n",p->w1);

		optimise_w2(p);

		if (rate(0,0,p) > p->rate0s0_target) 
			p->w1_max=p->w1;
		else
			p->w1_min=p->w1; 

	}

}

void optimise_w2(struct osb_2_dual_params *p)	// line0service 1 weight
{

	int e = 5;

	while(abs(rate(0,1,p) - p->rate0s1_target) > e) {
		p->w2=(p->w2_max+p->w2_min)/2;
		printf("w2=%lf\n",p->w2);

		optimise_w3(p);

		if (rate(0,1,p) > p->rate0s1_target) 
			p->w2_max=p->w2;
		else
			p->w2_min=p->w2; 

	}

}

void optimise_w3(struct osb_2_dual_params *p)		// line1service0 weight
{

	int e = 5;
	while(abs(rate(1,0,p) - p->rate1s0_target) > e) {
		p->w3=(p->w3_max+p->w3_min)/2;
		printf("w3=%lf\n",p->w3);

		optimise_w4(p);

		if (rate(1,0,p) > p->rate1s0_target) 
			p->w3_max=p->w3;
		else
			p->w3_min=p->w3; 

	}

}

void optimise_w4(struct osb_2_dual_params *p)	// line1service 1 weight
{

	int e = 5;

	while(abs(rate(1,1,p) - p->rate1s1_target) > e) {
		p->w4=(p->w4_max+p->w4_min)/2;
		printf("w4=%lf\n",p->w4);

		optimise_l1(p);

		if (rate(1,1,p) > p->rate1s1_target) 
			p->w4_max=p->w4;
		else
			p->w4_min=p->w4; 

	}

}
*/

#ifdef notdef

void optimise_w(struct osb_2_dual_params *p)
{
	int e=20;

	while(1) {
		p->w1=(p->w1_max+p->w1_min)/2;
                p->w2=(p->w2_max+p->w2_min)/2;
		p->w3=(p->w3_max+p->w3_min)/2;
                p->w4=(p->w4_max+p->w4_min)/2;

		optimise_l1(p);

		printf("w1=%.8lf w2=%.8lf w3=%.8lf w4=%.8lf l1=%4.2lf l2=%4.2lf\n",p->w1,p->w2,p->w3,p->w4,p->l1,p->l2);
		/*if (RATE_HIGH(0,0) && RATE_HIGH(0,1) && RATE_HIGH(1,0) && RATE_HIGH(1,1)) {
			printf("All rates are too high\n");
			//adjust_p(p);
		}*/

			
		if (RATE_GOOD(0,1) && RATE_GOOD(1,1)) {
			printf("rates are good!\n");
			return;
		}

		if (rate(0,1,p,true) > p->rate_target[0][1])	{
			p->w1_max=p->w1;
			//p->w3_min=p->w3;
		}
		else {
			p->w1_min=p->w1;
			//p->w3_max=p->w3;
		}
		if (rate(1,1,p,true) > p->rate_target[1][1])	{
			p->w2_max=p->w2;
			//p->w4_min=p->w4;
		}
		else {
			p->w2_min=p->w2;
			//p->w4_max=p->w4;
		}
		/*if (rate(1,0,p) > p->rate_target[1][0])	
			p->w3_max=p->w3;
		else
			p->w3_min=p->w3;
		if (rate(1,1,p) > p->rate_target[1][1])	
			p->w4_max=p->w4;
		else
			p->w4_min=p->w4;*/
/*
		if (RATE_GOOD(0,1) && RATE_HIGH(1,1)) {
			printf("rate 0,1 is good but 1,1 is too high!\n");
			p->w4_min=p->w4;	// increase weight of b2s0
		}
		if (RATE_GOOD(1,1) && RATE_HIGH(0,1)) {
			printf("rate 1,1 is good but 0,1 is too high!\n");
			p->w3_min=p->w3;	// increase weight of b1s0
		}	
		
		if (RATE_GOOD(0,1) && RATE_LOW(1,1)) {
			printf("rate 0,1 is good but 1,1 is too low!\n");
			p->w4_max=p->w4;	// decrease weight of b2s0
			//return;
		}
		if (RATE_GOOD(1,1) && RATE_LOW(0,1)) {
			printf("rate 1,1 is good but 0,1 is too low!\n");
			p->w3_max=p->w3;	// decrease weight of b1s0
		}
*/

		tot_pow(0,p,true);
		tot_pow(1,p,true);
	
	}

}

#endif


void optimise_w(struct isb_2_dual_params *p)
{
	int e=15;
	
	p->w1=(p->w1_max+p->w1_min)/2;
        p->w2=(p->w2_max+p->w2_min)/2;
	p->w3=(p->w3_max+p->w3_min)/2;
        p->w4=(p->w4_max+p->w4_min)/2;

	while(1) {
		
		p->w1_min=0.0;
		p->w1_max=1.0;
		while(!RATE_GOOD(0,1,true)) {
			p->w1=(p->w1_max+p->w1_min)/2;
			printf("Currently adjusting w1\n");
			printf("w1=%.8lf w2=%.8lf w3=%.8lf w4=%.8lf l1=%4.2lf l2=%4.2lf\n",p->w1,p->w2,p->w3,p->w4,p->l1,p->l2);
			optimise_l1(p);
			if (RATE_HIGH(0,1,true))
				p->w1_max=p->w1;
			else
				p->w1_min=p->w1;

		}
		printf("w1 converged\n");
		p->w2_min=0.0;
		p->w2_max=1.0;
		while(!RATE_GOOD(1,1,true)) {
			p->w2=(p->w2_max+p->w2_min)/2;
			printf("Currently adjusting w2\n");
			printf("w1=%.8lf w2=%.8lf w3=%.8lf w4=%.8lf l1=%4.2lf l2=%4.2lf\n",p->w1,p->w2,p->w3,p->w4,p->l1,p->l2);
			optimise_l1(p);

			if(RATE_HIGH(1,1,true))
				p->w2_max=p->w2;
                        else
                                p->w2_min=p->w2;
		}
		printf("w2 converged\n");
		/*
		p->w3_min=0.0;
                p->w3_max=1.0;
                while(!RATE_GOOD(0,0,true)) {
                        p->w3=(p->w3_max+p->w3_min)/2;
                        printf("Currently adjusting w3\n");
                        printf("w1=%.8lf w2=%.8lf w3=%.8lf w4=%.8lf l1=%4.2lf l2=%4.2lf\n",p->w1,p->w2,p->w3,p->w4,p->l1,p->l2);
                        optimise_l1(p);

                        if(RATE_HIGH(0,0,true))
                                p->w3_max=p->w3;
                        else
                                p->w3_min=p->w3;
                }
                printf("w3 converged\n");
	
		p->w4_min=0.0;
		p->w4_max=1.0;
		while(!RATE_GOOD(1,0,true)) {
			p->w4=(p->w4_max+p->w4_min)/2;
			printf("Currently adjusting w4\n");
			printf("w1=%.8lf w2=%.8lf w3=%.8lf w4=%.8lf l1=%4.2lf l2=%4.2lf\n",p->w1,p->w2,p->w3,p->w4,p->l1,p->l2);
			optimise_l1(p);

			if(RATE_HIGH(1,0,true))
				p->w4_max=p->w4;
                        else
                                p->w4_min=p->w4;
		}
		printf("w4 converged\n");
		*/
		
		if (RATE_GOOD(0,1,false) && RATE_GOOD(1,1,false)) {
			printf("All good!\n");
			return;
		}
		/*
		if (ALL_RATES_GOOD) {
                        printf("All good!\n");
                        return;
                }
		*/
		//e=e-e/4;
	
	}
}
/*
void adjust_p(struct osb_2_dual_params *p)
{

	while (1) {
		p->p0_budget=(p->p0_budget_max+p->p0_budget_min)/2;
		p->p1_budget=(p->p1_budget_max+p->p1_budget_min)/2;
		
		optimise_l1(p);

		if (ALL_RATES_GOOD) {
			return;
		}

		

	}


}
*/

void optimise_l1(struct isb_2_dual_params *p)
{

	double pow,last;

	p->l1=1.0;
	p->l1_min=0.0;
	do {
		p->l1*=2;
		optimise_l2(p);
	} while(tot_pow(0,p,false) > p->p0_budget);
	p->l1_max=p->l1;
	//printf("First value of l1 to meet power constraint = %lf\n",p->l1_max);

	while(1) {
		p->l1=(p->l1_max+p->l1_min)/2;
		
		optimise_l2(p);

		pow=tot_pow(0,p,false);

		if (pow > p->p0_budget) {
			p->l1_min=p->l1;
		}
		else {
			p->l1_max=p->l1;
			if (pow == last) {
			//if (fabs(pow-last) < 1e-2)
				/*printf("l1 converged\n");
				tot_pow(0,p,true);
                                tot_pow(1,p,true);
                                rate(0,0,p);
                                rate(0,1,p);
                                rate(1,0,p);
                                rate(1,1,p);*/
				return;
			}
		}
		last=pow;
	}
}


void optimise_l2(struct isb_2_dual_params *p)
{

	double pow,last;

	p->l2=1.0;
	p->l2_min=0.0;
	do {
		p->l2*=2;
		optimise_p2(p);
	} while(tot_pow(1,p,false) > p->p1_budget);
	p->l2_max=p->l2;	
	//printf("First value of l2 to meet power constraint = %lf\n",p->l2_max);

	while(1) {
		p->l2=(p->l2_max+p->l2_min)/2;
		optimise_p2(p);

		pow = tot_pow(1,p,false);	

		if (pow > p->p1_budget) {
			p->l2_min=p->l2;
		}
		else {
			p->l2_max=p->l2;
			if (pow==last) {
			//if (fabs(pow-last) < 1e-2)
				/*printf("l2 converged\n");
				tot_pow(0,p,true);
				tot_pow(1,p,true);
				rate(0,0,p);
        			rate(0,1,p);
			        rate(1,0,p);
			        rate(1,1,p);*/
				return;
			}
		}
		last=pow;
	}

}


void optimise_p(struct isb_2_dual_params *p)
{

	int i,k=0;
	int b1,b2,b1max=0,b2max=0,b1max_last,b2max_last;
	int s0,s1,s0max=0,s1max=0,s0max_last,s1max_last;
	int services=2;
	double lk_max=0.0,lk=0.0;

	//printf("w1=%.8lf w2=%.8lf w3=%.8lf w4=%.8lf l1=%4.2lf l2=%4.2lf\n",p->w1,p->w2,p->w3,p->w4,p->l1,p->l2);

	for (i=0;i<DMTCHANNELS;i++) {
		//i=1;
		lk_max=0.0;
		lk=0.0;
		k=0;
		b1max=0;
		b2max=0;
		s0max=0;
		s1max=0;
		
		b1max_last=0;
                b2max_last=0;
                s0max_last=0;
                s1max_last=0;
		while (1) {
			k++;
			for (s0=0,lk=0.0,lk_max=0.0;s0<services;s0++) {
				for (b1=0;b1<=MAXBITSPERTONE;b1++) {
					//printf("b1=%d\tb2=%d\t",b1,b2);
					lk = l_k(b1,b2max,s0,s1max,i,p);
					//printf("lk=%lf\t",lk);
					//printf("p1=%lf\tp2=%lf\n",watts_to_dbmhz(p->p_mat[0][b1][b2][i])
					//			 ,watts_to_dbmhz(p->p_mat[1][b1][b2][i]));
					if (lk > lk_max) {	// we have found the max of l_k!
					//printf("new value of lk_max = %lf\n",lk_max);
						lk_max=lk;
						p->b_mat[0][i] = b1;
						//p->b_mat[1][i] = b2;
						p->s_mat[0][i] = s0;
						//p->s_mat[1][i] = s1;
						b1max=b1;
						//b2max=b2;
						s0max=s0;
						//s1max=s1;
					}
				}
			}
			//printf("iteration %d b1max = %d b2 = %d s0max = %d s1 = %d lk_max = %lf\n",k,b1max,b2max,s0max,s1max,lk_max);
			for (s1=0,lk=0.0,lk_max=0.0;s1<services;s1++) {
				for (b2=0;b2<=MAXBITSPERTONE;b2++) {
					//printf("b1=%d\tb2=%d\t",b1,b2);
					lk = l_k(b1max,b2,s0max,s1,i,p);
					//printf("lk=%lf\t",lk);
					//printf("p1=%lf\tp2=%lf\n",watts_to_dbmhz(p->p_mat[0][b1][b2][i])
					//			 ,watts_to_dbmhz(p->p_mat[1][b1][b2][i]));
					if (lk > lk_max) {	// we have found the max of l_k!
					//printf("new value of lk_max = %lf\n",lk_max);
						lk_max=lk;
						//p->b_mat[0][i] = b1;
						p->b_mat[1][i] = b2;
						//p->s_mat[0][i] = s0;
						p->s_mat[1][i] = s1;
						//b1max=b1;
						b2max=b2;
						//s0max=s0;
						s1max=s1;
					}
				}
			}
			//printf("iteration %d b1max = %d b2max = %d s0max = %d s1max = %d lk_max = %lf\n",k,b1max,b2max,s0max,s1max,lk_max);
			if ((b1max == b1max_last) && (b2max == b2max_last) && (s0max == s0max_last) && (s1max == s1max_last)) {
				//printf("yay! isb_dual converged!\n");
				//printf("iteration %d b1max = %d b2max = %d s0max = %d s1max = %d lk = %lf\n",k,b1max,b2max,s0max,s1max,lk_max);
				break;
			}
			b1max_last=b1max;
			b2max_last=b2max;
			s0max_last=s0max;
			s1max_last=s1max;
				
		}
		if (i==0);
			//printf("iteration %d b1max = %d b2 = %d s0max = %d s1 = %d lk = %lf\n",k,b1max,b2max,s0max,s1max,lk);
		/*printf("max value of lk on channel %d was %lf at b1 = %d b2 = %d p1 = %lf p2 = %lf\n",i,lk_max
						,p->b_mat[0][i],p->b_mat[1][i]
						,watts_to_dbmhz(p->p_mat[0][b1max][b2max][i]),watts_to_dbmhz(p->p_mat[1][b1max][b2max][i]));*/
		//exit(1);
	}
	//exit(1);
	
	//print_b_and_p(b,p);	
}	

void optimise_p2(struct isb_2_dual_params *p)
{

	int i,k=0;
	int b1,b2,b1max=0,b2max=0,b1max_last,b2max_last,_b1max,_b2max;
	int s0,s1,s0max=0,s1max=0,s0max_last,s1max_last;
	int services=2;
	double lk_max=0.0,lk=0.0;
	double lk1_max=0.0,lk1=0.0;

	//printf("w1=%.8lf w2=%.8lf w3=%.8lf w4=%.8lf l1=%4.2lf l2=%4.2lf\n",p->w1,p->w2,p->w3,p->w4,p->l1,p->l2);

	for (i=0;i<DMTCHANNELS;i++) {
		//i=1;		
		b1max=0;	// per channel bits that maximise lagrangian
		b2max=0;	// ditto
		s0max=0;	// per channel service that maximises lagrangian
		s1max=0;
		
		lk1_max=0.0;	// to find max lagragian per channel
		lk1=0.0;
		for (s0=0;s0<services;s0++) {
			for (s1=0;s1<services;s1++) {
				k=0;
				while (1) {
					k++;
					lk_max=0.0;	// to find max lagragian per service combo
					lk=0.0;					
					_b1max=0;	// per (service combo && channel) bits that maximise lagrangian
					_b2max=0;	// ditto
					for (b1=0;b1<=MAXBITSPERTONE;b1++) {
						lk = l_k(b1,_b2max,s0,s1,i,p);
						if (lk >lk_max) {
							lk_max=lk;
							_b1max=b1;
						}
					}
					for (b2=0;b2<=MAXBITSPERTONE;b2++) {
						lk=l_k(_b1max,b2,s0,s1,i,p);
						if (lk>lk_max) {
							lk_max=lk;
							_b2max=b2;
						}
					}
					if ((_b1max == b1max_last) && (_b2max == b2max_last)) {
						//printf("iterative bit maximisation loop converged after %d iterations\n",k);
						break;
					}
					b1max_last=_b1max;
					b2max_last=_b2max;
				}
				lk1=l_k(_b1max,_b2max,s0,s1,i,p);
				if (lk1>lk1_max) {
					lk1_max=lk1;
					s0max=s0,
					s1max=s1;
					b1max=_b1max;
					b2max=_b2max;
				}
			}
		}
		p->b_mat[0][i] = b1max;
		p->b_mat[1][i] = b2max;
		p->s_mat[0][i] = s0max;
		p->s_mat[1][i] = s1max;
		//printf("b1max=%d\tb2max=%d\ts0max=%d\ts1max=%d\tchannel=%d %lf\n",b1max,b2max,s0max,s1max,i,lk1_max);
		//exit(1);
	}
	//exit(1);
	//print_b_and_p(b,p);	

}

double l_k(int b1,int b2,int s0,int s1,int channel,struct isb_2_dual_params *p)
{

	int i;
	int b1s0=0;
	int b1s1=0;
	int b2s0=0;
	int b2s1=0;

	if (s0==0) {
		b1s0=b1;
		b1s1=0;
	}
	else if (s0==1) {
		b1s1=b1;
		b1s0=0;
	}
	if (s1==0) {
		b2s0=b2;
		b2s1=0;
	}
	else if (s1==1) {
		b2s1=b2;
		b2s0=0;
	}
		
	/*return p->w1*b1s0 + (1-p->w1)*b2s0 + p->w2*b1s1 + (1-p->w2)*b2s1 
		- p->l1*p->p_mat[0][s0][s1][b1][b2][channel]
		- p->l2*p->p_mat[1][s0][s1][b1][b2][channel];
        */

	return p->w3*b1s0 + p->w1*b1s1 + p->w4*b2s0 + p->w2*b2s1 
		- p->l1*p->p_mat[0][s0][s1][b1][b2][channel]
                - p->l2*p->p_mat[1][s0][s1][b1][b2][channel];

		
	//return p->w*b1 + (1-p->w)*b2 - p->l1*p->p_mat[0][b1][b2][channel] - p->l2*p->p_mat[1][b1][b2][channel];

	//return b1 + b2 - l1*p_mat[0][b1][b2][channel] - l2*p_mat[1][b1][b2][channel];
}

void isb_init_p_matrix(double P_MAT)
{

	int k,b1,b2;
	int s0,s1;
	int services=2;

	for (k=0;k<DMTCHANNELS;k++) {
		for (b1=0;b1<=MAXBITSPERTONE;b1++) {
			for (b2=0;b2<=MAXBITSPERTONE;b2++) {
				for (s0=0;s0<services;s0++) {
					for (s1=0;s1<services;s1++) {
						psd(p_mat,b1,b2,s0,s1,k);
					}
				}
			}
		}
	}


}

void psd(double P_MAT,int b1,int b2,int s0,int s1,int channel)
{

	int i=0;
	p_mat[0][s0][s1][b1][b2][channel] = 0.0;
	p_mat[1][s0][s1][b1][b2][channel] = 0.0;
	double p_last0;
	double p_last1;
	struct line* current = list_head;
	double gain0;
	double gain1;
	double diff=0.0,last_diff=0.0;
	int diverge=0;

	current=get_line(0);
	gain0 = current->gain[channel];
	current=get_line(1);
	gain1 = current->gain[channel];


	//printf("b1 = %d\tb2 = %d channel = %d\n",b1,b2,channel);
	//getchar();
	do {
		i++;
		last_diff=diff;
		p_last0=p_mat[0][s0][s1][b1][b2][channel];
		p_last1=p_mat[1][s0][s1][b1][b2][channel];
		//for (i=0;i<7;i++) {
		_psd(0,s0,s1,b1,b2,p_mat,channel,gain0);
		_psd(1,s0,s1,b1,b2,p_mat,channel,gain1);
		diff = p_mat[0][s0][s1][b1][b2][channel] - p_last0;
		if (diff > last_diff) {
			//printf("looks like this is diverging\n");
			diverge++;
			if (diverge == 3) {
				p_mat[0][s0][s1][b1][b2][channel] = DBL_MAX;
				p_mat[1][s0][s1][b1][b2][channel] = DBL_MAX;
				break;
			}
		}
		//printf("p1 = %lf\tp2 = %lf\n",watts_to_dbmhz(p[0][channel]),watts_to_dbmhz(p[1][channel]));
		//}
		//printf("p0 = %e p_last0 = %e p1 = %e p_last1 = %e\n",p[0][channel],p_last0,p[1][channel],p_last1);	
	} while ((fabs(p_last0-p_mat[0][s0][s1][b1][b2][channel]) > 1e-8) && (fabs(p_last1-p_mat[1][s0][s1][b1][b2][channel]) > 1e-8));// && i<15);
	//printf("p1 = %lf\tp2 = %lf\t",watts_to_dbmhz(p[0][channel]),watts_to_dbmhz(p[1][channel]));
	//printf("took %d iterations to converge\n",i);
	//if (diverge==3)
		//printf("psd solution diverged for b1=%d b2=%d channel =%d\n",b1,b2,channel);
}

void _psd(int line_id,int s0,int s1,int b1,int b2,double P_MAT,int channel,double gain)
{
	double gamma_hat;
	int b;

	if (line_id == 0) {	// we're finding the psd on line0
		b=b1;
		if (s0==0)	// the service on line0 is service 0, i.e. low BER
			gamma_hat=pow(10,(GAMMA0+MARGIN)/10);
		else 
			gamma_hat=pow(10,(GAMMA1+MARGIN)/10);
	}
	else { 
		b=b2;
		if (s1==0)	// the service on line1 is service 0
			gamma_hat=pow(10,(GAMMA0+MARGIN)/10);
		else
			gamma_hat=pow(10,(GAMMA1+MARGIN)/10);
	}
	p_mat[line_id][s0][s1][b1][b2][channel] = (pow(2,(double)b)-1) * gamma_hat/_cnr(line_id,p_mat,channel,gain,b1,b2,s0,s1);

}

double _cnr(int line_id,double P_MAT,int channel,double gain,int b1,int b2,int s0,int s1)
{
	//struct line *current = list_head;
	double noise;
	static int calls=0;
	int i,k;
	static double alien_xtalk_array[2][DMTCHANNELS];
	static double bk_n;

	if (calls++ == 0) {
		bk_n = dbmhz_to_watts(-140);
		for (i=0;i<lines;i++) {
			for (k=0;k<DMTCHANNELS;k++) {
				alien_xtalk_array[i][k] = alien_xtalk(i,k);
			}
		}
	}
		
	
	noise = fsan_sum(osb_fext(line_id,p_mat,channel,b1,b2,s0,s1),alien_xtalk_array[line_id][channel]) + bk_n;
	//current=get_line(line_id);

	return gain/noise;	

}

double osb_fext(int line_id,double P_MAT,int channel,int b1,int b2,int s0,int s1)
{

	int i;
	double xtalk_gain;
	double noise=0.0;
	
	for (i=0;i<lines;i++) {
		if (i != line_id) {
			xtalk_gain = *(channel_matrix + channel + (i * DMTCHANNELS) + (lines * line_id * DMTCHANNELS));
			noise += xtalk_gain * p_mat[i][s0][s1][b1][b2][channel];
		}
	}	

	return noise;
}

int rate(int line_id,int service,struct isb_2_dual_params *p,bool print)
{

	int rate=0,i;

	for (i=0;i<DMTCHANNELS;i++) {
		if (p->s_mat[line_id][i] == service)
			rate += p->b_mat[line_id][i];
	}

	if (print)
		printf("current rate on line %d service %d is %d\n",line_id,service,rate);

	return rate;

}

double tot_pow(int line_id,struct isb_2_dual_params *p,bool print)
{
	
	double tot=0.0;
	int i;
	
	for(i=0;i<DMTCHANNELS;i++) {
		tot += p->p_mat[line_id][p->s_mat[0][i]][p->s_mat[1][i]][p->b_mat[0][i]][p->b_mat[1][i]][i];
	}

	if (print)
		printf("current power on line %d is %lf\n",line_id,tot);

	return tot;
	
}

void isb_dual_init_lines(struct isb_2_dual_params *p)
{

	int i,k;
	struct line* current=list_head;

	for (i=0;i<lines;i++) {
		current = get_line(i);	
		current->is_dual = 1;	// dual qos line
		current->gamma[0] = GAMMA0;
		current->gamma[1] = GAMMA1;

		for (k=0;k<DMTCHANNELS;k++) {
			current->service[k] = p->s_mat[i][k];
			current->b[k] = p->b_mat[i][k];
			current->psd[k] = watts_to_dbmhz(p->p_mat[i][p->s_mat[0][k]][p->s_mat[1][k]][p->b_mat[0][k]][p->b_mat[1][k]][k]);
		}
	}
}
