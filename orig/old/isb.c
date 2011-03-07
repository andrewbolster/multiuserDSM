#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "multiuser_load.h"
#include <string.h>
#include <float.h>
#include <stdbool.h>
#define P_MAT p_mat[2][MAXBITSPERTONE+1][MAXBITSPERTONE+1][DMTCHANNELS]
#define B_MAT b_mat[2][DMTCHANNELS]

static void optimise_w(struct isb_2_params *p);
static void optimise_l1(struct isb_2_params *p);
static void optimise_l2(struct isb_2_params *p);
static void isb_init_params(struct isb_2_params *p);
static void optimise_p(struct isb_2_params *p);
static double l_k(int b1,int b2,int channel,struct isb_2_params *p);
static void osb_init_p_matrix(double P_MAT);
static void psd(double P_MAT,int b1,int b2,int channel);
static void _psd(int line_id,int b1,int b2,double P_MAT,int channel,double gain);
static double _cnr(int line_id,double P_MAT,int channel,double gain,int b1,int b2);
static double osb_fext(int line_id,double P_MAT,int channel,int b1,int b2);
static int rate(int line_id,int B_MAT);
static double tot_pow(int line_id,struct isb_2_params *p);

struct isb_2_params {

	double p_mat[2][MAXBITSPERTONE+1][MAXBITSPERTONE+1][DMTCHANNELS];
	int b_mat[2][DMTCHANNELS];
	double w;
	double w_max;
	double w_min;
	double l1;
	double l1_min;
	double l1_max;
	double l2;
	double l2_min;
	double l2_max;
	int rate0_target;
	double p0_budget;
	double p1_budget;
};

struct isb_2_params *isb_2()
{

	struct isb_2_params *p;



	p = (struct isb_2_params *)malloc(sizeof(struct isb_2_params));

	if (p==NULL) {
		printf("malloc sucks\n");
		exit(1);
	}

	
	isb_init_params(p);

	optimise_w(p);

	rate(0,p->b_mat);
	rate(1,p->b_mat);

	return p;
}


void optimise_w(struct isb_2_params *p)
{

	int e = 5;

	while(abs(rate(0,p->b_mat) - p->rate0_target) > e) {
		p->w=(p->w_max+p->w_min)/2;
		printf("w=%lf\n",p->w);

		optimise_l1(p);

		if (rate(0,p->b_mat) > p->rate0_target) 
			p->w_max=p->w;
		else
			p->w_min=p->w; 

	}

	//optimise_l1(p);
}

/*
void optimise_w(struct osb_2_params *p)
{

	int e = 5;
	int w1;
	int w2;
	double w_grad;
	double delta=0.01;

	p->w=0.5;
	optimise_l1(p);
	w1=rate(0,p->b_mat);
	printf("first point = %d\n",w1);

	while(abs(w1 - p->rate0_target) > e) {
		//p->w=(p->w_max+p->w_min)/2;
		//printf("w=%lf\n",p->w);
		p->w+=delta;
		optimise_l1(p);
		w2=rate(0,p->b_mat);
		printf("second point = %d\n",w2);
		p->w-=delta;

		w_grad=((double)w2-(double)w1)/delta;
		printf("gradient = %lf\n",w_grad);
		
		p->w=p->w-(w1-p->rate0_target)/w_grad;
		printf("next guess = %lf\n",p->w);
		optimise_l1(p);
		w1=rate(0,p->b_mat); 

	}

}
*/
void optimise_l1(struct isb_2_params *p)
{

	double pow,last;

	p->l1=1.0;
	p->l1_min=0.0;
	do {
		p->l1*=2;
		optimise_l2(p);
	} while(tot_pow(0,p) > p->p0_budget);
	p->l1_max=p->l1;
	//printf("First value of l1 to meet power constraint = %lf\n",p->l1_max);

	while(1) {
		p->l1=(p->l1_max+p->l1_min)/2;
		
		optimise_l2(p);

		pow=tot_pow(0,p);

		if (pow > p->p0_budget) {
			p->l1_min=p->l1;
		}
		else {
			p->l1_max=p->l1;
			if (pow == last)
			//if (fabs(pow-last) < 1e-2)
				return;
		}
		last=pow;
	}
}


void optimise_l2(struct isb_2_params *p)
{

	double pow,last;
	int calls=0;

	p->l2=1.0;
	p->l2_min=0.0;
	do {
		p->l2*=2;
		optimise_p(p);
		calls++;
	} while(tot_pow(1,p) > p->p1_budget);
	p->l2_max=p->l2;	
	//printf("First value of l2 to meet power constraint = %lf\n",p->l2_max);

	while(1) {
		p->l2=(p->l2_max+p->l2_min)/2;
		optimise_p(p);
		calls++;

		pow = tot_pow(1,p);	

		if (pow > p->p1_budget) {
			p->l2_min=p->l2;
		}
		else {
			p->l2_max=p->l2;
			if (pow==last) {
			//if (fabs(pow-last) < 1e-2)
				//printf("Root found in %d calls to optimise_p\n",calls);
				//exit(0);
				return;
			}
		}
		last=pow;
	}

}


/*
void optimise_l2(struct osb_2_params *p)
{

	double pow,last;
	int calls=0;
	double l2_last;

	p->l2=1.0;

	while (1) {
		optimise_p(p);
		pow=tot_pow(1,p);
		printf("power = %lf\n",pow);
		p->l2=p->l2 + (pow - p->p1_budget)/(p->l2-l2_last);

		if (pow==last && (pow < p->p1_budget))
			break;

		last=pow;
	}

}
*/
/*
void optimise_l2(struct osb_2_params *p)
{

	double pow,last;
	double l2_g1;
	double l2_g2;
	double l2_grad;
	double delta=1;
	int calls=0;
	bool delta_changed=false;
	bool debug=true;
	FILE *fp;
	p->l2=1.0;

	
	while(1) {
		optimise_p(p);
		calls++;				
		l2_g1 = tot_pow(1,p);			// find first point on pow curve
		if (debug)
			printf("first point on p curve, l2=%lf\tp=%lf\n",p->l2,l2_g1);

	
		p->l2 += delta;				// add delta to l2
		optimise_p(p);
		calls++;				
		l2_g2 = tot_pow(1,p);			// find second point on pow curve	
		if (debug)
			printf("second point on p curve, l2=%lf\tp=%lf\n",p->l2,l2_g2);
		p->l2 -= delta;				// remove delta from l2
		
		l2_grad=(l2_g2-l2_g1)/delta;		// gradient of pow with respect to l2

//		if (l2_grad == 0) {
//			delta_changed=true;
//			delta*=2;
//			printf("l2_grad was 0, so increasing delta to %lf\n",delta);
//			continue;
//		}
//
//		if (delta_changed) {
//			delta=0.1;
//			delta_changed=false;	
//		}

		if (debug)
			printf("gradient = %lf\n",l2_grad);

		p->l2 = p->l2 - (l2_g1-p->p1_budget)/l2_grad;		// next guess at root x_n+1=x_n-f(x_n)/f'(x_n)
		if (debug)
			printf("next guess at root = %lf\n",p->l2);

		optimise_p(p);
		calls++;				
		pow = tot_pow(1,p);			
		if (debug)
			printf("Value of p at next guess = %lf\n",pow);		
		
		if (pow==last) {
		//if (fabs(pow-last) < 1e-2)
			printf("Root found in %d calls to optimise_p\n",calls);
			printf("w=%4.2lf l1=%4.2lf l2=%4.2lf\n\n\n",p->w,p->l1,p->l2);
			exit(0);
			return;
		}
		last=pow;
	}


//	fp=fopen("f_l2.txt","w");
//	
//	if (fp==NULL)
//		exit(20);
//
//	for (p->l2=1;p->l2<200;p->l2+=0.1) {
//		optimise_p(p);
//		pow=tot_pow(1,p);
//
//		fprintf(fp,"%lf\t%lf\n",p->l2,p->p0_budget-pow);
//	}
//
//	fclose(fp);
//
//	exit(0);

}
*/
void optimise_p(struct isb_2_params *p)
{

	int i,k=0;
	int b1=0,b2=0,b1max=0,b2max=0,b1max_last,b2max_last;
	double lk_max=0.0,lk=0.0;

	//printf("w=%4.2lf l1=%4.2lf l2=%4.2lf\n",p->w,p->l1,p->l2);

	for (i=0;i<DMTCHANNELS;i++) {
		//i=1;
		lk_max=0.0;
		lk=0.0;
	 	while(1) {	
			for (b1=0,lk=0.0,lk_max=0.0;b1<=MAXBITSPERTONE;b1++) {
				//printf("b1=%d\tb2=%d\t",b1,b2);
				lk = l_k(b1,b2max,i,p);
				//printf("lk=%lf\t",lk);
				//printf("p1=%lf\tp2=%lf\n",watts_to_dbmhz(p->p_mat[0][b1][b2][i])
				//			 ,watts_to_dbmhz(p->p_mat[1][b1][b2][i]));
				if (lk > lk_max) {	// we have found the max of l_k!
					//printf("new value of lk_max = %lf\n",lk_max);
					lk_max=lk;
					p->b_mat[0][i] = b1;
					b1max=b1;
				}
			}
			//printf("iteration %d b1max = %d b2 = %d lk=%lf\n",k,b1max,b2max,lk_max);
			for (b2=0,lk=0.0,lk_max=0.0;b2<=MAXBITSPERTONE;b2++) {
				lk = l_k(b1max,b2,i,p);
				if (lk > lk_max) {
					lk_max=lk;
					p->b_mat[1][i] = b2;
					b2max=b2;
				}
			}
			//printf("iteration %d b1 = %d b2max = %d lk=%lf\n",k,b1max,b2max,lk_max);
			k++;
			if ((b2max == b2max_last) && (b1max == b1max_last)) {
				//printf("isb converged hooray! b1max=%d b2max=%d lk_max=%lf\n",b1max,b2max,lk_max);
				break;
			}
			b1max_last=b1max;
			b2max_last=b2max;
	
		}
		//exit(1);
		/*printf("channel=%d lkmax= %lf b1 = %d b2 = %d p1 = %lf p2 = %lf\n",i,lk_max
						,p->b_mat[0][i],p->b_mat[1][i]
						,watts_to_dbmhz(p->p_mat[0][b1max][b2max][i]),watts_to_dbmhz(p->p_mat[1][b1max][b2max][i]));*/
	}
	
	//print_b_and_p(b,p);	
}	

double l_k(int b1,int b2,int channel,struct isb_2_params *p)
{

	int i;
		
	return p->w*b1 + (1-p->w)*b2 - p->l1*p->p_mat[0][b1][b2][channel] - p->l2*p->p_mat[1][b1][b2][channel];

	//return 1.9*b1 + 0.1*b2 - p->l1*p->p_mat[0][b1][b2][channel] - p->l2*p->p_mat[1][b1][b2][channel];
}

void isb_init_params(struct isb_2_params *p)
{

	memset(p->b_mat,0,sizeof(int)*lines*DMTCHANNELS);
	memset(p->p_mat,0,sizeof(double)*lines*(MAXBITSPERTONE+1)*(MAXBITSPERTONE+1)*DMTCHANNELS);

	p->w_min=0.0;
	p->w_max=1.0;

	p->l1_min=0.0;
	p->l1_max=1.0;

	p->l2_min=0.0;
	p->l2_max=1.0;

	p->rate0_target=800;

	p->p0_budget=0.1;
	p->p1_budget=0.1;

	osb_init_p_matrix(p->p_mat);

}

void osb_init_p_matrix(double P_MAT)
{

	int k,b1,b2;

	for (k=0;k<DMTCHANNELS;k++) {
		for (b1=0;b1<=MAXBITSPERTONE;b1++) {
			for (b2=0;b2<=MAXBITSPERTONE;b2++) {
				psd(p_mat,b1,b2,k);
			}
		}
	}


}

void psd(double P_MAT,int b1,int b2,int channel)
{

	int i=0;
	p_mat[0][b1][b2][channel] = 0.0;
	p_mat[1][b1][b2][channel] = 0.0;
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
		p_last0=p_mat[0][b1][b2][channel];
		p_last1=p_mat[1][b1][b2][channel];
		//for (i=0;i<7;i++) {
		_psd(0,b1,b2,p_mat,channel,gain0);
		_psd(1,b1,b2,p_mat,channel,gain1);
		diff = p_mat[0][b1][b2][channel] - p_last0;
		if (diff > last_diff) {
			//printf("looks like this is diverging\n");
			diverge++;
			if (diverge == 3) {
				p_mat[0][b1][b2][channel] = DBL_MAX;
				p_mat[1][b1][b2][channel] = DBL_MAX;
				break;
			}
		}
		//printf("p1 = %lf\tp2 = %lf\n",watts_to_dbmhz(p[0][channel]),watts_to_dbmhz(p[1][channel]));
		//}
		//printf("p0 = %e p_last0 = %e p1 = %e p_last1 = %e\n",p[0][channel],p_last0,p[1][channel],p_last1);	
	} while ((fabs(p_last0-p_mat[0][b1][b2][channel]) > 1e-8) && (fabs(p_last1-p_mat[1][b1][b2][channel]) > 1e-8));// && i<15);
	//printf("p1 = %lf\tp2 = %lf\t",watts_to_dbmhz(p[0][channel]),watts_to_dbmhz(p[1][channel]));
	//printf("took %d iterations to converge\n",i);
	//if (diverge==3)
		//printf("psd solution diverged for b1=%d b2=%d channel =%d\n",b1,b2,channel);
}

void _psd(int line_id,int b1,int b2,double P_MAT,int channel,double gain)
{
	double gamma_hat=pow(10,(GAMMA+MARGIN)/10);
	int b;

	if (line_id == 0)
		b=b1;
	else 
		b=b2;
	
	p_mat[line_id][b1][b2][channel] = (pow(2,(double)b)-1) * gamma_hat/_cnr(line_id,p_mat,channel,gain,b1,b2);

}

double _cnr(int line_id,double P_MAT,int channel,double gain,int b1,int b2)
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
		
	
	noise = fsan_sum(osb_fext(line_id,p_mat,channel,b1,b2),alien_xtalk_array[line_id][channel]) + bk_n;
	//current=get_line(line_id);

	return gain/noise;	

}

double osb_fext(int line_id,double P_MAT,int channel,int b1,int b2)
{

	int i;
	double xtalk_gain;
	double noise=0.0;
	
	for (i=0;i<lines;i++) {
		if (i != line_id) {
			xtalk_gain = *(channel_matrix + channel + (i * DMTCHANNELS) + (lines * line_id * DMTCHANNELS));
			noise += xtalk_gain * p_mat[i][b1][b2][channel];
		}
	}	

	return noise;
}

int rate(int line_id,int B_MAT)
{

	int rate=0,i;

	for (i=0;i<DMTCHANNELS;i++) {
		rate += b_mat[line_id][i];
	}

	printf("current rate on line %d is %d\n",line_id,rate);

	return rate;

}

double tot_pow(int line_id,struct isb_2_params *p)
{
	
	double tot=0.0;
	int i;
	
	for(i=0;i<DMTCHANNELS;i++) {
		tot += p->p_mat[line_id][p->b_mat[0][i]][p->b_mat[1][i]][i];
	}
	//printf("current power on line %d is %lf\n",line_id,tot);

	return tot;
	
}
