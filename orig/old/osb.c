#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "multiuser_load.h"
#include <string.h>
#include <float.h>
#include <stdbool.h>
#include "osb_2.h"
#include "psd_vector.h"
//#define P_MAT p_mat[2][MAXBITSPERTONE+1][MAXBITSPERTONE+1][DMTCHANNELS]
//#define B_MAT b_mat[2][DMTCHANNELS]
static void optimise_w(struct osb_2_params *p);
static void optimise_l1(struct osb_2_params *p);
static void optimise_l2(struct osb_2_params *p);
static void optimise_p(struct osb_2_params *p);
static double l_k(int b1,int b2,int channel,struct osb_2_params *p);
static void osb_init_p_matrix(double P_MAT);
//static void psd(double P_MAT,int b1,int b2,int channel);
//static void _psd(int line_id,int b1,int b2,double P_MAT,int channel,double gain);
//static double _cnr(int line_id,double P_MAT,int channel,double gain,int b1,int b2);
//static double osb_fext(int line_id,double P_MAT,int channel,int b1,int b2);
static int rate(int line_id,int B_MAT);
static double tot_pow(int line_id,struct osb_2_params *p);
//static void greedy_lk_search(int,struct osb_2_params *,int * b_max,bool);
static bool p_debug=true;
static psd_vector *psder;
struct osb_2_params *osb_2(struct osb_2_params *p)
{
	if (p==NULL) {
		printf("malloc sucks\n");
		exit(1);
	}
	psder = new psd_vector;
	optimise_w(p);
	rate(0,p->b_mat);
	rate(1,p->b_mat);
	p->osb_init_lines();
	return p;
}
/*
struct rate_pair osb_rate_region(double w)
{
	struct rate_pair ret;
	struct osb_2_params *p;
	p = (struct osb_2_params*)malloc(sizeof(struct osb_2_params));
	if (p==NULL) {
		printf("malloc sucks\n");
		exit(1);
	}
	osb_init_params(p);
	p->w=w;
	optimise_l1(p);
	ret.r1 = rate(0,p->b_mat);
	ret.r2 = rate(1,p->b_mat);
	free(p);
	return ret;
}
*/
void optimise_w(struct osb_2_params *p)
{
	if (p->rate1_target == NOT_SET) {
		while(abs(rate(0,p->b_mat) - p->rate0_target) > p->e) {
			p->w=(p->w_max+p->w_min)/2;
			//printf("w=%lf\n",p->w);
			optimise_l1(p);
			if (rate(0,p->b_mat) > p->rate0_target) 
				p->w_max=p->w;
			else
				p->w_min=p->w;
		}
	}
	else if (p->rate0_target == NOT_SET) {	
		while(abs(rate(1,p->b_mat) - p->rate1_target) > p->e) {
			p->w=(p->w_max+p->w_min)/2;
			optimise_l1(p);
			if (rate(1,p->b_mat) > p->rate1_target)
                                p->w_min=p->w;
                        else
                                p->w_max=p->w;
		}
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
void optimise_l1(struct osb_2_params *p)
{
	double pow,last=DBL_MAX;	// compiler warning
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
void optimise_l2(struct osb_2_params *p)
{
	double pow,last=DBL_MAX;	// compiler warning
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
void optimise_p(struct osb_2_params *p)
{
	int i;
	int b1,b2,b1max,b2max;
	double lk_max=0.0,lk=0.0;
	int *b_max = new int[lines];
	if (p_debug) 
		printf("w=%20.18lf l1=%4.2lf l2=%4.2lf\n",p->w,p->l1,p->l2);
	for (i=0;i<DMTCHANNELS;i++) {
		//i=1;
		lk_max=0.0;
		lk=0.0;
		for (b1=0;b1<=MAXBITSPERTONE;b1++) {
			//b1=13;
			for (b2=0;b2<=MAXBITSPERTONE;b2++) {
			//b2=15;
				//printf("%d\t%d\t",b1,b2);
				lk = l_k(b1,b2,i,p);
				//printf("%lf\n",lk);
				//printf("p1=%lf\tp2=%lf\n",watts_to_dbmhz(p->p_mat[0][b1][b2][i])
				//			 ,watts_to_dbmhz(p->p_mat[1][b1][b2][i]));
				if (lk > lk_max) {	// we have found the max of l_k!
					//printf("new value of lk_max = %lf\n",lk_max);
					lk_max=lk;
					p->b_mat[0][i] = b1;
					p->b_mat[1][i] = b2;
					b1max=b1;
					b2max=b2;
				}
			}
		}
		//greedy_lk_search(i,p,b_max,false);
		//printf("Max value of l(k) was at %d %d\n",b1max,b2max);
	/*	
		if (b1max != b_max[0] || b2max != b_max[1]) {
			printf("greedy search got it wrong tone %d\n",i);
			print_vector(b_max,"b_max");
			printf("b1max = %d b2max = %d\n",b1max,b2max);
			printf("l1 = %lf l2 =%lf\n",p->l1,p->l2);
			getchar();
			greedy_lk_search(i,p,b_max,true);
		}
	*/	
		p->b_mat[0][i] = b_max[0];
		p->b_mat[1][i] = b_max[1];
		/*printf("channel=%d lk_max= %lf at b1 = %d b2 = %d p1 = %lf p2 = %lf\n",i,lk_max
						,p->b_mat[0][i],p->b_mat[1][i]
						,watts_to_dbmhz(p->p_mat[0][b1max][b2max][i]),watts_to_dbmhz(p->p_mat[1][b1max][b2max][i]));*/
	}
	//exit(1);
	//print_b_and_p(b,p);	
}	
#ifdef NOTDEF
void greedy_lk_search(int tone,struct osb_2_params *o,int *max_b,bool debug)
{
	double cost[lines];
	double min_p[lines];
	double min_cost=DBL_MAX;
	int min_user;
	double max_lk=-DBL_MAX;
	double lk;
	double *p = new double [lines];
	double *old_p = new double [lines];
	double *g = new double [lines];
	int *b = new int [lines];
	double *b_d = new double [lines];
	//int *max_b = new int[lines];
	for (int user=0;user<lines;user++) {
		b[user] = 0;			// start at zero bit vector
		b_d[user] = 0;
		g[user] = line_array[user]->gamma[line_array[user]->service[tone]];
	}	
/*
	while (1) {
		min_cost=DBL_MAX;
		//calculate_psd_vector(b,g,tone,p);
		psder->calc(b,g,tone,p);
		if (debug) {
			print_vector(b,"old b");
			print_vector(p,"old p");
		}
		memcpy(old_p,p,sizeof(double)*lines);
		for (int user=0;user<lines;user++) {
			if (b[user] == MAXBITSPERTONE) {
				if (debug)
					printf("line %d is at max bits\n",user);
				cost[user] = DBL_MAX;
				continue;
			}
			cost[user]=0;
			b[user]++;
			//calculate_psd_vector(b,g,tone,p);
			psder->calc(b,g,tone,p);
			if (debug) {
				print_vector(b,"new b");
				print_vector(p,"new_p");			
			}
			for(int user1=0;user1<lines;user1++) {
				if (p[user1] < 0) {
					cost[user] = DBL_MAX;
					break;
				}
				//cost[user]+=(p[user1]-old_p[user1]);
			}
			cost[user]+=o->l1*(p[0]-old_p[0]);
			cost[user]+=o->l2*(p[1]-old_p[1]);
			if (debug) {
				printf("cost[%d] = %g\n",user,cost[user]);
				getchar();
			}
			if (cost[user] < min_cost) {
				min_cost=cost[user];
				min_user=user;	
				memcpy(min_p,p,sizeof(double)*lines);
			}	
			b[user]--;
		}
		if (debug)
			printf("min cost was on user %d\n",min_user);
		b[min_user]++;
		lk = o->w*b[0] + (1-o->w)*b[1] - o->l1*min_p[0] - o->l2*min_p[1];
		if (debug) {
			printf("value of lk having increased b[%d] = %lf\n",min_user,lk);
			print_vector(b,"b_vec");
		}
		if (lk > max_lk) {
			max_lk=lk;
		//	printf("still going up\n");
			memcpy(max_b,b,sizeof(int)*lines);
		}
		else {
		//	printf("last value was lower\n");
			break;
		}
	}
	*/
	while (1) {
		min_cost=DBL_MAX;
		/*
		//calculate_psd_vector(b,g,tone,p);
		psder->calc(b,g,tone,p);	// calculate current p vector given current b vector
		//
		if (debug) {
			print_vector(b,"old b");
			print_vector(p,"old p");
		}
		memcpy(old_p,p,sizeof(double)*lines);	// save this value in old_p
		*/
		for (int user=0;user<lines;user++) {
			// set all p values to reference value
			for (int user1=0;user1<lines;user1++) {
				p[user1] = minus_40_dbmhz_watts;
				b_d[user] = 0;				
			}
			// if b[user] == 15 cant add any more bits, so set it's cost to DBL_MAX
			if (b[user] == MAXBITSPERTONE) {						
				if (debug)
					printf("line %d is at max bits\n",user);
				cost[user] = DBL_MAX;
				continue;
			}
			b_d[user]=(double)b[user];
			calculate_b_vector_one_bit_set(p,9.95,tone,user,b_d);
			memcpy(old_p,p,sizeof(double)*lines);	// save this value in old_p
			if (debug) {		
				print_vector(b_d,"old b");
	                        print_vector(p,"old p");
			}
			// initialise cost to zero, and add a bit to user
			cost[user]=0;
			//calculate_psd_vector(b,g,tone,p);
			//psder->calc(b,g,tone,p);
			b_d[user]++;
			calculate_psd_vector(b_d,g,tone,p);
			if (debug) {
				print_vector(b_d,"b");
	                        print_vector(p,"p");
			}
			for(int user1=0;user1<lines;user1++) {
				if (p[user1] < 0) {
					cost[user] = DBL_MAX;
					break;
				}
				//cost[user]+=(p[user1]-old_p[user1]);
			}
			cost[user]+=(p[0]-old_p[0]);
			cost[user]+=(p[1]-old_p[1]);
			if (debug) {
				printf("cost[%d] = %g\n",user,cost[user]);
				getchar();
			}
			if (cost[user] < min_cost) {
				min_cost=cost[user];
				min_user=user;	
				//memcpy(min_p,p,sizeof(double)*lines);
			}	
		}
		if (debug)
			printf("min cost was on user %d\n",min_user);
		b[min_user]++;
		psder->calc(b,g,tone,p);
		lk = o->w*b[0] + (1-o->w)*b[1] - o->l1*p[0] - o->l2*p[1];
		if (debug) {
			printf("value of lk having increased b[%d] = %lf\n",min_user,lk);
			print_vector(b,"b_vec");
		}
		if (lk > max_lk) {
			max_lk=lk;
			if (debug)
				printf("still going up\n");
			memcpy(max_b,b,sizeof(int)*lines);
		}
		else {
			if (debug)
				printf("last value was lower\n");
			break;
		}
	}
}
#endif
double l_k(int b1,int b2,int channel,struct osb_2_params *p)
{
	//int i;
	return p->w*b1 + (1-p->w)*b2 - p->l1*p->p_mat[0][b1][b2][channel] - p->l2*p->p_mat[1][b1][b2][channel];
	//return 1.9*b1 + 0.1*b2 - p->l1*p->p_mat[0][b1][b2][channel] - p->l2*p->p_mat[1][b1][b2][channel];
}
void osb_2_params::osb_2_params_init()
{
	//memset(b_mat,0,sizeof(int)*lines*DMTCHANNELS);
	//memset(p_mat,0,sizeof(double)*lines*(MAXBITSPERTONE+1)*(MAXBITSPERTONE+1)*DMTCHANNELS);
	w_min=0.0;
	w_max=1.0;
	l1_min=0.0;
	l1_max=1.0;
	l2_min=0.0;
	l2_max=1.0;
	rate0_target=NOT_SET;
	rate1_target=NOT_SET;
	p0_budget=0.110;
	p1_budget=0.110;
	e=NOT_SET;
}
osb_2_params::osb_2_params()
{
	osb_2_params_init();
	osb_init_p_matrix(p_mat);
}
void osb_init_p_matrix(double P_MAT)
{
	int k,b1,b2;
	int *b = new int [lines];
	double *g = new double [lines];
	double *p = new double [lines];	
	g[0]=GAMMA;
	g[1]=GAMMA;
	for (k=0;k<DMTCHANNELS;k++) {
		for (b1=0;b1<=MAXBITSPERTONE;b1++) {
			for (b2=0;b2<=MAXBITSPERTONE;b2++) {
				b[0]=b1;
				b[1]=b2;
				calculate_psd_vector(b,g,k,p);
				p_mat[0][b1][b2][k]=p[0];	
				p_mat[1][b1][b2][k]=p[1];
				if (p[0] < 0 || p[1] < 0) {
					p_mat[0][b1][b2][k]=DBL_MAX;
                                	p_mat[1][b1][b2][k]=DBL_MAX;	
				}
			}
		}
	}
}
/*
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
*/
/*
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
*/
/*
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
	noise = osb_fext(line_id,p_mat,channel,b1,b2) + alien_xtalk_array[line_id][channel] + bk_n;
	//current=get_line(line_id);
	return gain/noise;	
}
*/
/*
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
*/
int rate(int line_id,int B_MAT)
{
	int rate=0,i;
	for (i=0;i<DMTCHANNELS;i++) {
		rate += b_mat[line_id][i];
	}
	if (p_debug)
		printf("current rate on line %d is %d\n",line_id,rate);
	return rate;
}
void osb_2_params::osb_init_lines()
{
	int i,k;
        struct line* current=list_head;
        for (i=0;i<lines;i++) {
                current = get_line(i);
                current->is_dual = 0;   // dual qos line
                strcpy(current->loading_algo,"OSB");
		//current->alien_xtalk_model=ALIEN_XTALK_MODEL;         //FIXME
                for (k=0;k<DMTCHANNELS;k++) {
                        current->b[k] = b_mat[i][k];
                        current->psd[k] = watts_to_dbmhz(p_mat[i][b_mat[0][k]][b_mat[1][k]][k]);
                }
        }
}
double tot_pow(int line_id,struct osb_2_params *p)
{
	double tot=0.0;
	int i;
	for(i=0;i<DMTCHANNELS;i++) {
		tot += p->p_mat[line_id][p->b_mat[0][i]][p->b_mat[1][i]][i];
	}
	//printf("current power on line %d is %lf\n",line_id,tot);
	return tot;
}
