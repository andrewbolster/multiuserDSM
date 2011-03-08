#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "multiuser_load.h"
static void sort(double *,int *,int *);
static int min_index(double *,int);
static double am_p_min_fm(double *p,double *cnr,int *b,int n,int rate,double gamma,int half);
static double dbw_load(double *p,double *cnrs,int *b,int slice,int rate1,int rate2,double gamma1,double gamma2);
#define TOP_TONES 0	// higher snr service
#define BOTTOM_TONES 1	// lower snr service
#define HIGH 1
#define LOW 0
void dbw_cbr_exhaustive(int line_id,double gamma1,double gamma2, int rate1, int rate2)
{
	struct line *current;
	int i,j,min_slice;
	double cnrs[DMTCHANNELS];
	int _sort[DMTCHANNELS];
	int _desort[DMTCHANNELS];
	double p[DMTCHANNELS];		// local array to power energy results from minimisation routine
	double min_p[DMTCHANNELS];
	double p1,p2;
	double p_tot;
	double p_min_tot = 1e100;	// this is the minimum power total for all channels, set to large value to start
	int b[DMTCHANNELS];
	int min_b[DMTCHANNELS];
	int s0=0,s1=0;
	current = get_line(line_id);
	current->is_dual = 1;           // remember that this line is dual qos!
        current->gamma[0] = gamma1;     // and remember gamma values! - higher gamma i.e lower BER
        current->gamma[1] = gamma2;     // lower gamma i.e. higher BER
        strcpy(current->loading_algo,"dbw_cbr_exhaustive");
	memset(p,0,(sizeof(double) * DMTCHANNELS));
	sort(current->cnr,_sort,_desort);	// create sort and desort arrays for cnrs for line
	for (i=0;i<DMTCHANNELS;i++) {
		cnrs[i] = current->cnr[_sort[i]];	// copy sorted cnrs to local array
	}
	for (i=1;i<DMTCHANNELS;i++) {
                if (i*MAXBITSPERTONE < rate1) {
                        //printf("%d channels cannot support rate of %d\n",i,rate1);
                        continue;
                }
                if ((DMTCHANNELS-i)*MAXBITSPERTONE < rate2) {
                        //printf("%d channels cannot support rate of %d\n",(DMTCHANNELS-i),rate1);
                        continue;
                }
                p1 = am_p_min_fm(p,cnrs,b,i,rate1,gamma1,TOP_HALF);		//top half of channels, i.e. highest SNRs
		p2 = am_p_min_fm(&p[i],&cnrs[i],&b[i],(DMTCHANNELS-i),rate2,gamma2,BOTTOM_HALF);
		//printf("slice = %d\tpower = %lf\n",i,p1+p2);
		for (j=0,p_tot=0.0;j<DMTCHANNELS;j++) {
                        p_tot += p[j];
                }
		//printf("p_tot = %lf\ti = %d\n",p_tot,i);
		if (p_tot < p_min_tot) {
                        p_min_tot = p_tot;
                        min_slice = i;
                        memcpy(min_p,p,(sizeof(double) * DMTCHANNELS));
                        memcpy(min_b,b,(sizeof(double) * DMTCHANNELS));
                        //printf("New p_min, i = %d\n",i);
                }
	}
	printf("p_min = %lf\tslice point = %d\n",p_min_tot,min_slice);
/*
	for (i=0;i<DMTCHANNELS;i++) {
                //printf("b[%d] = %d\tp[%d] = %4.2e\n",i,min_b[i],i,min_p[i]);
                current->b[i] = min_b[_desort[i]];
                current->psd[i] = watts_to_dbmhz(min_p[_desort[i]]);
        }
*/
	for (i=0;i<DMTCHANNELS;i++) {
                //printf("b[%d] = %d\tp[%d] = %4.2e\n",i,b[i],i,p[i]);
                if (i < min_slice) {
                        current->service[_sort[i]] = 0;
                        s0+=min_b[i];
                }
                else if (i >= min_slice) {
                        current->service[_sort[i]] = 1;
                        s1+=min_b[i];
                }
                current->b[i] = min_b[_desort[i]];
                current->psd[i] = watts_to_dbmhz(min_p[_desort[i]]);
        }
        if (s0 != rate1) {
                printf("couldnt reach target data rate for service 0 rate = %d target_rate = %d\n",s0,rate1);
                //exit(1);
        }
        if (s1 != rate2) {
                printf("couldnt reach target data rate for service 1 rate = %d target_rate = %d\n",s1,rate2);
                //exit(1);
        }
	printf("s0 rate = %d\ts1 rate = %d\n",s0,s1);
}
int dbw_cbr(int line_id,double gamma1,double gamma2, int rate1, int rate2)
{
	struct line *current;
	int i,j,iters;
	double cnrs[DMTCHANNELS];
	int _sort[DMTCHANNELS];
	int _desort[DMTCHANNELS];
	double p[DMTCHANNELS];		// local array to power energy results from minimisation routine
	//double p1,p2;
	double p_low,p_high,p_min;
	double p_total;
	int b[DMTCHANNELS];
	int low = rate1/MAXBITSPERTONE+1;
	int high = DMTCHANNELS-1-(rate2/MAXBITSPERTONE);
	int step = (high - low)/4;
	int slice = DMTCHANNELS/2;
	int dir,prev_dir;
	int s0=0;
	int s1=0;
	current = get_line(line_id);
	current->is_dual = 1;		// remember that this line is dual qos!
	current->gamma[0] = gamma1;	// and remember gamma values! - higher gamma i.e lower BER
	current->gamma[1] = gamma2;	// lower gamma i.e. higher BER
	memset(p,0,(sizeof(double) * DMTCHANNELS));
	sort(current->cnr,_sort,_desort);	// create sort and desort arrays for cnrs for line
	for (i=0;i<DMTCHANNELS;i++) {
                cnrs[i] = current->cnr[_sort[i]];       // copy sorted cnrs to local array
        }
	//p_high = dbw_load(p,cnrs,b,slice,rate1,rate2,gamma1,gamma2);
	iters = 0;
	while (1) {
		iters++;
		high = slice + step;
		low = slice - step;
		p_high = dbw_load(p,cnrs,b,high,rate1,rate2,gamma1,gamma2);
		p_low = dbw_load(p,cnrs,b,low,rate1,rate2,gamma1,gamma2);
		printf("high = %d p_high =%lf\tlow = %d p_low = %lf\tstep = %d\n",high,p_high,low,p_low,step);
		if ((p_high > p_min) && (p_low > p_min) && (step == 1)) {
			//printf("done\n");
			break;
		}
		//prev_dir = dir;
		if (p_high < p_low) {
			slice = high;
			p_min = p_high;
			dir = HIGH;
		}
		else {
			slice = low;
			p_min = p_low;
			dir = LOW;
		}
		//if (dir == LOW)
		//	printf("Direction = down\n");
		//if (dir == HIGH)
		//	printf("Direction = up\n");
		if (step != 1) 		
			step = (high - low)/4;
	}
	p_min = dbw_load(p,cnrs,b,slice,rate1,rate2,gamma1,gamma2);	//although we know the slice, have to recalculate p and b for correct slice point
	printf("p_min = %lf mW\tslice point = %d\n\n",p_min * 1e3,slice);
/*
	for (i=0;i<DMTCHANNELS;i++) {
		double m,g;
		int s;
		if (i < slice) {
			g = gamma1;
			m = 10*log10(p[i]*cnrs[i]/(pow(2,b[i]) - 1)) - gamma1 + C_G;
			s = 0;
		}
		else if (i >=slice) {
			g = gamma2;
			m = 10*log10(p[i]*cnrs[i]/(pow(2,b[i]) - 1)) - gamma2 + C_G;
			s = 1;
		}
		printf("p[%d] = %lf\tb[%d] = %d\tcnr[%d] = %4.2e\tg = %4.2lf\tm = %4.2lf\ttone = %d\tservice = %d\n",i,watts_to_dbmhz(p[i]),i,b[i],i,cnrs[i],g,m,_sort[i],s);
	}
*/
	for (i=0;i<DMTCHANNELS;i++) {
		//printf("b[%d] = %d\tp[%d] = %4.2e\n",i,b[i],i,p[i]);
		if (i < slice) {
			current->service[_sort[i]] = 0;
			s0+=b[i];
		}
		else if (i >= slice) {
			current->service[_sort[i]] = 1;
			s1+=b[i];
		}
		current->b[i] = b[_desort[i]];
		current->psd[i] = watts_to_dbmhz(p[_desort[i]]);
	}
	if (s0 != rate1) {
		printf("couldnt reach target data rate for service 0 rate = %d target_rate = %d\n",s0,rate1);
		//exit(1);
	}
	if (s1 != rate2) {
		printf("couldnt reach target data rate for service 1 rate = %d target_rate = %d\n",s1,rate2);
		//exit(1);
	}
	/*for (i=0;i<DMTCHANNELS;i++) {
		printf("service[%d] = %d\n",i,current->service[i]);
	}*/
	printf("s0 rate = %d\ts1 rate = %d\n",s0,s1);
	return iters;	// number of iterations
}
static double dbw_load(double *p,double *cnrs,int *b,int slice,int rate1,int rate2,double gamma1,double gamma2)
{
	double p1;
	double p2;
	p1 = am_p_min_fm(p,cnrs,b,slice,rate1,gamma1,TOP_TONES);              //top half of channels, i.e. highest SNRs
        p2 = am_p_min_fm(&p[slice],&cnrs[slice],&b[slice],(DMTCHANNELS-slice),rate2,gamma2,BOTTOM_TONES);	//bottom half, lower SNRs
	//printf("Slice = %d\tp = %lf\n\n",slice,p1+p2);
	return p1+p2;
}
static double am_p_min_fm(double *p,double *cnr,int *b,int n,int rate,double gamma,int half)
{
        int i,tone,b_total=0,slice;
        int tone_full[n];             // set to 1 if tone is at MAXBITSPERTONE
        //double gamma = 9.8;
        double margin = MARGIN;
        double c_g = C_G;
        double gamma_hat = pow(10,(gamma+margin-c_g)/10);
        double delta_p[n];
        double p_tot = 0.0,p_budget=0.1;        // 0.1 W = 100mW
	char str[15];
	for (i=0;i<n;i++) {
		b[i] = 0;		// zero all bits before loading
		tone_full[i] = 0;	// all tones are not full
		p[i] = 0;		// zero all powers before loading
	}
	if (half == TOP_TONES)
		strcpy(str,"TOP TONES");
	else if (half == BOTTOM_TONES)
		strcpy(str,"BOTTOM TONES");
	for (i=0;i<n;i++) {		// just calculate full delta_p array once then update one member at a time inside loop
        	if (!tone_full[i]) {
                	delta_p[i] = (pow(2,(b[i] + 1))-1) * gamma_hat/cnr[i] - (pow(2,b[i])-1) * gamma_hat/cnr[i];
                }
                else {
                        delta_p[i] = 100;
                }
        }
        while (1) {
                if ((tone = min_index(delta_p,n)) == -1) {
                        printf("All tones are full!\n");
			printf("rate = %d\n",b_total);
                        break;
                }
		/*
		if (tone == 100) {
                        double p,p1;
                        printf("Added a bit on tone %d\n",tone);
			printf("cnr = %lf\tgamma = %lf\thalf = %s\n",cnr[tone],gamma_hat,str);
                        p = (pow(2,b[tone])-1) * gamma_hat/cnr[tone];
                        p1 = (pow(2,(b[tone]+1))-1) * gamma_hat/cnr[tone];
                        printf("b = %d\t\tp = %e = %lf\n",b[tone],p,watts_to_dbmhz(p));
                        printf("new_b = %d\tnew_p = %e = %lf\n",b[tone]+1,p1,watts_to_dbmhz(p1));
			printf("margin = %4.2lf\n\n",10*log10(p1*cnr[tone]/(pow(2,b[tone]+1) - 1)) - gamma + C_G);
                }*/
		/*
		{
			double p,p1;
			p1 = (pow(2,(b[tone]+1))-1) * gamma_hat/cnr[tone];
			if ((10*log10(p1*cnr[tone]/(pow(2,b[tone]+1) - 1))) - gamma + C_G < 2.9) {
				printf("Added a bit on tone %d\n",tone);
                        	printf("cnr = %lf\tgamma = %lf\thalf = %s\n",cnr[tone],gamma_hat,str);
				p = (pow(2,b[tone])-1) * gamma_hat/cnr[tone];
				p1 = (pow(2,(b[tone]+1))-1) * gamma_hat/cnr[tone];
				printf("b = %d\t\tp = %e = %lf\n",b[tone],p,watts_to_dbmhz(p));
				printf("new_b = %d\tnew_p = %e = %lf\n",b[tone]+1,p1,watts_to_dbmhz(p1));
				printf("margin = %4.2lf\n\n",10*log10(p1*cnr[tone]/(pow(2,b[tone]+1) - 1)) - gamma + C_G);
			}
		}*/
                b[tone]++;
		b_total++;
		p_tot += delta_p[tone];
                p[tone] += delta_p[tone];
		if (b_total == rate) {
			break;
		}
                //printf("Added a bit on tone %d\n",tone);
                if (b[tone] == MAXBITSPERTONE)
                        tone_full[tone] = 1;
                if (p_tot > p_budget) {
			b[tone]--;
			b_total--;
			p_tot -= delta_p[tone];
			p[tone] -= delta_p[tone];
                        break;
		}
		if (!tone_full[tone]) {
			delta_p[tone] = (pow(2,(b[tone] + 1))-1) * gamma_hat/cnr[tone] - (pow(2,b[tone])-1) * gamma_hat/cnr[tone];
                }
                else {
                        delta_p[tone] = 100;
                }
        }
	if (half == TOP_TONES)
		slice = n;
	else if (half == BOTTOM_TONES)
		slice = DMTCHANNELS - n;
	if (b_total != rate) {
		//printf("could not reach target data rate for p_budget = %lfb_total = %d\trate = %d\tslice = %d\n",p_budget,b_total,rate,slice);
	}
	return p_tot;
}
static void sort(double g_n[], int *chanorder, int *desort)
{
        int temp1,temp2,i,j;
        bool swapped = true;
        double g_n1[DMTCHANNELS];
        memcpy(g_n1,g_n,sizeof(g_n1));
        for (i=0;i<DMTCHANNELS;i++) {
                *(chanorder + i) = i;
        }
        while (swapped == true) {
                swapped = false;
                for (i=0;i<DMTCHANNELS-1;i++) {
                        if (g_n1[i] < g_n1[i+1]) {
                                temp1 = g_n1[i];
                                g_n1[i] = g_n1[i+1];
                                g_n1[i+1] = temp1;
                                temp2 = *(chanorder+i);
                                *(chanorder+i) = *(chanorder+i+1);
                                *(chanorder+i+1) = temp2;
                                swapped = true;
                        }
                }
        }
        /*for (i=0;i<DMTCHANNELS;i++) {
                printf("b_hat[%d] = %d\t",i,b_hat[i]);
        }*/
        for (i=0,j=0;i<DMTCHANNELS;) {
                if (chanorder[j] == i) {
                        desort[i] = j;
                        i++;
                        j=0;
                }
                else {
                        j++;
                        continue;
                }
        }
        /*      
        for (i=0;i<DMTCHANNELS;i++) {
                printf("chanorder[%d] = %d\tchanorder[desort[%d]] = %d\n",i,*(chanorder+i),i,chanorder[desort[i]]);
        }*/
}
static int min_index(double *delta_p,int n)
{
        int i,min_index=-1;
        double min = 1e100;
        for (i=0;i<n;i++) {
                if (delta_p[i] < min) {
                        min = delta_p[i];
                        min_index=i;
                }
        }
        return min_index;
}
