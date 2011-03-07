#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "multiuser_load.h"

static void sort(double *,int *,int *);
static int min_index(double *,int);
static int am_p_min_fm(double *p,double *cnr,int *b,int n,int rate,double gamma,int half,double *p_used);
static int dbw_vbr_load(double *p,double *cnrs,int *b,int slice,int rate2,double gamma1,double gamma2);

#define TOP_TONES 0	// higher snr service
#define BOTTOM_TONES 1	// lower snr service
#define HIGH 1
#define LOW 0

static double tot_p_budget;

void dbw_vbr_cbr_exhaustive(int line_id,double gamma1,double gamma2, int rate2,double p_budget)
{

	struct line *current;
	int i,j,slice,max_slice,s0,s1;
	double cnrs[DMTCHANNELS];
	int _sort[DMTCHANNELS];
	int _desort[DMTCHANNELS];
	double p[DMTCHANNELS];		// local array to power energy results from minimisation routine
	double min_p[DMTCHANNELS];
	double p0,p1;
	double p_tot;
	double p_min_tot = 1e100;	// this is the minimum power total for all channels, set to large value to start
	int b[DMTCHANNELS];
	int min_b[DMTCHANNELS];
	int rate_top=0,rate_bot=0,max_rate=0;
	int rate_bot_got=0;
	double p_used;

	tot_p_budget=p_budget;		// uh oh terrible hack

	current = get_line(line_id);
	
	current->is_dual = 1;           // remember that this line is dual qos!
        current->gamma[0] = gamma1;     // and remember gamma values! - higher gamma i.e lower BER
        current->gamma[1] = gamma2;     // lower gamma i.e. higher BER
	strcpy(current->loading_algo,"dbw_vbr_cbr_exhaustive");
	

	memset(p,0,(sizeof(double) * DMTCHANNELS));

	sort(current->cnr,_sort,_desort);	// create sort and desort arrays for cnrs for line

	for (i=0;i<DMTCHANNELS;i++) {
		cnrs[i] = current->cnr[_sort[i]];	// copy sorted cnrs to local array
	}


	for (slice=1;slice<DMTCHANNELS;slice++) {

                if ((DMTCHANNELS-slice)*MAXBITSPERTONE < rate2) {
                        //printf("%d channels cannot support rate of %d\n",(DMTCHANNELS-i),rate2);
                        continue;
                }

 		rate_bot = am_p_min_fm(&p[slice],&cnrs[slice],&b[slice],(DMTCHANNELS-slice),rate2,gamma2,BOTTOM_TONES,&p_used); //bottom half, lower SNRs
	
		if (rate_bot != rate2) {
			//printf("slice %d could not support rate %d\n",slice,rate2);
			continue;
		}

		rate_bot_got=1;

        	rate_top = am_p_min_fm(p,cnrs,b,slice,0,gamma1,TOP_TONES,&p_used);              //top half of channels, i.e. highest SNRs VBR service

		//printf("slice = %d\trate_top = %d\n",slice,rate_top);

		//printf("slice = %d\tpower = %lf\n",i,p1+p2);

		for (j=0,p_tot=0.0;j<DMTCHANNELS;j++) {
                        p_tot += p[j];
                }
		//printf("p_tot = %lf\ti = %d\n",p_tot,i);

		if (rate_top > max_rate) {
			max_rate=rate_top;
                        max_slice = slice;
                        memcpy(min_p,p,(sizeof(double) * DMTCHANNELS));
                        memcpy(min_b,b,(sizeof(double) * DMTCHANNELS));
                        //printf("New p_min, i = %d\n",i);
                }

	
	}

	if (rate_bot_got == 0) {
		printf("Couldnt achieve rate 2 target for any slice point\n");
		memset(current->b,0,sizeof(int)*DMTCHANNELS);
		memset(current->psd,0,sizeof(double)*DMTCHANNELS);
		return;
	}

	//printf("max rate = %d\tslice point = %d\n",max_rate,max_slice);

	p0=0.0;
	p1=0.0;
	s0=0;
	s1=0;
	for (i=0;i<DMTCHANNELS;i++) {
                //printf("b[%d] = %d\tp[%d] = %4.2e\n",i,b[i],i,p[i]);
                if (i < max_slice) {
                        current->service[_sort[i]] = 0;
                        s0+=min_b[i];
                        p0+=min_p[i];
                }
                else if (i >= max_slice) {
                        current->service[_sort[i]] = 1;
                        s1+=min_b[i];
                        p1+=min_p[i];

                }
                current->b[i] = min_b[_desort[i]];
                current->psd[i] = watts_to_dbmhz(min_p[_desort[i]]);
        }


        printf("service 0 rate = %d service 1 rate = %d\n",s0,s1);
        //printf("service 0 power = %lf service 1 power = %lf total power = %lf\n",p0,p1,(p0+p1));



}

int dbw_vbr_cbr(int line_id,double gamma1,double gamma2, int rate2)
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
	int low=0;				
	int high=DMTCHANNELS-1-rate2/MAXBITSPERTONE;	// check this as slice_max will be truncated
	int step = (high-low)/4;
	int slice=(high-low)/2;
	int rate_high,rate_low,max_rate;
	int dir,prev_dir;
	int s0=0,s1=0;
	double p0=0.0,p1=0.0;

	current = get_line(line_id);

	current->is_dual = 1;		// remember that this line is dual qos!
	current->gamma[0] = gamma1;	// and remember gamma values! - higher gamma i.e lower BER
	current->gamma[1] = gamma2;	// lower gamma i.e. higher BER
	strcpy(current->loading_algo,"dbw_vbr_cbr");


	memset(p,0,(sizeof(double) * DMTCHANNELS));

	sort(current->cnr,_sort,_desort);	// create sort and desort arrays for cnrs for line

	for (i=0;i<DMTCHANNELS;i++) {
                cnrs[i] = current->cnr[_sort[i]];       // copy sorted cnrs to local array
        }


	iters = 0;
	while (1) {
		high = slice + step;
		low = slice - step;
		//printf("High = %d\n",high);
		//printf("Low = %d\n",low);

		rate_high = dbw_vbr_load(p,cnrs,b,high,rate2,gamma1,gamma2);
		rate_low = dbw_vbr_load(p,cnrs,b,low,rate2,gamma1,gamma2);
		
		//printf("rate_high = %d\n",rate_high);
		//printf("rate_low = %d\n",rate_low);


		if ((rate_high <= max_rate) && (rate_low <= max_rate) && (step == 1)) {
			break;
		}

		if (rate_high > rate_low) {
			slice = high;
			max_rate=rate_high;
		}
		else {
			slice = low;
			max_rate=rate_low;
		}
		if (step != 1)
			step = (high-low)/4;


	}

	dbw_vbr_load(p,cnrs,b,slice,rate2,gamma1,gamma2);	// we know the slice point but need to recalculate p[] and b[] for this point

	printf("max_rate = %d\tslice point = %d\n",max_rate,slice);



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
		if (i==slice-1)
			printf("slice point\n");
	}
*/

	for (i=0;i<DMTCHANNELS;i++) {
		//printf("b[%d] = %d\tp[%d] = %4.2e\n",i,b[i],i,p[i]);
		if (i < slice) {
			current->service[_sort[i]] = 0;
			s0+=b[i];
			p0+=p[i];
		}
		else if (i >= slice) {
			current->service[_sort[i]] = 1;
			s1+=b[i];
			p1+=p[i];

		}
		current->b[i] = b[_desort[i]];
		current->psd[i] = watts_to_dbmhz(p[_desort[i]]);
	}


	printf("service 0 rate = %d service 1 rate = %d\n",s0,s1);
	printf("service 0 power = %lf service 1 power = %lf total power = %lf\n",p0,p1,(p0+p1));

	/*for (i=0;i<DMTCHANNELS;i++) {
		printf("service[%d] = %d\n",i,current->service[i]);
	}*/
	return iters;	// number of iterations


}

static int dbw_vbr_load(double *p,double *cnrs,int *b,int slice,int rate2,double gamma1,double gamma2)
{

	int rate_top;
	int rate_bot;
	double p_used;		//power used by CBR service

        rate_bot = am_p_min_fm(&p[slice],&cnrs[slice],&b[slice],(DMTCHANNELS-slice),rate2,gamma2,BOTTOM_TONES,&p_used);	//bottom half, lower SNRs
	
	rate_top = am_p_min_fm(p,cnrs,b,slice,0,gamma1,TOP_TONES,&p_used);              //top half of channels, i.e. highest SNRs VBR service


	//printf("Slice = %d\tp = %lf\n\n",slice,p1+p2);

	if (rate_bot == -1)
		return -1;		// current slice point could not support CBR service

	return rate_top;

}


static int am_p_min_fm(double *p,double *cnr,int *b,int n,int rate,double gamma,int half,double *p_used)
{

        int i,tone,b_total=0,slice;
        int tone_full[n];             // set to 1 if tone is at MAXBITSPERTONE
        //double gamma = 9.8;
        double margin = MARGIN;
        double c_g = C_G;
        double gamma_hat = pow(10,(gamma+margin-c_g)/10);
        double delta_p[n];
        double p_tot = 0.0,p_budget;        // 0.1 W = 100mW
	char str[15];

	for (i=0;i<n;i++) {
		b[i] = 0;		// zero all bits before loading
		tone_full[i] = 0;	// all tones are not full
		p[i] = 0;		// zero all powers before loading
	}
 	

	if (half == TOP_TONES) {
		strcpy(str,"TOP TONES");
		p_budget=tot_p_budget-*p_used;
	}
	else if (half == BOTTOM_TONES) {
		strcpy(str,"BOTTOM TONES");
		p_budget=tot_p_budget;
	}
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

                b[tone]++;
		b_total++;
		p_tot += delta_p[tone];
                p[tone] += delta_p[tone];

		
		if ((half == BOTTOM_TONES) && (b_total == rate)) {
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

	if (half == TOP_TONES) {
		slice = n;
		if (b_total != rate) {
			//printf("could not reach target data rate for p_budget = %lf\nb_total = %d\trate = %d\tslice = %d\n",p_budget,b_total,rate,slice);
			b_total=-1;	// return negative if could not reach target rate
		}
	}
	else if (half == BOTTOM_TONES) {
		slice = DMTCHANNELS - n;
		*p_used=p_tot;
	
	}


	return b_total;

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


