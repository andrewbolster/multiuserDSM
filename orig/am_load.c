#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "multiuser_load.h"
int am_load_ra(int line_id, double p_budget,int half)
{
	int i,tone,b_total=0;
	int tone_full[DMTCHANNELS];		// set to 1 if tone is at MAXBITSPERTONE
	double gamma = GAMMA;
	double margin = MARGIN;
	double c_g = C_G;
	double gamma_hat = pow(10,(gamma+margin-c_g)/10);
	double delta_p[DMTCHANNELS];
	double p[DMTCHANNELS];
	double p_tot = 0.0;
	double min_psd = -60.0;
	struct line* current;
	double mask = dbmhz_to_watts(-20);
	memset(p,0,DMTCHANNELS*sizeof(double));
	current = get_line(line_id);
	current->is_dual = 0;	// just single qos
	current->gamma[0] = gamma;
	strcpy(current->loading_algo,"rate adaptive");
	for (i=0;i<DMTCHANNELS;i++) {
		current->b[i] = 0;	// set all bits to zero
		tone_full[i] = 0;	// all tones are not full
	}
/*	if ((half == TOP_HALF) && DMTCHANNELS == 512) {
		for (i=0;i<DMTCHANNELS/2;i++) {
			tone_full[i] = 1;
		}
	}
	else if ((half == BOTTOM_HALF) && DMTCHANNELS == 512) {
		for (i=DMTCHANNELS/2;i<DMTCHANNELS;i++) {
			tone_full[i] = 1;
		}
	}*/
	for (i=0;i<DMTCHANNELS;i++) {		// calc delta_p array once then update each member individually
        	if (!tone_full[i]) {
			delta_p[i] = (pow(2,(current->b[i] + 1))-1) * gamma_hat/current->cnr[i] - (pow(2,current->b[i])-1) * gamma_hat/current->cnr[i];
		}
		else {
			delta_p[i] = 100;
		}
	}
	while (1) {
		if ((tone = find_min(delta_p)) == -1) {
			printf("All tones are full!\n");
			break;
		}
		current->b[tone]++;
		p_tot += delta_p[tone];
		p[tone] += delta_p[tone];
		/*if (tone == 0) {
			printf("Added a bit on tone %d line %d\n",tone,current->line_id);
			printf("cnr = %4.2e\n",current->cnr[tone]);
			printf("b = %d\n",current->b[tone]);
			printf("psd = %4.2lf\n",watts_to_dbmhz((pow(2,current->b[tone])-1) * (gamma_hat/current->cnr[tone])));
			printf("margin = %4.3lf\n",10*log10(current->cnr[tone]*(pow(2,current->b[tone])-1) * gamma_hat/current->cnr[tone]/(pow(2,current->b[tone]) - 1)) - current->gamma[0]);
		}*/
		if (current->b[tone] == MAXBITSPERTONE)
			tone_full[tone] = 1;
		if (p_tot > p_budget) {
			current->b[tone]--;
			p_tot -= delta_p[tone];
			break;
		}
		else if (p[tone] > mask) {
			current->b[tone]--;
			p[tone]-=delta_p[tone];
			p_tot -= delta_p[tone];
			tone_full[tone]=1;
		}
		if (!tone_full[tone]) {
                        delta_p[tone] = (pow(2,(current->b[tone] + 1))-1) * gamma_hat/current->cnr[tone] - (pow(2,current->b[tone])-1) * gamma_hat/current->cnr[tone];
                }
                else {
                        delta_p[tone] = 100;
                }
	}
	for (i=0;i<DMTCHANNELS;i++) {
		b_total += current->b[i];
		current->psd[i] = watts_to_dbmhz((pow(2,current->b[i])-1) * (gamma_hat/current->cnr[i]));
		if ((current->psd[i] < min_psd) && current->b[i] != 0)
			current->psd[i] = min_psd;
	}
	//printf("line_id = %d\trate = %d\n",line_id,b_total);
	for (i=0;i<DMTCHANNELS;i++) {
                current->service[i] = 0;
        }
	return b_total;	
}
int am_load_fm(int line_id,int rate,double p,int half)
{
        int i,tone,b_total=0;
        int tone_full[DMTCHANNELS];             // set to 1 if tone is at MAXBITSPERTONE
        double gamma = GAMMA;
        double margin = MARGIN;
        double c_g = C_G;
        double gamma_hat = pow(10,(gamma+margin-c_g)/10);
        double delta_p[DMTCHANNELS];
        double p_tot = 0.0,p_budget=p;        // 0.1 W = 100mW/FRAMERATE
        //double min_psd = -60.0;
        struct line* current;
        current = get_line(line_id);	
	current->is_dual = 0;		// just single qos level
	current->gamma[0] = gamma;
	strcpy(current->loading_algo,"fixed margin");
	if (rate == 0) {
		for (i=0;i<DMTCHANNELS;i++) {
			current->b[i]=0;
			current->psd[i]=0;
		}
		return 0;
	}
        for (i=0;i<DMTCHANNELS;i++) {
                current->b[i] = 0;      // set all bits to zero
                tone_full[i] = 0;       // all tones are not full
        }
	if ((half == TOP_HALF) && DMTCHANNELS == 512) {
                for (i=0;i<DMTCHANNELS/2;i++) {
                        tone_full[i] = 1;
                }
        }
        else if ((half == BOTTOM_HALF) && DMTCHANNELS == 512) {
                for (i=DMTCHANNELS/2;i<DMTCHANNELS;i++) {
                        tone_full[i] = 1;
                }
        }
	for (i=0;i<DMTCHANNELS;i++) {
		if (!tone_full[i]) {
			delta_p[i] = (pow(2,(current->b[i] + 1))-1) * gamma_hat/current->cnr[i] - (pow(2,current->b[i])-1) * gamma_hat/current->cnr[i];
		}
		else {
			delta_p[i] = 100;
		}
	}
        while (1) {
                if ((tone = find_min(delta_p)) == -1) {
                        printf("All tones are full!\n");
                        break;
                }
		/*if (tone == 224) {
			double p,p1;
			printf("cnr = %lf\tgamma = %lf\n",current->cnr[tone],gamma_hat);
			p = (pow(2,current->b[tone])-1) * gamma_hat/current->cnr[tone];
			p1 = (pow(2,(current->b[tone]+1))-1) * gamma_hat/current->cnr[tone];
			printf("Added a bit on tone %d\n",tone);
			printf("b = %d\t\tp = %e = %lf\n",current->b[tone],p,watts_to_dbmhz(p));
			printf("new_b = %d\tnew_p = %e = %lf\n",current->b[tone]+1,p1,watts_to_dbmhz(p1));
		}*/
                current->b[tone]++;
		b_total++;
                p_tot += delta_p[tone];
		if (b_total == rate)
			break;
                //printf("Added a bit on tone %d\n",tone);
                if (current->b[tone] == MAXBITSPERTONE)
                        tone_full[tone] = 1;
                if (p_tot > p_budget) {
                        current->b[tone]--;
			b_total--;
			p_tot -= delta_p[tone];
			break;
		}
		if (!tone_full[tone]) {
                        delta_p[tone] = (pow(2,(current->b[tone] + 1))-1) * gamma_hat/current->cnr[tone] - (pow(2,current->b[tone])-1) * gamma_hat/current->cnr[tone];
                }
                else {
                        delta_p[tone] = 100;
                }
        }
        for (i=0;i<DMTCHANNELS;i++) {
                current->psd[i] = watts_to_dbmhz((pow(2,current->b[i])-1) * (gamma_hat/current->cnr[i]));
		//printf("psd[%d] = %lf\n",i,current->psd[i]);
		//if ((current->psd[i] < min_psd) && current->b[i] != 0)
                        //current->psd[i] = min_psd;
        }
	for (i=0;i<DMTCHANNELS;i++) {
		current->service[i] = 0;
	}
        //printf("line_id = %d energy = %lf mW rate = %d\n",line_id,p_tot*1e3,b_total);
	if (b_total != rate)
		printf("could not reach target data rate\n");
	return b_total;
}
int find_min(double *delta_p)
{
        int i,min_index=-1;
        double min = 1e100;
        for (i=0;i<DMTCHANNELS;i++) {
                if (delta_p[i] < min) {
                        min = delta_p[i];
                        min_index=i;
                }
        }
        return min_index;
}
int no_fext_ra(int line_id, double p_budget,double *waterfill_level, int* active_channels)  // This should be a method of line but its not cos i cant be arsed
{
	struct line* current;
	current = get_line(line_id);
	int i,tone,b_total=0;
	int tone_full[DMTCHANNELS];		// set to 1 if tone is at MAXBITSPERTONE
	double gamma_hat = UNdB((current->gamma[0]+MARGIN));
	double delta_p[DMTCHANNELS];
	int b[DMTCHANNELS];
	double p_tot = 0.0;
	if (current->is_dual) {
		printf("Tried to call no_fext_ra on a dual qos line\n");
		exit(2);
	}
	check_line_sanity(current);	
	for (i=0;i<DMTCHANNELS;i++) {
		b[i]=0;
		tone_full[i] = 0;	// all tones are not full
	}
	for (i=0;i<DMTCHANNELS;i++) {		// calc delta_p array once then update each member individually
        	if (!tone_full[i]) {
			delta_p[i] = (pow(2,(b[i] + 1))-1) * gamma_hat/(current->gain[i]/(current->alien_xtalk[i]+current->bkn)) - (pow(2,b[i])-1) * gamma_hat/(current->gain[i]/(current->alien_xtalk[i]+current->bkn));
		}
		else {
			delta_p[i] = 100;
		}
	}
	while (1) {
		if ((tone = find_min(delta_p)) == -1) {
			printf("All tones are full!\n");
			break;
		}
		b[tone]++;
		p_tot += delta_p[tone];
/*		
		printf("Added a bit on tone %d line %d\n",tone,line_id);
		//printf("cnr = %4.2e\n",current->cnr[tone]);
		printf("b = %d\n",b[tone]);
		printf("psd = %4.2lf\n",watts_to_dbmhz(gamma_hat*(pow(2,b[tone])-1)/(current->gain[tone]/(current->alien_xtalk[tone]+current->bkn))));
		//printf("margin = %4.3lf\n",10*log10(current->cnr[tone]*(pow(2,current->b[tone])-1) * gamma_hat/current->cnr[tone]/(pow(2,current->b[tone]) - 1)) - current->gamma[0]);
		getchar();
*/	
		if (b[tone] == MAXBITSPERTONE)
			tone_full[tone] = 1;
		if (p_tot > p_budget) {
			b[tone]--;
			p_tot -= delta_p[tone];
			break;
		}
		if (!tone_full[tone]) {
                        delta_p[tone] = (pow(2,(b[tone] + 1))-1) * gamma_hat/(current->gain[tone]/(current->alien_xtalk[tone]+current->bkn)) - (pow(2,b[tone])-1) * gamma_hat/(current->gain[tone]/(current->alien_xtalk[tone]+current->bkn));
                }
                else {
                        delta_p[tone] = 100;
                }
	}
	*active_channels=0;
	for (i=0;i<DMTCHANNELS;i++) {
		b_total += b[i];
		if (b[i] > 0)
			(*active_channels)++;
	}
	*waterfill_level=p_tot/(*active_channels);
	if (b_total == MAXBITSPERTONE * DMTCHANNELS) {
		printf("Line %d can load the max rate, you might have convergence issues\n",line_id);
	}
	char tag[100];
	sprintf(tag,"no_fext_b_user_%d",line_id);
	dump_int_vector(b,tag);
	return b_total;	
}
double no_fext_ra_frac(int line_id, double p_budget)  // This should be a method of line but its not cos i cant be arsed
{
	struct line* current;
	current = get_line(line_id);
	int i,tone;
	double b_total=0;
	double bit_inc=0.01;
	int tone_full[DMTCHANNELS];		// set to 1 if tone is at MAXBITSPERTONE
	double gamma_hat = UNdB((current->gamma[0]+MARGIN));
	double delta_p[DMTCHANNELS];
	double b[DMTCHANNELS];
	double p_tot = 0.0;
	if (current->is_dual) {
		printf("Tried to call no_fext_ra on a dual qos line\n");
		exit(2);
	}
	check_line_sanity(current);	
	for (i=0;i<DMTCHANNELS;i++) {
		b[i]=0;
		tone_full[i] = 0;	// all tones are not full
	}
	for (i=0;i<DMTCHANNELS;i++) {		// calc delta_p array once then update each member individually
        	if (!tone_full[i]) {
			delta_p[i] = (pow(2,(b[i] + bit_inc))-1) * gamma_hat/(current->gain[i]/(current->alien_xtalk[i]+current->bkn)) - (pow(2,b[i])-1) * gamma_hat/(current->gain[i]/(current->alien_xtalk[i]+current->bkn));
		}
		else {
			delta_p[i] = 100;
		}
	}
	while (1) {
		if ((tone = find_min(delta_p)) == -1) {
			printf("All tones are full!\n");
			break;
		}
		b[tone]+=bit_inc;
		p_tot += delta_p[tone];
/*		
		printf("Added a bit on tone %d line %d\n",tone,line_id);
		//printf("cnr = %4.2e\n",current->cnr[tone]);
		printf("b = %d\n",b[tone]);
		printf("psd = %4.2lf\n",watts_to_dbmhz(gamma_hat*(pow(2,b[tone])-1)/(current->gain[tone]/(current->alien_xtalk[tone]+current->bkn))));
		//printf("margin = %4.3lf\n",10*log10(current->cnr[tone]*(pow(2,current->b[tone])-1) * gamma_hat/current->cnr[tone]/(pow(2,current->b[tone]) - 1)) - current->gamma[0]);
		getchar();
*/	
		if (b[tone] >= MAXBITSPERTONE)
			tone_full[tone] = 1;
		if (p_tot > p_budget) {
			b[tone]-=bit_inc;
			p_tot -= delta_p[tone];
			break;
		}
		if (!tone_full[tone]) {
                        delta_p[tone] = (pow(2,(b[tone] + bit_inc))-1) * gamma_hat/(current->gain[tone]/(current->alien_xtalk[tone]+current->bkn)) - (pow(2,b[tone])-1) * gamma_hat/(current->gain[tone]/(current->alien_xtalk[tone]+current->bkn));
                }
                else {
                        delta_p[tone] = 100;
                }
	}
	for (i=0;i<DMTCHANNELS;i++) {
		b_total += b[i];
	}
	if (b_total >= MAXBITSPERTONE * DMTCHANNELS) {
		printf("Line %d can load the max rate, you might have convergence issues\n",line_id);
	}
	return b_total;	
}
int waterfill(int user,double p_budget,int *_b,double *_p,double *noise)
{
	struct line* current;
	current = get_line(user);
	int i,tone,b_total=0;
	int tone_full[DMTCHANNELS];		// set to 1 if tone is at MAXBITSPERTONE
	double gamma_hat = UNdB((current->gamma[0]+MARGIN));
	double delta_p[DMTCHANNELS];
	double p[DMTCHANNELS];
	double n[DMTCHANNELS];
	int b[DMTCHANNELS];
	double p_tot = 0.0;
	memset(p,0,sizeof(double)*DMTCHANNELS);
	memset(b,0,sizeof(int)*DMTCHANNELS);
	if (noise==NULL) {
		memset(n,0,sizeof(double)*DMTCHANNELS);
	}
	else {
		memcpy(n,noise,sizeof(double)*DMTCHANNELS);
	}
	if (current->is_dual) {
		printf("Tried to call no_fext_ra on a dual qos line\n");
		exit(2);
	}
	check_line_sanity(current);	
	for (i=0;i<DMTCHANNELS;i++) {
		b[i]=0;
		tone_full[i] = 0;	// all tones are not full
	}
	for (i=0;i<DMTCHANNELS;i++) {		// calc delta_p array once then update each member individually
        	if (!tone_full[i]) {
			delta_p[i] = (pow(2,(b[i] + 1))-1) * gamma_hat/(current->gain[i]/(current->alien_xtalk[i]+current->bkn+n[i])) - (pow(2,b[i])-1) * gamma_hat/(current->gain[i]/(current->alien_xtalk[i]+current->bkn+n[i]));
		}
		else {
			delta_p[i] = 100;
		}
	}
	while (1) {
		if ((tone = find_min(delta_p)) == -1) {
			printf("All tones are full!\n");
			break;
		}
		p[tone]+=delta_p[tone];
		b[tone]++;
		p_tot += delta_p[tone];
/*		
		printf("Added a bit on tone %d line %d\n",tone,line_id);
		//printf("cnr = %4.2e\n",current->cnr[tone]);
		printf("b = %d\n",b[tone]);
		printf("psd = %4.2lf\n",watts_to_dbmhz(gamma_hat*(pow(2,b[tone])-1)/(current->gain[tone]/(current->alien_xtalk[tone]+current->bkn))));
		//printf("margin = %4.3lf\n",10*log10(current->cnr[tone]*(pow(2,current->b[tone])-1) * gamma_hat/current->cnr[tone]/(pow(2,current->b[tone]) - 1)) - current->gamma[0]);
		getchar();
*/	
		if (b[tone] == MAXBITSPERTONE)
			tone_full[tone] = 1;
		if (p_tot > p_budget) {
			b[tone]--;
			p_tot -= delta_p[tone];
			p[tone] -= delta_p[tone];
			break;
		}
		if (!tone_full[tone]) {
                        delta_p[tone] = (pow(2,(b[tone] + 1))-1) * gamma_hat/(current->gain[tone]/(current->alien_xtalk[tone]+current->bkn+n[i])) - (pow(2,b[tone])-1) * gamma_hat/(current->gain[tone]/(current->alien_xtalk[tone]+current->bkn+n[i]));
                }
                else {
                        delta_p[tone] = 100;
                }
	}
	for (i=0;i<DMTCHANNELS;i++) {
		b_total += b[i];
	}
	if (b_total == MAXBITSPERTONE * DMTCHANNELS) {
		printf("Line %d can load the max rate, you might have convergence issues\n",user);
	}
	char tag[100];
	sprintf(tag,"no_fext_b_user_%d",user);
	dump_int_vector(b,tag);
	if (_b != NULL) {
		memcpy(_b,b,sizeof(int)*DMTCHANNELS);
	}
	if (_p != NULL) {
		memcpy(_p,p,sizeof(double)*DMTCHANNELS);
	}
	return b_total;	
}
