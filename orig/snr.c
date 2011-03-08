#include "multiuser_load.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
int check_all_margins(double gap)
{
	struct line* current = list_head;
	int i,k;
	//int service;
	int ret = 0;
	double m;
	for (i=0;i<lines;i++) {
		current = get_line(i);
		if (current->is_dual == 0) {		
			for (k=0;k<DMTCHANNELS;k++) {
				if ((m = 10*log10(current->snr[k]/(pow(2,current->b[k]) - 1)) - current->gamma[0] + C_G) < MARGIN - gap) {
					/*printf("Tone %d Line %d is out of margin = %.15lf\n",k,i,m);
					for (int j=0;j<lines;j++) {
						current=get_line(j);
						printf("%d\t%4.2lf\t%4.2lf\n",current->b[k],current->psd[k],current->gamma_m[k]);
					}*/
					ret = 1;
					return ret;
				}
			}
		}
		else if (current->is_dual == 1) {
			for (k=0;k<DMTCHANNELS;k++) {
				if ((m = 10*log10(current->snr[k]/(pow(2,current->b[k]) - 1)) - current->gamma[current->service[k]] + C_G) < MARGIN - gap) {
					printf("Tone %d Line %d is out of margin = %.15lf\n",k,i,m);
					printf("gamma = %4.2lf\tservice = %d\n",current->gamma[current->service[k]],current->service[k]);
					ret = 1;
					return ret;
				}
			}
		}	
	}
	return ret;
}
void calculate_snr()
{
	int i,k,j;
	struct line* current = list_head;
	double snr[DMTCHANNELS];
	double noise;
	//FILE *fp;
	//char fn[255];
	// snr on tone k and line i
	for (i=0;i<lines;i++) {
		current = get_line(i);
		/*sprintf(fn,"line%d_snrs.txt",i);
		if ((fp = fopen(fn,"w")) == NULL) {
			printf("Cant open file for writing\n");
			exit(1);
		}		
		fprintf(fp,"#tone\tsnr\tpsd\n");*/
		check_line_sanity(current);
		current->p_total = 0.0;
		current->b_total = 0;
		for (j=0;j<8;j++) {
			current->rate[j] = 0;
			current->_rate[j] = 0;
		}
		for (k=0;k<DMTCHANNELS;k++) {
			//snr[k] = (dbmhz_to_watts(current->psd[k]) * current->gain[k]) / (calc_fext_noise(i,k) + 1e-14)
			//noise = fsan_sum(calc_fext_noise(i,k,lines,channel_matrix),alien_xtalk(i,k)) + dbmhz_to_watts(current->background_noise);	// -140dbM/Hz noise floor
			noise = calc_fext_noise(i,k,lines,channel_matrix,current) + alien_xtalk(i,k)  + dbmhz_to_watts(current->background_noise);	// -140dbM/Hz noise floor
			snr[k] = dbmhz_to_watts(current->psd[k]) * current->gain[k]/noise;	// snr calculation	
			current->snr[k] = snr[k];						// store snr in struct line
			current->cnr[k] = current->gain[k]/noise;				// calculate cnr
			current->gamma_m[k] = 10*log10(current->snr[k]/(pow(2,current->b[k]) - 1)) - current->gamma[current->service[k]] + C_G;
			if (current->b[k] != 0)
				current->symerr[k] = symerr(current->snr[k],current->b[k]);
			else
				current->symerr[k] = 0;
			//fprintf(fp,"%d\t%lf dB\t%lf\n",k,10*log10(snr[k]),current->psd[k]);
			current->p_total += dbmhz_to_watts(current->psd[k]);
			current->b_total += current->b[k];
			current->rate[current->service[k]] += current->b[k];
			if (current->is_frac)	
				current->_rate[current->service[k]] += current->_b[k];
		}
		//fclose(fp);
	}
}
double calc_fext_noise(int line_id, int channel, int lines, double *channel_matrix,struct line* victim)
{
	struct line *current;
	int i;
	double xtalk_gain;
	double noise[lines];
	double n=0;
	for (i=0;i<lines;i++) {
		noise[i] = 0;
		if (i != line_id) {
			//xtalk_gain = *(channel_matrix + channel + (i * DMTCHANNELS) + (lines * line_id * DMTCHANNELS));
			xtalk_gain = get_xtalk_gain(i,line_id,channel);
			current = get_line(i);						// i is line_id of xtalking line
			#ifdef FSAN
			noise[i] = dbmhz_to_watts(current->psd[channel]) * xtalk_gain;		// for use with fsan_sum of all xtalkers
			#endif
			#ifndef FSAN
			n += dbmhz_to_watts(current->psd[channel]) * xtalk_gain; 		// direct summation
			//victim->xtalk_noise[i][channel] = dbmhz_to_watts(current->psd[channel])*xtalk_gain;
			#endif
			//printf("xtalk gain between from line %d into %d on channel %d = %e\n",i,line_id,channel,xtalk_gain);
			//printf("psd on line %d on channel %d = %lf dBm/Hz = %e watts\n",i,channel,current->psd[channel],dbmhz_to_watts(current->psd[channel]));
			//printf("fext noise from line %d on channel %d = %e watts\n",i,channel,dbmhz_to_watts(current->psd[channel]) * xtalk_gain);
		}
	}
	#ifdef FSAN
	n = fsan_sum_i(noise);		// change this for direct summation
	#endif
	return n;	// watts
}	
/*
Assume co-located xtalker with alien noise profile given in vdsl functional requirements doc
Using model A, assume FEXT only from the NT side, ie XA.NT.A model (not really sure this is right, prob should include NEXT from alien sources)
*/
double alien_xtalk(int line_id, int channel)
{
	struct line *current;
	double xtalk_gain;	
	double noise_psd;
	double h2;
	//double freq = 2156.25 + 4312.5 * channel;
	double freq = 140156.25 + 4312.5 * channel;
	double noise = 0;
	current = get_line(line_id);
	h2 = insertion_loss(2,current->length/1000,freq);
	//xtalk_gain = fext(freq,current->length/1000,h2);	// fext gain from fictional co-located xtalker
	xtalk_gain = h2 * fext(freq,current->length/1000);	// fext gain from fictional co-located xtalker
	noise_psd = psd_alien_xtalk(freq,NT,current->alien_xtalk_model,FEXT);		
	//noise = xtalk_gain * dbmhz_to_watts(noise_psd);
	noise = UNdB(-70) * dbmhz_to_watts(noise_psd);
	return noise;	//watts
}
/* this is not quite right but not far off, should be linear with a log freq scale */
double psd_alien_xtalk(double freq, int end, int model, int direction)
{
        double psd=0.0;		// shut up compiler
	end=0;			// shut up compiler
	direction=0;
	switch(model) {
	case MODEL_A:
		if (freq <= 50e3) {
			psd = -22.2;
		}
		else if ((freq > 50e3) && (freq <= 75e3)) {
			psd = -22.2 + (freq-50e3)/(75e3-50e3)*(-30.6 - (-22.2));
		}
		else if ((freq > 75e3) && (freq <= 100e3)) {
			psd = -30.6 + (freq-75e3)/(100e3-75e3)*(-34.2 - (-30.6));
		}
		else if ((freq > 100e3) && (freq <= 292e3)) {
			psd = -34.2 + (freq-100e3)/(292e3-100e3)*(-35.3 - (-34.2));
		}
		else if ((freq > 292e3) && (freq <=400e3)) {
			psd = -35.3 + (freq-292e3)/(400e3-292e3)*(-43.7 - (-35.3));
		}
		else if ((freq > 400e3) && (freq <= 1104e3)) {
			psd = -43.7 + (freq-400e3)/(1104e3-400e3)*(-52.6 - (-43.7));
		}
		else if ((freq > 1104e3) && (freq <= 2500e3)) {
			psd = -52.6 + (freq-1104e3)/(2500e3-1104e3)*(-99.6 - (-52.6));
		}
		else if ((freq > 2500e3) && (freq <= 3637e3)) {
			psd = -99.6  + (freq-2500e3)/(3637e3-2500e3)*(-111.3 - (-99.6));
		}
		else if ((freq > 3637e3) && (freq <= 30000e3)) {
			psd = -111.3 + (freq-3637e3)/(30000e3-3637e3)*(-111.5 - (-111.3));
		}
		return psd;
	case NONE:
		return 0;
	default:
		printf("oops whats going on?\n");
		exit(2);
	}
}
void set_psd_all_lines(double psd,int lines)
{
	int i,j;
	struct line *current;
	for (i=0;i<lines;i++) {         // set all psds to -36.5dbm/Hz just for testing
                current = get_line(i);
                for(j=0;j<DMTCHANNELS;j++) {
                        current->psd[j] = psd;
                }
        }
}
void set_psd_on_line(int line_id,double psd)
{
	struct line* current;
	int i;
	if ((current = get_line(line_id)) == NULL) {
		printf("Line_id %d does not exist\n",line_id);
	}
	for (i=0;i<DMTCHANNELS;i++) {
		current->psd[i] = psd;
	}
}
double fsan_sum(double pow1, double pow2)
{
	double kn = 1/0.6;
	double fsan_sum;
	pow1 /= CHANNEL_BANDWIDTH;	// convert from watts to watts/Hz
	pow2 /= CHANNEL_BANDWIDTH;
	pow1 = pow(pow1,kn);
	pow2 = pow(pow2,kn);		// p^kn
	fsan_sum = pow((pow1 + pow2),1/kn) * CHANNEL_BANDWIDTH;	// back to watts
	return fsan_sum;
}
double fsan_sum_i(double *p)
{
        double kn = 1/0.6;
        double fsan_sum;
	double pow_p[lines];
	double pow_sum=0;
	int i;
	for (i=0;i<lines;i++) {
        	p[i] /= CHANNEL_BANDWIDTH;      // convert from watts to watts/Hz
	}
        for (i=0;i<lines;i++) {            // p^kn
		pow_p[i] = pow(p[i],kn);
	}
	for (i=0;i<lines;i++) {
		pow_sum += pow_p[i];
	}
        fsan_sum = pow((pow_sum),1/kn) * CHANNEL_BANDWIDTH; // back to watts
        return fsan_sum;
}
double dbmhz_to_watts(double psd)
{
	double watts;
	if (psd == 0)			// using 0dBm/Hz to represent off state as this psd level is unrealistic anyway
		return 0;
	watts = pow(10,psd/10) * 1e-3 * CHANNEL_BANDWIDTH;
	return watts;
}
double watts_to_dbmhz(double e)
{
	double psd;
	psd = 10*log10((e*1e3)/CHANNEL_BANDWIDTH);
	return psd;
}
void margin_tighten(int line_id)
{
	struct line *current;
	int k;
	double p_new=0.0,p=0.0;
	double gamma_hat;
	current = get_line(line_id);
	for (k=0;k<DMTCHANNELS;k++) {
		gamma_hat = pow(10,(current->gamma[current->service[k]]+MARGIN-C_G)/10);
		if ((current->gamma_m[k] - MARGIN) > 0.01) { // FIXME to use target margin from line struct,Anyway, this tones margin is more than 0.1dB greater than target
			p_new = (pow(2,current->b[k]) - 1) * (gamma_hat/current->cnr[k]);
//			printf("tone = %d p_old = %lf p_new = %lf\n",k,dbmhz_to_watts(current->psd[k]),p_new);
			p+=dbmhz_to_watts(current->psd[k]) - p_new;     // running total of collected power due to reducing margin
			current->psd[k] = watts_to_dbmhz(p_new);	// assign new psd value
		}
	}
//	printf("After one pass, we collected %lf watts due to excess margin\n",p);
	for (k=0;k<DMTCHANNELS;k++) {
		gamma_hat = pow(10,(current->gamma[current->service[k]]+MARGIN-C_G)/10);
		if ((MARGIN - current->gamma_m[k]) > 0.01) { // FIXME to use target margin from line struct,Anyway, this tones margin is more than 0.1dB less than target
			p_new = (pow(2,current->b[k]) - 1) * (gamma_hat/current->cnr[k]);
//			printf("tone = %d p_old = %lf p_new = %lf\n",k,dbmhz_to_watts(current->psd[k]),p_new);
			p -= p_new - dbmhz_to_watts(current->psd[k]);
			if (p < 0) {		// not enough power to bring this tone up to margin so delete a bit
				p += p_new - dbmhz_to_watts(current->psd[k]);	// add unused power back
				current->b[k]--;				// delete bit
//				printf("Deleted a bit on tone %d\n",k);
				p_new = (pow(2,current->b[k]) - 1) * (gamma_hat/current->cnr[k]); // calc new power
				p += dbmhz_to_watts(current->psd[k]) - p_new;			  // add gained power!
			}		
			current->psd[k] = watts_to_dbmhz(p_new);				// assign new psd value		
		}
	}
	//printf("p left over = %lf\n",p);
	//exit(0);
	//printf("margin tightened on line %d\n",line_id);
}
void check_line_sanity(struct line * current)
{
	for (int i=0;i<DMTCHANNELS;i++) {
		if (current->b[i] < 0) {
			printf("WTF??\n");
			getchar();
		}
		if (current->gain[i] <= 0) {
			printf("ORLY??\n");
			getchar();
		}
	}
	if (current->background_noise == NOT_SET || current->bkn == NOT_SET) {
		printf("Background noise level is not set on line %d!\n",current->line_id);
		exit(1);
	}
	if (current->alien_xtalk_model == NOT_SET) {
		printf("Alien xtalk model is not set on line %d!\n",current->line_id);
		exit(1);
	}
}
void set_background_noise(double level)
{
	int i;
	struct line *current;
	for (i=0;i<lines;i++) {
		current=get_line(i);
		current->background_noise=level;
		current->bkn=dbmhz_to_watts(level);
	}
}
void set_alien_xtalk_model(int model)
{
	int i;
	struct line *current;
	for (i=0;i<lines;i++) {
		current=get_line(i);
		current->alien_xtalk_model=model;
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			current->alien_xtalk[tone]=alien_xtalk(i,tone);
		}
	}
}
/*
double get_xtalk_gain(int xtalker, int victim, int channel)
{
	return *(channel_matrix + channel + (xtalker * DMTCHANNELS) + (lines * victim * DMTCHANNELS));
}
double get_channel_gain(int line_id, int channel)
{
	return *(channel_matrix + channel + (line_id * DMTCHANNELS) + (lines * line_id * DMTCHANNELS));
}
*/
