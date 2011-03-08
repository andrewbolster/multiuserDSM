#include "multiuser_load.h"
#include "util.h"
#include "osb_new.h"
#include <cfloat>
static double min_step=1000;
double min_w_step=1e-4;
osb::osb()
{
	/*for (int user=0;user<lines;user++) 
		l[user]=0;
	*/
	s_w=min_w_step;
	mask=dbmhz_to_watts(-36.5);
	p_tol=0.005;
}
int osb::run()
{
	s_l=min_step;
	double prev_p_distance = DBL_MAX;
	//double min_p_distance = DBL_MAX;
	double p_distance = DBL_MAX;
	double prev_w_distance = DBL_MAX;
	//double min_w_distance = DBL_MAX;
	double w_distance = DBL_MAX;
	//double *best_l = new double [lines];
	//double *best_w = new double [lines];
	//bool first_run=true;
	bool reset_l=true;
	bool reset_w=true;
	while ( !(converged()) ) {		
		//psd->print_cache_size();		
		prev_p_distance = p_distance;	
		prev_w_distance = w_distance;
		if (reset_l) {
			prev_p_distance=DBL_MAX;
		}
		reset_l=false;
		if (reset_w) {
			prev_w_distance=DBL_MAX;
		}
		reset_w=false;
		for (int user=0;user<lines;user++) {
			last_rate[user]=current_rate[user];
			last_pow[user]=current_pow[user];
		}
		optimise_p();
		for (int user=0;user<lines;user++) {
			current_rate[user]=rate(user);
			current_pow[user]=tot_pow(user);
		}
		p_distance = calc_p_distance();
		printf("previous l distance was %lf\n",prev_p_distance);
		printf("current l distance is %lf\n",p_distance);
		printf("current  lstep is %lf\n",s_l);
		update_l();
/*	
		if (p_distance <= prev_p_distance) {
			if (p_distance <= min_p_distance) {
				min_p_distance=p_distance;	
				for (int user=0;user<lines;user++)
					best_l[user]=l[user];
			}
			update_l();
			s_l *= 2;
			printf("l step size increased to %lf\n",s_l);
		}
		else { 
			printf("Starting a new l trajectory\n");
			for (int user=0;user<lines;user++) 
				l[user] = best_l[user];
			p_distance=min_p_distance;
			s_l = min_step;
			reset_l=true;
			continue;
		}
*/
/*
		w_distance = calc_w_distance();
		printf("previous distance was %lf\n",prev_w_distance);
		printf("current distance is %lf\n",w_distance);
		printf("current step is %lf\n",s_w);
		if (w_distance <= prev_w_distance) {
			if (w_distance <= min_w_distance) {
				min_w_distance=w_distance;	
				for (int user=1;user<lines;user++)
					best_w[user]=w[user];
			}
			update_w();
			s_w *= 2;
			printf("w step size increased to %lf\n",s_w);
		}
		else { 
			printf("Starting a new w trajectory\n");
			for (int user=1;user<lines;user++) 
				w[user] = best_w[user];
			w_distance=min_w_distance;
			s_w = min_w_step;
			reset_w=true;
		}
*/
		/*
		if (!first_run) {
			int osc=0;
			for (int user=0;user<lines;user++) {
				if (pow_oscillating(user)) {
					osc++;
				}
				//printf("user %d's power is oscillating\n",user);
				if (osc==6) {
					s_l/=2;
					printf("Reducing step size to %lf\n",s_l);
					break;
				}
			}
		}
		*/
		//update_l();
		//update_w();
		//
		//getchar();
		//first_run=false;
	}
	init_lines();
	calculate_snr();
	return 0;
}
double osb::calc_w_distance()
{
	double sum=0.0;
	for (int user=1;user<lines;user++) {
		sum += pow((rate_target[user] - current_rate[user]),2);
	}
	return sqrt(sum);
}
void osb::update_w()
{
	for (int user=1;user<lines;user++) {
		double update = s_w*(rate_target[user]-current_rate[user]);
		/*if ( is_oscillating(user) && abs(current_rate[user] - rate_target[user]) > e) {
			printf("looks like w%d is oscillating\n",user);
			s_w[user] /= 2;
			//getchar();
		}*/
		if (w[user] + update < 0)
			w[user] = w[user];
		else
			w[user] = w[user] + update;
	}
}
void osb::optimise_p()
{
	double lk,lk_max;
	int * b_max = new int [lines];
	int * b_vec = new int [lines];
	double * p = new double [lines];
	print_vector(l,"l");
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		lk=0.0;
		lk_max=-DBL_MAX;	
		vector_next(b_vec,lines,MAXBITSPERTONE,true);	// reset b_vec
		//printf("Starting per tone op\n");
		while(vector_next(b_vec,lines,MAXBITSPERTONE,false)) {
			lk = osb_lk(b_vec,tone);
			if (lk > lk_max) {
				lk_max=lk;
				for (int user=0;user<lines;user++) {
					b_max[user] = b_vec[user];
				}
			}	
		}
		//printf("Finished one tone\n");
        	psd->calc(b_max,g,tone,p);
		for (int user=0;user<lines;user++) {
			b[tone][user] = b_max[user];	
			_p[tone][user] = p[user];		
		}
	}
	delete[] b_max;
	delete[] b_vec;
	delete[] p;
}
double osb::osb_lk(int *b,int tone)
{
        double b_sum=0,p_sum=0;
        double *p = new double [lines];
        for (int user=0;user<lines;user++) {
                b_sum+=b[user]*w[user];
        }
        psd->calc(b,g,tone,p);
        for (int user=0;user<lines;user++) {
                if (p[user] < 0 && b[user] > 0 || p[user] > mask) {
			delete[] p;
                        return -DBL_MAX;
		}
                p_sum+=l[user]*p[user];
                last_psd[user]=p[user];
        }
        delete[] p;
        return b_sum-p_sum;
}
void osb::init_lines()
{
        struct line *current;
        for (int user=0;user<lines;user++) {
                current=get_line(user);
                current->is_dual=0;
                strcpy(current->loading_algo,"OSB");
                for (int tone=0;tone<DMTCHANNELS;tone++) {
                        current->b[tone] = b[tone][user];
                        current->psd[tone] = watts_to_dbmhz(_p[tone][user]);
                }
        }
}
