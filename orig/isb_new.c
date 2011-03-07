#include "multiuser_load.h"
#include "isb_new.h"
#include "psd_vector.h"
#include <cmath>
#include <cstdlib>
#include <cfloat>

extern bool isb_debug;
bool lk_debug=false;

double min_step=1e-6;
double min_s_l=1000;

static double mask;

isb::isb()
{

	rate_target = new int [lines];
	current_rate = new int [lines];
	last_rate = new int [lines];

	current_pow = new double [lines];
	last_pow = new double [lines];

	w = new double [lines];
	w_1 = new double [lines];
	w_2 = new double [lines];

	s_w = new double [lines];

	sl = new double [lines];

	l = new double [lines];
	l_min = new double [lines];
	l_max = new double [lines];

	last_psd = new double [lines];

	b = new int *[DMTCHANNELS];
	_p = new double *[DMTCHANNELS];
	g = new double [lines];

	for (int user=0;user<lines;user++) {
		rate_target[user] = NOT_SET;
		g[user]=9.95;
		s_w[user]=0.1;
		sl[user]=min_step;
	}

	w[0]=1;
	
	for (int user=0;user<lines;user++) {
		w[user]=1;
		l[user]=10000;
		//l[user]=0;	
	}

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		b[tone] = new int [lines];
		_p[tone] = new double [lines];
		//calculate_b_vector_flat_psd(dbmhz_to_watts(-40),g,tone,b[tone]);
		
		//printf("channel %d\n",tone);
		for (int user=0;user<lines;user++) {
			//printf("b[%d] = %d\n",user,b[tone][user]);	
			b[tone][user]=0;
		}
		//getchar();
	}

	p_budget=0.110;
	p_tol=0.001;

	s_l=min_s_l;

	psd = new psd_vector;
	//psd->_caching_on=false;

	mask = dbmhz_to_watts(-30);
}

isb::~isb(){

	delete psd;
}

int isb::run()
{

	double prev_p_distance = DBL_MAX;
	double min_p_distance = DBL_MAX;
	double p_distance = DBL_MAX;
	double *best_l = new double [lines];
	bool first_run=true;
	bool reset=true;

	while ( !(converged()) ) {		
		//psd->print_cache_size();		
		prev_p_distance = p_distance;
		
		if (reset) {
			prev_p_distance=DBL_MAX;
		}
		reset=false;

		for (int user=0;user<lines;user++) {
			last_rate[user]=current_rate[user];
			last_pow[user]=current_pow[user];
		}
		optimise_p();
		for (int user=0;user<lines;user++) {
			current_rate[user]=rate(user);
			current_pow[user]=tot_pow(user);
		}
/*
		p_distance = calc_p_distance();
	
		printf("previous distance was %lf\n",prev_p_distance);
		printf("current distance is %lf\n",p_distance);
		printf("current step is %lf\n",s_l);
	
		if (p_distance <= prev_p_distance) {
			if (p_distance <= min_p_distance) {
				min_p_distance=p_distance;	
				for (int user=0;user<lines;user++)
					best_l[user]=l[user];
			}
			update_l();
			s_l *= 2;
			printf("step size increased to %lf\n",s_l);
		}
		else { 
			printf("Starting a new trajectory\n");
			for (int user=0;user<lines;user++) {
				l[user] = best_l[user];
			}
			p_distance = min_p_distance;
			s_l = min_s_l;
			reset=true;
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
		update_l();
		//update_w();
		//
		//getchar();
		//first_run=false;
	}

	init_lines();

	calculate_snr();
	
	return 0;

}

double isb::calc_p_distance()
{

	double sum=0.0;

	for (int user=0;user<lines;user++) {
		sum += pow((p_budget - tot_pow(user)),2);
	}

	return sqrt(sum);

}

bool isb::converged()
{
	//if (all_rates_good() && all_powers_good() && (p_budget - tot_pow(0) < 0.001) && (tot_pow(0) < p_budget))
	//if (all_rates_good() && all_powers_good())
	if (all_powers_within_tol())
		return true;
	else
		return false;
}




void isb::optimise_p()
{

	double lk,lk_max;
	bool converged;
	int b_max;

	if (isb_debug) {
		print_weights();
		print_ls();
	}
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		int *b_last = new int [lines];
		for (int user=0;user<lines;user++) {
			b[tone][user]=0;
		}
		while(1) {
			for (int user=0;user<lines;user++) {
				b_last[user]=b[tone][user];
			}	
			for (int user=0;user<lines;user++) {
				for (b[tone][user]=0,lk=0.0,lk_max=-DBL_MAX;b[tone][user]<=MAXBITSPERTONE;b[tone][user]++) {
					lk = l_k(b,tone);
					if (lk>lk_max) {
						if (lk_debug) {
							printf("New lk_max %4.2lf on tone %d found at\n",lk,tone);
							for (int user1=0;user1<lines;user1++) {
								printf("%d ",b[tone][user1]);
							}
							printf("\n");
							getchar();
						}
						lk_max=lk;
						b_max=b[tone][user];
					}
				}
				b[tone][user]=b_max;
			}
	
			converged=true;
			for (int user=0;user<lines;user++) {
				if (b_last[user]!=b[tone][user]) {
					converged=false;
					break;
				}
			}
	
			if (converged) {
				double *p = new double [lines];
				//printf("Convergence!\n");
				//calculate_psd_vector(b_last,g,tone,p);
				psd->calc(b_last,g,tone,p);
				for (int user=0;user<lines;user++) 
					_p[tone][user]=p[user];
				break;
				delete[] p;
			}

		}
		delete[] b_last;

	}

	//printf("optimise p returned\n");

}


double isb::l_k(int **b,int tone)
{

	double b_sum=0,p_sum=0;
	double *p = new double [lines];
	int *_b = new int [lines];

	for (int user=0;user<lines;user++) {
		b_sum+=b[tone][user]*w[user];
		_b[user]=b[tone][user];
	}


	//calculate_psd_vector(_b,g,tone,p);
	psd->calc(_b,g,tone,p);
	
	for (int user=0;user<lines;user++) {
		if (p[user] < 0 || p[user] > mask) {
			//p[user]=DBL_MAX;
			delete[] p;
			delete[] _b;
			return -DBL_MAX;
		}
		p_sum+=l[user]*p[user];
		last_psd[user]=p[user];
	}

	delete[] p;
	delete[] _b;

	return b_sum-p_sum;

}


void isb::update_w()
{
	double s=0.001;

	for (int user=1;user<lines;user++) {
		double update = s*(rate_target[user]-current_rate[user]);
		w_2[user]=w_1[user];
		w_1[user]=w[user];
		
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

void isb::update_l()
{

	for (int user=0;user<lines;user++) {
		double pow = current_pow[user];
		double update = s_l*(pow-p_budget);
		if (l[user] + update < 0)
			l[user] = l[user]*0.9;
		//	l[user] = 0;
		else {
			//if (fabs(p_budget - pow) < p_tol)
			//	l[user] = l[user] + update*0.5;
			//else
			l[user] = l[user] + update;
		}
	}

}

bool isb::rate_oscillating(int user)
{

	if ((current_rate[user] > rate_target[user] && last_rate[user] < rate_target[user]) || (current_rate[user] < rate_target[user] && last_rate[user] > rate_target[user]))
		return true;
	else
		return false;
}

bool isb::pow_oscillating(int user)
{

	//if ((current_pow[user] > p_budget+p_tol && last_pow[user] < p_budget-p_tol) || (current_pow[user] < p_budget-p_tol && last_pow[user] > p_budget+p_tol))
	if ((current_pow[user] > p_budget && last_pow[user] < p_budget) || (current_pow[user] < p_budget && last_pow[user] > p_budget))
		return true;
	else 
		return false;

}

bool isb::all_rates_good()
{

	rate(0);
	for (int user=1;user<lines;user++) {
		if (abs(rate_target[user] - current_rate[user]) > e) {
			return false;
		}
	}

	return true;

}

bool isb::all_powers_good()
{

	for (int user=0;user<lines;user++) {
		if (tot_pow(user) > p_budget) {
		//if ((p_budget - tot_pow(user)) > 0.01 && (tot_pow(user) < p_budget)) // 10mw
			return false;
		}
	}
	
	return true;

}

bool isb::all_powers_within_tol()
{

	for (int user=0;user<lines;user++) {
		if (l[user] < 1e-20) {
			continue;
		}
		if (fabs(tot_pow(user) - p_budget) > p_tol) {
		//if ((p_budget - tot_pow(user)) > 0.01 && (tot_pow(user) < p_budget)) // 10mw
			return false;
		}
	}
	
	return true;

}


int isb::rate(int user)
{

	int sum=0;

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		sum+=b[tone][user];
	}

	printf("rate of user %d = %d\n",user,sum);
	
	return sum;

}


double isb::tot_pow(int user)
{
	double sum=0;

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		sum+=_p[tone][user];
		
		if (_p[tone][user] < 0) {
			printf("Ooops! p on tone %d user %d is negative\n",tone,user);
			printf("b on tone %d is %d,%d\n",tone,b[tone][0],b[tone][1]);
			exit(0);
		}
		
	}

	printf("power used by user %d = %16.14lf\n",user,sum);
	
	return sum;

}

void isb::init_lines()
{

	struct line *current;

	for (int user=0;user<lines;user++) {
		current=get_line(user);
		current->is_dual=0;
		strcpy(current->loading_algo,"ISB");
		
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			current->b[tone] = b[tone][user];
			current->psd[tone] = watts_to_dbmhz(_p[tone][user]);	
		}
	}

}
/*
int isb::run()
{

	for (int user=1;user<lines;user++) {
		w[user]=(w_max[user]+w_min[user])/2;
	}

	for (int user=1;user<lines;user++) {

		if (isb_debug) {
			printf("Now adjusting w%d\n",user);
		}

		while(abs(rate(user) - rate_target[user]) > e) {
			w[user]=(w_max[user]+w_min[user])/2;
				
			optimise_l(0);

			if (rate(user) > rate_target[user]) {
				if (1) {
					printf("current rate of user %d is higher than rate target %d\n",user,rate_target[user]);
				}
				w_max[user]=w[user];
			}
			else {
				if (1) {
					printf("current rate of user %d is less than rate target %d\n",user,rate_target[user]);
				}
				w_min[user]=w[user];
			}
		}
	}

	init_lines();

	calculate_snr();


}
*/

void isb::optimise_l(int user)
{

	double pow,last;

	l[user]=1.0;
        l_min[user]=0.0;
		
        do {
		if (isb_debug)	{
			printf("Finding the max value of l%d\n",user);
		}

                l[user]*=2;
		if (user+1==lines)
			optimise_p();
		else
                	optimise_l(user+1);
        } while(tot_pow(user) >  p_budget);

        l_max[user]=l[user];

	if (user ==0) {
		printf("Max value of l0 is %lf\n",l[user]);
		printf("l1 = %lf\tl2 = %lf\n",l[1],l[2]);
		printf("power on line 0 is %lf\n",tot_pow(user));
		getchar();
	}

	if (isb_debug) {
		printf("Found max value of l%d = %4.2lf\n",user,l[user]);
	}


	while(1) {
		l[user]=(l_max[user]+l_min[user])/2;
		
		if (user+1==lines)
			optimise_p();
		else
			optimise_l(user+1);

		pow=tot_pow(user);

		if (pow > p_budget) {
			if (isb_debug) {
				printf("Total power on line %d = %6.4lf is greater than budget\n",user,tot_pow(user));
			}
			l_min[user]=l[user];
		} 
		else {
			if (isb_debug) {
				printf("Total power on line %d = %6.4lf is less than budget\n",user,tot_pow(user));
			}
			l_max[user]=l[user];
			
			if (pow==last) {
				if (isb_debug) {	
					printf("optimise l%d returned pow = %6.4lf\n",user,tot_pow(user));
				}
				return;
			}

		}
		last=pow;

	}

}
