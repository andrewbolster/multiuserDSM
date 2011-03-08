#include "multiuser_load.h"
#include "isb3g_frac.h"
#include <cassert>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include "util.h"
static double per_tone_opt_time=0;
static double total_opt_time=0;
static double n_opt=0;
static bool sol_debug=false;
static bool list_debug=false;
static void *opt_p_top_half(void *p);
static void *opt_p_bottom_half(void *p);
static void *opt_p_thread(void * p);
struct thread_args{
	isb3g_frac *o;
	int lower;
	int upper;	
	int thread_id;
};
isb3g_frac::isb3g_frac()
{
	_log = fopen("isb3g_frac.log","w");
	if (_log == NULL) {
		printf("Cannot open log file\n");
		exit(2);
	}
	_b = new double*[DMTCHANNELS];
	_p = new double*[DMTCHANNELS];
	_g = new double*[DMTCHANNELS];
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		_b[tone] = new double[lines];
		_p[tone] = new double[lines];
		_g[tone] = new double[lines];
	}
	_l = new double[lines];
	_best_l = new double[lines];
	_w = new double[lines];
	_p_distance=DBL_MAX;
	_prev_p_distance=DBL_MAX;
	_min_p_distance=DBL_MAX;
	_l_evals=0;
	_p_budget=new double[lines];
	for (int user=0;user<lines;user++) {
		_l[user]=0;
		_w[user]=1;
		_p_budget[user]=0.110;
	}
	//_psd = new psd_vector*[ISB_THREADS];
	/*
	for (int n=0;n<ISB_THREADS;n++) {
		_psd[n] = new psd_vector;
	}*/
	_min_sl=500;
	_p_tol=0.015;
	_dynamic_lambda=false;
	_spectral_mask=dbmhz_to_watts(-20);
	_bit_inc=1;
	// initialisation of _b and _p to shut up valgrind
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			_b[tone][user] = 0;
			_p[tone][user] = 0;
			_g[tone][user] = 9.95;
		}
	}
}
isb3g_frac::~isb3g_frac()
{
	fclose(_log);
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		delete[] _b[tone];
		delete[] _p[tone];
	}
	delete[] _b;
	delete[] _p;
	delete[] _l;
	delete[] _w;
	delete[] _p_budget;
	/*
	for (int n=0;n<ISB_THREADS;n++) {
		delete _psd[n];
	}
	delete[] _psd;
	*/
}
int isb3g_frac::run()
{
	if (_dynamic_lambda)
		_min_sl=1;
	_sl=_min_sl;		// set step size in here in case we changed the default value
	while(!converged() || _l_evals == 0) {
		print_vector(_l,"l");
		//print_vector_to_file(_l,"l",_log);
		optimise_p();
		if (_dynamic_lambda)
			update_sl();
		_l_evals++;
		//fprintf(_log,"Number of lambda evaluations = %d\n",_l_evals);
		update_l();
		//getchar();
		printf("Current distance = %.15lf\n",calc_p_distance());
	}
	printf("Number of lambda evaluations = %d\n",_l_evals);
	init_lines();
	calculate_snr();
	return 0;
}
int isb3g_frac::bisect_l()
{
	double *l_last = new double[lines];
	double *p_current = new double[lines];
	while (1) {	
		for (int user=0;user<lines;user++) {
			double l_min;
			double l_max;
			double pow,last;
			_l[user]=0;
			while(1) {
				//print_vector(_l,"_l");
				optimise_p();
				if (tot_pow(user) > _p_budget[user]) {
					if (_l[user] ==0) {
						_l[user]++;
					}
					else {
						_l[user]*=2;
					}
				}
				else 
					break;
			} 
			//printf("l_max[%d] = %lf\n",user,_l[user]);
			l_max=_l[user];
			l_min=_l[user]/2;
			while(1) {
				_l[user] = (l_max + l_min)/2;
				//print_vector(_l,"_l");
				optimise_p();
				pow = tot_pow(user);
				if (pow > _p_budget[user]) {
					l_min = _l[user];
				}
				else {
					l_max = _l[user];
				}
				if (pow==last) {
					//printf("we're done on user %d\n",user);
					break;
				}
				last=pow;
			}
		}
		print_vector(_l,"_l");
		print_vector(l_last,"l_last");
		for (int user=0;user<lines;user++) {
			p_current[user] = tot_pow(user);
		}
		print_vector(p_current,"p_current");
		bool done=true;
		for (int user=0;user<lines;user++) {
			if (l_last[user] != _l[user]) {
				done=false;
			}
		}
		if (done) {
			printf("The whole lot has converged\n");
			break;
		}
		else {
			printf("Still going\n");
		}
		memcpy(l_last,_l,sizeof(double)*lines);
		if (all_powers_within_tol()) {
			break;
		}
	}
	init_lines();
	calculate_snr();
	return 0;
}
void isb3g_frac::update_sl()
{
	static double min_p_distance_last;
	_p_distance = calc_p_distance();
	if (_p_distance <= _prev_p_distance) {		// new distance is better than last one
		printf("Last distance was %lf, current distance is %lf\n",_prev_p_distance,_p_distance);
		if (_p_distance <= _min_p_distance) {	// new distance is best so far
			printf("Min distance was %lf, new min distance is %lf\n",_min_p_distance,_p_distance);
			_min_p_distance = _p_distance;
			for (int user=0;user<lines;user++) {
				_best_l[user] = _l[user];	// record best l vector
			}
		}
		_sl*=2;
		printf("Step size increased to %lf\n",_sl);
	}
	else {	
		if (_min_p_distance == min_p_distance_last) {
			printf("Looks like we are stuck\n");
			_bit_inc/=2;
		}
		min_p_distance_last = _min_p_distance;
		printf("Starting a new trajectory\n");
		printf("_sl = %lf\n",_sl);
		for (int user=0;user<lines;user++) {
			_l[user] = _best_l[user];
		}
		_sl = _min_sl;
		_prev_p_distance = DBL_MAX;
		return;
	}	
	_prev_p_distance = _p_distance;
	return;
}
void isb3g_frac::optimise_p()
{
	struct timeval time_a,time_b;
	gettimeofday(&time_a,NULL);
	pthread_t th[ISB_THREADS];
	struct thread_args args[ISB_THREADS];
	for (int t=0;t<ISB_THREADS;t++) {
		args[t].o = this;
		args[t].lower = t *(DMTCHANNELS/ISB_THREADS);
		args[t].upper = (t+1) *(DMTCHANNELS/ISB_THREADS);
		args[t].thread_id=t;
		pthread_create(&(th[t]),NULL,opt_p_thread,(void *)&(args[t]));
	}
	for (int t=0;t<ISB_THREADS;t++) {
		pthread_join(th[t],NULL);
	}
	gettimeofday(&time_b,NULL);
	double t = ((time_b.tv_sec - time_a.tv_sec)*1e6 + (time_b.tv_usec - time_a.tv_usec))/1e6;
	printf("average time taken for optimise_p = %lf seconds\n",t/DMTCHANNELS);
}
void *opt_p_thread(void * p)
{
	thread_args *th = static_cast<thread_args *>(p);
	isb3g_frac *o = static_cast<isb3g_frac *>(th->o);	
	printf("Hello from thread %d!\n",th->thread_id);	
	//printf("My lower bound is %d\n",th->lower);	
	//printf("My upper bound is %d\n",th->upper);
	//printf("My psd pointer is %p\n",o->_psd[th->thread_id]);
	//getchar();	
	for (int tone=th->lower;tone<th->upper;tone++) {
		o->per_tone_optimise_p(tone);	
	}
	return (void *)NULL;	// silence warning
}
void isb3g_frac::per_tone_optimise_p(int tone)
{
	double lk,lk_max;
	bool converged;
	double b_max;
	//if (isb_debug) {
	//	print_weights();
	//	print_ls();
	//}
	double *b_last = new double [lines];
	for (int user=0;user<lines;user++) {
		_b[tone][user]=0;
	}
	while(1) {
		for (int user=0;user<lines;user++) {
			b_last[user]=_b[tone][user];
		}	
		for (int user=0;user<lines;user++) {
			for (_b[tone][user]=0,lk=0.0,lk_max=-DBL_MAX;_b[tone][user]<=MAXBITSPERTONE;_b[tone][user]+=_bit_inc) {
				lk = l_k(_b,tone);
				if (lk>lk_max) {
					/*if (lk_debug) {
						printf("New lk_max %4.2lf on tone %d found at\n",lk,tone);
						for (int user1=0;user1<lines;user1++) {
							printf("%d ",b[tone][user1]);
						}
						printf("\n");
						getchar();
					}*/
					lk_max=lk;
					b_max=_b[tone][user];
				}
			}
			_b[tone][user]=b_max;
		}
		converged=true;
		for (int user=0;user<lines;user++) {
			if (b_last[user]!=_b[tone][user]) {
				converged=false;
				break;
			}
		}
		if (converged) {
			double *p = new double [lines];
			//printf("Convergence!\n");
			calculate_psd_vector(b_last,_g[tone],tone,p);
			//psd->calc(b_last,_g[tone],tone,p);
			for (int user=0;user<lines;user++) 
				_p[tone][user]=p[user];
			break;
			delete[] p;
		}
	}
	delete[] b_last;
	//printf("optimise p returned\n");
}
double isb3g_frac::l_k(double **_b,int tone)
{
	double b_sum=0,p_sum=0;
	double *p = new double [lines];
	double *b = new double [lines];
	for (int user=0;user<lines;user++) {
		b_sum+=_b[tone][user]*_w[user];
		b[user]=_b[tone][user];
	}
	calculate_psd_vector(b,_g[tone],tone,p);
	//psd->calc(b,_g[tone],tone,p);
	for (int user=0;user<lines;user++) {
		if (p[user] < 0 || p[user] > _spectral_mask) {
			//p[user]=DBL_MAX;
			delete[] p;
			delete[] b;
			return -DBL_MAX;
		}
		p_sum+=_l[user]*p[user];
		//last_psd[user]=p[user];
	}
	delete[] p;
	delete[] b;
	return b_sum-p_sum;
}
void isb3g_frac::update_l()
{
	for (int user=0;user<lines;user++) {
                double pow = tot_pow(user);
                double update = _sl*(pow-_p_budget[user]);
                if (_l[user] + update < 0)
                        _l[user] = _l[user]*0.1;
		else
			_l[user] = _l[user] + update;
	}
}
bool isb3g_frac::converged()
{
	if (all_powers_within_tol())
                return true;
        else
                return false;
}
void isb3g_frac::init_lines()
{
        struct line *current;
        for (int user=0;user<lines;user++) {
                current=get_line(user);
                current->is_dual=0;
                current->is_frac=true;
                strcpy(current->loading_algo,"ISB3gfrac");
                for (int tone=0;tone<DMTCHANNELS;tone++) {
                        current->b[tone] = _b[tone][user];
                        current->_b[tone] = _b[tone][user];
                        current->psd[tone] = watts_to_dbmhz(_p[tone][user]);
                }
        }
}
bool isb3g_frac::all_powers_within_tol()
{
 	for (int user=0;user<lines;user++) {
		if (_l[user] < 1e-300) {	// if l equals 0 it is a non active power constraint
			continue;
		}
        	if (fabs(tot_pow(user) - _p_budget[user]) > _p_tol) {
			return false;
		}
	}
	return true;
}
double isb3g_frac::tot_pow(int user)
{
	double sum=0;
        for (int tone=0;tone<DMTCHANNELS;tone++) {
                sum+=_p[tone][user];
        }
        printf("power used by user %d = %16.14lf\n",user,sum);
        fprintf(_log,"power used by user %d = %16.14lf\n",user,sum);
        return sum;
}
double isb3g_frac::calc_p_distance()
{
        double sum=0.0;
        for (int user=0;user<lines;user++) {
                sum += pow((_p_budget[user] - tot_pow(user)),2);
        }
	printf("Current p_distance is %lf\n",sqrt(sum));
        return sqrt(sum);
}
