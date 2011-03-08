#include "multiuser_load.h"
#include "isb3g.h"
#include <cassert>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include "util.h"
#include "psd_cache.h"
#include <boost/threadpool.hpp>
#include <boost/bind.hpp>
#include <boost/mem_fn.hpp>
#include <boost/bind/mem_fn.hpp>
using namespace boost::threadpool;
static double per_tone_opt_time=0;
static double total_opt_time=0;
static double n_opt=0;
static bool sol_debug=false;
static bool list_debug=false;
static void *opt_p_top_half(void *p);
static void *opt_p_bottom_half(void *p);
static void *opt_p_thread(void * p);
struct thread_args{
	isb3g *o;
	int lower;
	int upper;	
	int thread_id;
	bool done;
};
isb3g::isb3g(const int t)
{
	_log = fopen("isb3g.log","w");
	if (_log == NULL) {
		printf("Cannot open log file\n");
		exit(2);
	}
	_b = new int*[DMTCHANNELS];
	_p = new double*[DMTCHANNELS];
	_g = new double*[DMTCHANNELS];
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		_b[tone] = new int[lines];
		_p[tone] = new double[lines];
		_g[tone] = new double[lines];
	}
	_l = new double[lines];
	_best_l = new double[lines];
	_w = new double[lines];
        _w_max = new double[lines];
        _w_min = new double[lines];
	_bisect_completed = new bool[lines];
	_p_distance=DBL_MAX;
	_prev_p_distance=DBL_MAX;
	_min_p_distance=DBL_MAX;
	_l_evals=0;
	_p_budget=new double[lines];
	_rate_targets = new int[lines];
	for (int user=0;user<lines;user++) {
		_l[user]=0;
		_w[user]=1/(double)lines;
		_w_max[user]=100;
		_w_min[user]=0;
		_p_budget[user]=0.110;
		_rate_targets[user]=NOT_SET;
	}
	_num_threads=t;
	/*
	_psd = new psd_vector*[_num_threads];
	for (int n=0;n<_num_threads;n++) {
		_psd[n] = new psd_vector;
	}
	*/
	_min_sl=500;
	_p_tol=0.015;
	_dynamic_lambda=false;
	_spectral_mask=dbmhz_to_watts(50);
	// initialisation of _b and _p to shut up valgrind
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			_b[tone][user] = 0;
			_p[tone][user] = 0;
			_g[tone][user] = 9.95;
		}
	}
	pthread_mutex_init(&_tone_ref_lock,NULL);
	_threadpool=true;
	_rate_tol=10;
}
isb3g::~isb3g()
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
	for (int n=0;n<_num_threads;n++) {
		delete _psd[n];
	}
	delete[] _psd;
	*/
}
int isb3g::run()
{
/*
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
*/
	pool tp(_num_threads);
	bool target = false;
	for (int user=0;user<lines;user++) {
                if (_rate_targets[user] != NOT_SET)
                        target=true;
        }
        if (target) {
		do {
			for(int user=0;user<lines;user++) {
                                if (_rate_targets[user] == NOT_SET)
                                        continue;
				printf("\n\n");
				printf("Bisecting rate on line %d\n",user);
				//getchar();
				int rates[lines];
				int rate_delta[lines];
				for (int user=0;user<lines;user++) {
					rates[user] = tot_rate(user);	
					rate_delta[user] = rates[user] - _rate_targets[user];
				}					
				//print_vector(_w,"w");
				//print_vector(_l,"l");
				//print_vector(rates,"rates");
				//print_vector(rate_delta,"rate_delta");
				//_w_max[user]=100000;	// reset for next bisection
				//_w_min[user]=0;
				bool converged_early=false;	
				bisect_l(&tp);
				printf("%lf\t%d\n",_w[user],tot_rate(user));	
				if (abs(tot_rate(user) - _rate_targets[user]) < _rate_tol) {
					converged_early=true;
					printf("Converged early on line %d\n",user);
					continue;
				}
				if (tot_rate(user) > _rate_targets[user]) {	// max found, now find min	
					printf("First try was max, now finding min\n");
					_w_max[user] = _w[user];
					_w[user]/=2;
					while (1) {
						bisect_l(&tp);
						printf("%lf\t%d\n",_w[user],tot_rate(user));	
						if (abs(tot_rate(user) - _rate_targets[user]) < _rate_tol) {
							converged_early=true;
							break;
						}
						if (tot_rate(user) > _rate_targets[user]) {	// another new max	
							printf("Found another max\n");
							_w_max[user] = _w[user];
							_w[user]/=2;
						}
						else if (tot_rate(user) < _rate_targets[user]) {	// min found
							printf("Found min, breaking out of loop\n");
							_w_min[user] = _w[user];
							break;
						}
					}
				}
				else if (tot_rate(user) < _rate_targets[user]) {	// min found, now find max
					printf("First try was min, now finding max\n");
					_w_min[user] = _w[user];
					_w[user]*=2;
					while (1) {
						bisect_l(&tp);
						printf("%lf\t%d\n",_w[user],tot_rate(user));	
						if (abs(tot_rate(user) - _rate_targets[user]) < _rate_tol) {
							converged_early=true;
							break;
						}
						if (tot_rate(user) > _rate_targets[user]) {		// max found
							printf("Found max, breaking out of loop\n");
							_w_max[user] = _w[user];
							break;
						}
						else if (tot_rate(user) < _rate_targets[user]) {	// another new min
							_w_min[user] = _w[user];
							_w[user]*=2;
							printf("Found another min\n");
						}
					}
				}
				if (converged_early) {
					printf("Converged early on line %d\n",user);
					continue;
				}
				printf("Starting bisection with max = %lf and min = %lf\n",_w_max[user],_w_min[user]);
				_w[user] = (_w_max[user]+_w_min[user])/2;
                                while(1) {
					int rate = tot_rate(user);
					//printf("\n\n");
					//printf("Before bisect_l:\n");
					//print_vector(_w,"w");
					//printf("Rate[%d] = %d\n",user,rate);
					bisect_l(&tp);
					rate = tot_rate(user);
					//printf("After bisect_l:\n");
					//print_vector(_w,"w");
					//printf("Rate[%d] = %d\n",user,rate);
					printf("%lf\t%d\n",_w[user],tot_rate(user));	
					//printf("\n\n");
					if (tot_rate(user) > _rate_targets[user]) {
                                                _w_max[user] = _w[user];
						//printf("New w_max found = %50.48lf\n",_w_max[user]);
                                        }
                                        else if (tot_rate(user) < _rate_targets[user]) {
                                                _w_min[user] = _w[user];
						//printf("New w_min found = %50.48lf\n",_w_min[user]);
                                        }
                                        if (abs(tot_rate(user) - _rate_targets[user]) < _rate_tol) {
                                                //printf("Rate converged on line %d\n",user);
                                                break;
                                        }
                                        _w[user] = (_w_max[user]+_w_min[user])/2;
				}
				for (int user=0;user<lines;user++) {
					rates[user] = tot_rate(user);	
					rate_delta[user] = rates[user] - _rate_targets[user];
				}					
				//print_vector(_w,"w");
				//print_vector(_l,"l");
				//print_vector(rates,"rates");
				//print_vector(rate_delta,"rate_delta");
				//printf("\n\n");
			}
		}while(!rates_converged());
	}
	else {
		bisect_l(&tp);
	}
	print_vector(_w,"w");
	init_lines();
	calculate_snr();
}
int isb3g::bisect_l(pool *tp)
{
	double *l_last = new double[lines];
	double *p_current = new double[lines];
	//pool tp(_num_threads);
	//pool tp;
	int iters=0;
	while (1) {	
		for (int user=0;user<lines;user++) {
			double l_min;
			double l_max;
			double pow,last;
			_bisect_completed[user]=false;
			/*
			_l[user]=0;
			while(1) {
				print_vector(_l,"_l");
				optimise_p(&tp);
				for (int user=0;user<lines;user++) {
					printf("%lf\t",tot_pow(user));
				}
				printf("\n");
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
			*/
			//print_vector(_l,"_l");
			//
			//printf("Bisecting power on line %d\n",user);
			double temp = _l[user];
			_l[user] = 0;
			optimise_p(tp);
			if (tot_pow(user) < _p_budget[user]) {
				//printf("Non active power constraint on line %d\n",user);
				continue;
			}
			_l[user]=temp;
			optimise_p(tp);
			if (tot_pow(user) > _p_budget[user]) {	// current l is min
				//printf("Found the initial min on line %d = %lf p = %lf\n",user,_l[user],tot_pow(user));
				l_min=_l[user];
				while(1) {		// now find max
					//print_vector(_l,"_l");
					optimise_p(tp);
					//for (int user=0;user<lines;user++) {
					//	printf("%lf\t",tot_pow(user));
					//}		
					//printf("\n");
					if (tot_pow(user) > _p_budget[user]) {
						if (_l[user] ==0) {
							_l[user]++;
						}
						else {
							_l[user]*=2;
						}
					}
					else {
						//printf("Now found the max on line %d = %lf p = %lf\n",user,_l[user],tot_pow(user));
						break;
					}
				}
				l_max=_l[user];
			}
			else {	// current l is max
				//printf("Found the initial max on line %d = %lf p = %lf\n",user,_l[user],tot_pow(user));
				l_max=_l[user];
				while(1) {		// now find min
					//print_vector(_l,"_l");
					optimise_p(tp);	
					//for (int user=0;user<lines;user++) {
					//	printf("%lf\t",tot_pow(user));
					//}	
					//printf("\n");
					if (tot_pow(user) < _p_budget[user]) {
						_l[user]/=2;
						if (_l[user] == 0)
							break;
					}
					else {
						//printf("Now found the min on line %d = %lf p = %lf\n",user,_l[user],tot_pow(user));
						break;
					}
				}
				l_min=_l[user];
			}
			while(1) {
				_l[user] = (l_max + l_min)/2;
				//print_vector(_l,"_l");
				optimise_p(tp);
				//printf("p[%d] = %lf\n",user,tot_pow(user));
				//for (int user=0;user<lines;user++) {
				//	printf("%lf\t",tot_pow(user));
				//}
				//printf("\n");
				pow = tot_pow(user);
				if (pow > _p_budget[user]) {
					l_min = _l[user];
				}
				else {
					l_max = _l[user];
				}
				//if (pow==last || (fabs(_p_budget[user]-pow) < _p_tol)) {
				if (pow==last) {
					//printf("we're done on user %d\n",user);
					_bisect_completed[user]=true;
					break;
				}
				/*	
				if (pow > _p_budget[user]) {
					if ((pow - _p_budget[user]) < _p_tol) {
						for (int user=0;user<lines;user++) {
							p_current[user] = tot_pow(user);
						}
						//print_vector(p_current,"p_current");
						break;
					}
				}
				if (pow < _p_budget[user]) {
					if ((_p_budget[user] - pow) < _p_tol) {
						for (int user=0;user<lines;user++) {
							p_current[user] = tot_pow(user);
						}
						//print_vector(p_current,"p_current");
						break;
					}
				}
				*/
				last=pow;
			}
		}
		iters++;
		//print_vector(_l,"_l");
		//print_vector(l_last,"l_last");
		for (int user=0;user<lines;user++) {
			p_current[user] = tot_pow(user);
		}
		//print_vector(p_current,"p_current");
		bool done=true;
		for (int user=0;user<lines;user++) {
			if (l_last[user] != _l[user]) {
				done=false;
			}
		}
		if (done) {
			//printf("l bisection converged\n");
			break;
		}
		else {
			//printf("l bisection still going\n");
		 	//printf("Current iterations = %d\n",iters);
		}
		if (all_powers_within_tol()) {
			//printf("l bisection didnt converge yet, but all powers are within tolerance\n");
			break;
		}
		memcpy(l_last,_l,sizeof(double)*lines);
	}
	//print_vector(_l,"l");
	for (int user=0;user<lines;user++) {
		p_current[user] = tot_pow(user);
	}
	//print_vector(p_current,"p_current");
	//printf("Iterations for this bisection = %d\n",iters);
	return 0;
}
void isb3g::update_sl()
{
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
		printf("Starting a new trajectory\n");
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
void isb3g::optimise_p(pool *tp)
{
	static int threads_scheduled=0;
	int rc;
	_tone_ref_count=DMTCHANNELS;
	if (_threadpool) {
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			int psd_vec = (int)(tone/(DMTCHANNELS/_num_threads));
			//printf("Scheduling tone %d\n",tone);
			//printf("active = %d pending = %d\n",tp->active(),tp->pending());
			//per_tone_optimise_p(tone,(psd_vector *)NULL);
			//if (tp->schedule(boost::bind(boost::mem_fn(&isb3g::per_tone_optimise_p),this,tone,(psd_vector *)NULL)) == false) {
			if (tp->schedule(boost::bind(boost::mem_fn(&isb3g::per_tone_optimise_p),this,tone,(psd_vector *)NULL)) == false) {
			//if (tp->schedule(boost::bind(boost::mem_fn(&isb3g::per_tone_optimise_p),this,tone,_psd[psd_vec])) == false) {
				printf("Fuck up occured in threadpool\n");
				exit(1);
			}
			//threads_scheduled++;
			//printf("Threads scheduled = %d\n",threads_scheduled);
		}
		//printf("Waiting for threads to terminate\n");
		tp->wait();
		//sleep(5);
		//printf("Threads terminated\n");
		//
		//
		/*	
		while (1) {
			pthread_mutex_lock(&_tone_ref_lock);		
			if (_tone_ref_count<=0) {
				pthread_mutex_unlock(&_tone_ref_lock);		
				break;
			}
			pthread_mutex_unlock(&_tone_ref_lock);		
		}
		*/
		//sleep(1);
	}
	else {
		struct timeval time_a,time_b;
		gettimeofday(&time_a,NULL);
		/*
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);		
		*/
		pthread_t th[_num_threads];
		struct thread_args args[_num_threads];
		for (int t=0;t<_num_threads;t++) {
			args[t].o = this;
			args[t].lower = t *(DMTCHANNELS/_num_threads);
			args[t].upper = (t+1) *(DMTCHANNELS/_num_threads);
			args[t].thread_id=t;
			args[t].done=false;
			rc = pthread_create(&(th[t]),NULL,opt_p_thread,(void *)&(args[t]));
			if (rc) {
				printf("pthread_created returned error code %d\n",rc);
				exit(-1);
			}
			threads_scheduled++;
			//printf("Threads scheduled = %d\n",threads_scheduled);
		}
		for (int t=0;t<_num_threads;t++) {
			rc = pthread_join(th[t],NULL);
			if (rc) {
				printf("pthread_join returned error code %d\n",rc);
				exit(-1);
			}
			//while(args[t].done==false);
			//pthread_detach(th[t]);
		}
		gettimeofday(&time_b,NULL);
		double t = ((time_b.tv_sec - time_a.tv_sec)*1e6 + (time_b.tv_usec - time_a.tv_usec))/1e6;
		printf("average time taken for branch and bound = %lf seconds\n",t/DMTCHANNELS);
	}
}
void *opt_p_thread(void * p)
{
	thread_args *th = static_cast<thread_args *>(p);
	isb3g *o = static_cast<isb3g *>(th->o);	
	printf("Hello from thread %d!\n",th->thread_id);	
	//printf("My lower bound is %d\n",th->lower);	
	//printf("My upper bound is %d\n",th->upper);
	//printf("My psd pointer is %p\n",o->_psd[th->thread_id]);
	//getchar();	
	for (int tone=th->lower;tone<th->upper;tone++) {
		o->per_tone_optimise_p(tone,o->_psd[th->thread_id]);	
	}
	th->done=true;
	//pthread_exit((void *)NULL);
	return (void *)NULL;	// silence warning
}
void isb3g::per_tone_optimise_p(int tone,psd_vector *psd)
{
	double lk,lk_max;
	bool converged;
	int b_max,b_max_2;
	static int dump=0;
	FILE *dump_file;
	//if (isb_debug) {
	//	print_weights();
	//	print_ls();
	//}
	//
	//printf("starting on tone %d\n",tone);
	int *b_last = new int [lines];
	for (int user=0;user<lines;user++) {
		_b[tone][user]=0;
	}
	while(1) {
		for (int user=0;user<lines;user++) {
			b_last[user]=_b[tone][user];
		}	
		for (int user=0;user<lines;user++) {
			for (_b[tone][user]=0,lk=0.0,lk_max=-DBL_MAX;_b[tone][user]<=MAXBITSPERTONE;_b[tone][user]++) {
				lk = l_k(_b,tone,psd);
				if (lk>lk_max) {
					if (0) {
						printf("New lk_max %4.2lf on tone %d found at\n",lk,tone);
						print_vector(_b[tone],"b[tone]");	
						getchar();
					}
					lk_max=lk;
					b_max=_b[tone][user];
				}
				//else {
				//	b_max=_b[tone][user]-1;		// my more effecient search
				//	break;
				//}
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
			//print_vector(b_last,"b_last");
			//calculate_psd_vector(b_last,g,tone,p);
			if (_threadpool)
				calculate_psd_vector(b_last,_g[tone],tone,p,cache);
			else
				psd->calc(b_last,_g[tone],tone,p);
			for (int user=0;user<lines;user++) 
				_p[tone][user]=p[user];
			delete[] p;
			break;
		}
	}
	delete[] b_last;
	/*	
	pthread_mutex_lock(&_tone_ref_lock);	
	_tone_ref_count--;
	pthread_mutex_unlock(&_tone_ref_lock);	
	*/
	//printf("done on tone %d\n",tone);
}
/*
void isb3g::eff_per_tone_opt_p(int tone,psd_vector *psd)
{
	double lk,lk_max;
	bool converged;
	int b_max;
	int *b_last = new int [lines];
	for (int user=0;user<lines;user++) {
		_b[tone][user]=0;
	}
	while(1) {
		for (int user=0;user<lines;user++) {
			b_last[user]=_b[tone][user];
		}	
		for (int user=0;user<lines;user++) {
			//for (_b[tone][user]=0,lk=0.0,lk_max=-DBL_MAX;_b[tone][user]<=MAXBITSPERTONE;_b[tone][user]++) {
			//	lk = l_k(_b,tone,psd);
			//	if (lk>lk_max) {
			//		lk_max=lk;
			//		b_max=_b[tone][user];
			//	}
			//}
			int dir=UP;
			int bit=0;
			int max_bit=MAXBITSPERTONE;
			int min_bit=0;
			_b[tone][user]=bit;
			lk = l_k(_b,tone,psd);
			while(1) {
				if (dir == UP) {
					bit = (bit+max_bit)/2;
				}
				else if (dir == DOWN) {
					bit = bit - (min_bit+bit)/2;
				}
				_b[tone][user]=bit;
	                        lk = l_k(_b,tone,psd);
				if (lk > lk_max) {
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
			//print_vector(b_last,"b_last");
			//calculate_psd_vector(b_last,g,tone,p);
			psd->calc(b_last,_g[tone],tone,p);
			for (int user=0;user<lines;user++) 
				_p[tone][user]=p[user];
			break;
			delete[] p;
		}
	}
	delete[] b_last;
}
*/
double isb3g::l_k(int **_b,int tone,psd_vector *psd)
{
	double b_sum=0,p_sum=0;
	double *p = new double [lines];
	int *b = new int [lines];
	for (int user=0;user<lines;user++) {
		b_sum+=_b[tone][user]*_w[user];
		b[user]=_b[tone][user];
	}
	/*calculate_psd_vector(_b,g,tone,p);*/
	if (_threadpool)
		calculate_psd_vector(b,_g[tone],tone,p,cache);
	else
		psd->calc(b,_g[tone],tone,p);
	for (int user=0;user<lines;user++) {
		if ((p[user] < 0) || (p[user] > _spectral_mask)) {
			/*
			p[user]=DBL_MAX;
			if (p[user] > _spectral_mask) {
				printf("Spectral mask broken\n");
				getchar();
			}
			*/
			delete[] p;
			delete[] b;
			return -DBL_MAX;
		}
		p_sum+=_l[user]*p[user];
		/*last_psd[user]=p[user];*/
	}
	delete[] p;
	delete[] b;
	return b_sum-p_sum;
}
void isb3g::update_l()
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
bool isb3g::converged()
{
	if (all_powers_within_tol())
                return true;
        else
                return false;
}
void isb3g::init_lines()
{
        struct line *current;
        for (int user=0;user<lines;user++) {
                current=get_line(user);
                current->is_dual=0;
                strcpy(current->loading_algo,"ISB3g");
                for (int tone=0;tone<DMTCHANNELS;tone++) {
                        current->b[tone] = _b[tone][user];
                        current->psd[tone] = watts_to_dbmhz(_p[tone][user]);
                }
        }
}
bool isb3g::all_powers_within_tol()
{
 	for (int user=0;user<lines;user++) {
		if (_l[user] < 1e-10) {	// if l equals 0 it is a non active power constraint
			printf("non active power constraint on line %d\n",user);
			continue;
		}
		if (_bisect_completed[user]) {
			printf("Bisect fully completed for line %d p = %lf\n",user,tot_pow(user));
			continue;
		}
        	/*if (fabs(tot_pow(user) - _p_budget[user]) > _p_tol) {
			return false;
		}*/
		double pow=tot_pow(user);
		if (pow > _p_budget[user]) {
			if ((pow - _p_budget[user]) > _p_tol) {
				printf("power greater than tol on line %d = %lf\n",user,pow);
				return false;
			}
		}
		if (pow < _p_budget[user]) {
			if ((_p_budget[user] - pow) > _p_tol) {
				printf("power too low on line %d = %lf\n",user,pow);
				return false;
			}
		}
	}
	printf("All lines within power constraints\n");
	return true;
}
double isb3g::tot_pow(int user)
{
	double sum=0;
        for (int tone=0;tone<DMTCHANNELS;tone++) {
                sum+=_p[tone][user];
        }
        //printf("power used by user %d = %16.14lf\n",user,sum);
        fprintf(_log,"power used by user %d = %16.14lf\n",user,sum);
        return sum;
}
double isb3g::calc_p_distance()
{
        double sum=0.0;
        for (int user=0;user<lines;user++) {
                sum += pow((_p_budget[user] - tot_pow(user)),2);
        }
	printf("Current p_distance is %lf\n",sqrt(sum));
        return sqrt(sum);
}
bool isb3g::rates_converged()
{
        if (all_rates_within_tol())
                return true;
        else
                return false;
}
bool isb3g::all_rates_within_tol()
{
        for (int user=0;user<lines;user++) {
                if (_rate_targets[user] == NOT_SET)
                        continue;
                if (abs(tot_rate(user) - _rate_targets[user]) > _rate_tol) {
			return false;
		}
		/*
                int rate = tot_rate(user);
		if (rate > _rate_targets[user]) {
			printf("rate on line %d was greater than target so its grand\n",user);
			continue;
		}
		else {
			if ((_rate_targets[user] - rate) > _rate_tol) {
                        	printf("Rate on line %d was %d which is not good enough\n",user,rate);
				return false;
			}
                }
		*/
        }
        return true;
}
int isb3g::tot_rate(int user)
{
        int sum=0;
        for (int tone=0;tone<DMTCHANNELS;tone++) {
                sum+=_b[tone][user];
        }
        //printf("Rate of user %d = %d\n",user,sum);
        return sum;
}
