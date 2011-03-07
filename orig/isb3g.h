#ifndef ISB3g
#define ISB3g

#include "multiuser_load.h"
#include "psd_vector.h"
#include <cfloat>
#include <pthread.h>

//#define N_THREADS 1

#include <boost/threadpool.hpp>
#include <boost/bind.hpp>
#include <boost/mem_fn.hpp>
#include <boost/bind/mem_fn.hpp>

using namespace boost::threadpool;


class isb3g {

	public:
	isb3g(const int);
	~isb3g();

	int run();
	int bisect_l(pool *);

	void per_tone_optimise_p(int tone,psd_vector *);	
	
	psd_vector** _psd;
	
	double _p_tol;
	double _spectral_mask;
	double *_w;
	double _min_sl;
	bool _dynamic_lambda;
	
	double *_p_budget;

	int _num_threads;

	bool _threadpool;

	int *_rate_targets;
	int _rate_tol;

        double *_w_max;
        double *_w_min;
	
	
	private:
	void optimise_p(pool *);
	bool converged();
	void update_l();
	void update_sl();
	void init_lines();
	bool all_powers_within_tol();
	double tot_pow(int);
	double calc_p_distance();

	double l_k(int **b,int tone,psd_vector *);

	FILE *_log;
	int **_b;
	double **_p;
	double **_g;

	double *_l;
	double *_best_l;

	int _l_evals;

	bool *_bisect_completed;

	pthread_mutex_t _tone_ref_lock;
	int _tone_ref_count;

	double _p_distance;
	double _prev_p_distance;
	double _min_p_distance;

	double _sl;

	int tot_rate(int);
	bool all_rates_within_tol();
	bool rates_converged();

};



#endif 
