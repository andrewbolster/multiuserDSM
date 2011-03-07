#ifndef ISB3g_FRAC
#define ISB3g_FRAC

#include "multiuser_load.h"
#include "psd_vector.h"
#include <cfloat>

//#define N_THREADS 1


class isb3g_frac {

	public:
	isb3g_frac();
	~isb3g_frac();

	int run();
	int bisect_l();

	void per_tone_optimise_p(int tone);	
	
	//psd_vector** _psd;
	
	double _p_tol;
	double _spectral_mask;
	double *_w;
	double _min_sl;
	bool _dynamic_lambda;
	double _bit_inc;
	
	private:
	void optimise_p();
	bool converged();
	void update_l();
	void update_sl();
	void init_lines();
	bool all_powers_within_tol();
	double tot_pow(int);
	double calc_p_distance();

	double l_k(double **b,int tone);

	FILE *_log;
	double **_b;
	double **_p;
	double **_g;

	double *_l;
	double *_best_l;

	int _l_evals;

	double *_p_budget;

	double _p_distance;
	double _prev_p_distance;
	double _min_p_distance;

	double _sl;


};



#endif 
