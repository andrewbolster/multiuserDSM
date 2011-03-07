#ifndef MIPB_FRAC_H
#define MIPB_FRAC_H

#include "psd_vector.h"

class mipb_frac {

	public:
	mipb_frac();
	~mipb_frac();
	
	int run();
	double *_p_budget;

	bool _spectral_mask_on;
	double _spectral_mask;
	
	bool _graph_loading;
	bool _greedy;

	double _bit_inc;
	double *_wxt;
	double *_w;	// per user static weighting
	double _rate_tol;

	int *_rate_targets;

	private:

	double *_wp;	// per user dynamic power update weight
	double **_cost;
	double ***_delta_p;
	double **_b;
	double **_p;
	double *_p_used;
	double _p_ave;
	double _last_bit_cost;
	double *_b_total;
	int **_F;
	double _total_bits;
	bool _all_tones_full;

	double *_w_min;
	double *_w_max;

	//psd_vector *_psd;

	void init_cost_matrix(void);
	bool total_power_constraint_broken(int,int);
	bool spectral_mask_constraint_broken(int,int);
	int min_cost(int*,int*);
	int update_power(int,int);
	int update_wp();
	double cf(int,int);
	double var_cf(int,int);
	int calc_delta_p(int,int);
	void init_lines(void);
	void write_current_stats(int,int);
	void write_stats_totals();

	void dump_cost_matrix();
	void load();

	void reset_data();

	FILE *_wp_graph;

};


#endif
