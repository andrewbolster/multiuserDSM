#ifndef MULTIUSER_WEIGHTED_H
#define MULTIUSER_WEIGHTED_H

#include "multiuser_greedy.h"
#include "psd_vector.h"

/*
struct xtalker_metric
{
	int line_id;
	double lost_bits;
	double percentage;

	friend bool operator<(const xtalker_metric &left,const xtalker_metric &right)
        {
                return left.lost_bits < right.lost_bits;
        }

        friend bool operator>(const xtalker_metric &left,const xtalker_metric &right)
        {
                return left.lost_bits > right.lost_bits;
        }
	


};
*/


class multiuser_weighted: public multiuser_greedy
{
	public:
	multiuser_weighted();
	~multiuser_weighted();
	int run();
	int run_a();
	double calc_delta_p_cost(int,int);	
	int min_cost(int *,int *);
	double weighting(int);	
	double **flanagans_constant;
	double ***_lost_bits;
	double **_actual_lost_bits;
	struct xtalker_metric **_xtalker_table;
	double *w;
	double *w_max;
	double *w_min;
	bool graph_loading;
	double *a;
	char cf_string0[50];
	char cf_string1[50];
	double *linear_grad;
	double *linear_int;
	bool cost_function_debug;
	double ***linear_xt_coeffs;
	bool *user_re_init;
	double *ave_bits;
	double *ave_gain;
	double **bits;
	double **ave_xtalk;
	double *max_met;
	int iter_count;
	int e;
//	psd_vector *psd;

	bool rate_targets_off;
	bool w_updates_off;


	private:
	void set_loading_algo();
	double cost_function(int,int,double*,double*,int*);
	double cost_function1(int,int,double*,double*,int*,double*);
	double ya_cost_function(int,int,double*,double*,int*,double*);
	double linear_cost_function(int,int,double*,double*,int*,double*);
	double xtalk_only_cost_function(int,int,double*,double*,int*);
	double last_cf(int,int,double *,double*,int*,double*);

	bool all_rates_good();	
	bool finished();	

	bool total_power_constraint_broken(int,int);
	bool spectral_mask_constraint_broken(int,int);

	int update_power(int,int);

	void print_lost_bits();
	void sort_lost_bits();
	int calc_lost_bits();
	void update_w();

	void write_current_stats(int,int);
	void write_mw_stats(int,char *);

	double relative_interline_goodness(int,int);
	double relative_interline_goodness_channel(int,int,int);
	double relative_intraline_goodness(int,int);
	void calc_ave_bits(void);
	void calc_ave_xtalk(void);
	double f(double,int,int,int);


	void reset_data()
	{

		all_tones_full=false;
		total_bits=0;
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			for (int user=0;user<lines;user++) {
				if (tone==0) {
					b_total[user]=0;
					p_used[user]=0;
					user_frozen[user]=false;
					user_re_init[user]=false;
				}
				F[tone][user]=0;
				cost[tone][user]=0;
				b[tone][user]=0;
				p[tone][user]=0;
				for (int user1=0;user1<lines;user1++) {
					if (tone==0) {
						_actual_lost_bits[user][user1]=0;
					}
					_delta_p[tone][user][user1]=0;
				}
			}
		}
	}

	


};

#endif
