#ifndef MG
#define MG

#include "multiuser_load.h"
#include "psd_vector.h"
#include <cstring>

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

class multiuser_greedy{

	public:
	multiuser_greedy();
	virtual ~multiuser_greedy();
	virtual int run(void);
	int *b_target;

	private:
	int **F;
	//double **delta_p;
	double **cost;
	double **normalised_cost;
	double ***_delta_p;
	int **b;
	int *b_total;
	int total_bits;
	double **p;
	double *p_used;
	bool all_tones_full;
	double p_budget;
	double p_error;
	double spectral_mask;
	bool spectral_mask_on;
	bool calc_delta_p_debug;
	bool min_cost_debug;
	bool cost_function_debug;
	bool *user_frozen;
	
	psd_vector *psd;	

	virtual double calc_delta_p_cost(int,int);
	virtual int min_cost(int*,int*);
	bool all_lines_under_budget();
	void init_lines();
	virtual void init_cost_matrix(void);
	virtual void set_loading_algo();
	double power_used(int);
	virtual int update_power(int,int);
	virtual bool total_power_constraint_broken(int,int);
	virtual bool spectral_mask_constraint_broken(int,int);
	void freeze_user(int);
	void print_cost_matrix(int,int);
	void calc_normalised_cost_matrix();
	virtual void reset_data()
	{

		all_tones_full=false;
		total_bits=0;
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			for (int user=0;user<lines;user++) {
				if (tone==0) {
					b_total[user]=0;
					p_used[user]=0;
				}
				F[tone][user]=0;
				cost[tone][user]=0;
				b[tone][user]=0;
				p[tone][user]=0;
				for (int user1=0;user1<lines;user1++) {
					_delta_p[tone][user][user1]=0;
				}
			}
		}
	}

	friend class multiuser_weighted;
	friend class multiuser_new;
	friend class multiuser_fairer;
	friend class bbl_multiuser_greedy;
	friend class bbl_mg_dual;
};


class bbl_multiuser_greedy: public multiuser_greedy
{
	
	friend class bbl_mg_dual;
};

/*
class multiuser_weighted: public multiuser_greedy
{
	public:
	multiuser_weighted();
	double calc_delta_p_cost(int,int);	
	double weighting(int);	
	double flanagans_constant;

	private:
	void set_loading_algo();
	double cost_function(int,int,double*,double*,int*);
	double cost_function1(int,int,double*,double*,int*,double*);
};
*/


#endif
