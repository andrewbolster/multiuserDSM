#ifndef ISB_H

#define ISB_H 1

#include "psd_vector.h"

class isb {

	public:

	isb();
	~isb();
	virtual int run();
	int *rate_target;
	int e;
	double p_budget;	
	double p_tol;

	double *w;
	
	private:
	
	int *current_rate;
	int *last_rate;
	
	double *current_pow;
	double *last_pow;

	psd_vector *psd;

	double *w_1;
	double *w_2;

	double *s_w;

	double s_l;
	double *sl;

	double *l;
	double *l_min;	
	double *l_max;	

	double *last_psd;

	double **_p;
	int **b;
	double *g;

	void optimise_l(int);
	virtual void optimise_p();
	double l_k(int **,int);
	int rate(int);
	double tot_pow(int);
	virtual void init_lines(void);

	void update_w();
	void update_l();
	bool all_rates_good();	
	bool all_powers_good();	
	bool all_powers_within_tol();	
	bool converged();	
	bool rate_oscillating(int);
	bool pow_oscillating(int);
	
	double calc_p_distance();
	
	void print_weights() 
	{
		for(int user=0;user<lines;user++) {
			printf("w[%d] = %8.6e\n",user,w[user]);
		}

		printf("\n");
		
	}

	void print_ls()
	{
		for (int user=0;user<lines;user++) {
			printf("l[%d] = %14.12e\n",user,l[user]);
		}
		
		printf("\n");
	}

	friend class osb;


};

#endif
