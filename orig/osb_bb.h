#ifndef OSB_BB
#define OSB_BB

#include "multiuser_load.h"
#include "psd_vector.h"
#include <cfloat>

//#define N_THREADS 1

class region {
	
	public:

	region() {
		b_max = new int[lines];
		b_min = new int[lines];
		next = NULL;
		last = NULL;
	}

	~region() {
		//printf("Region deleted\n");
		//this->print_region();
		delete[] b_max;
		delete[] b_min;
	}

	void print_region() {
		printf("[");
		for (int user=0;user<lines;user++) {
			printf("%d-%d,",this->b_min[user],this->b_max[user]);
		}
		printf("]\n");
		printf("Upper bound = %.50lf\nLower bound = %.50lf\n",upper,lower);
	}

	void calc_upper_bound(double *w,double *l,double mask,double *p_min) {
		double b_sum=0.0;
		double p_sum=0.0;

		for (int user=0;user<lines;user++) {
			b_sum+=w[user]*b_max[user];
		}

		for (int user=0;user<lines;user++) {
			if (p_min[user] < 0 || p_min[user] > mask) {
				upper = -DBL_MAX;
				return;
			}
			p_sum+=l[user]*p_min[user];
		}

		upper=b_sum-p_sum;

	} 
	
	void calc_lower_bound(double *w,double *l,double mask,double *p_max) {
		double b_sum=0.0;
		double p_sum=0.0;

		for (int user=0;user<lines;user++) {
			b_sum+=w[user]*b_min[user];
		}

		for (int user=0;user<lines;user++) {
			if (p_max[user] < 0 || p_max[user] > mask) {
				lower = -DBL_MAX;
				return;
			}
			p_sum+=l[user]*p_max[user];
		}

		lower=b_sum-p_sum;

	}

	static bool region_is_within_region(const region * inner, const region * outer) {
		for (int i=0;i<lines;i++) {
			if ((outer->b_max[i] >= inner->b_max[i]) && (outer->b_min[i] <= inner->b_min[i]))
				;
			else
				return false;	
		}
		return true;
	}

	static void copy(region *dest,region *src) {
		for (int user=0;user<lines;user++) {
			dest->b_max[user] = src->b_max[user];
			dest->b_min[user] = src->b_min[user];
		}
		dest->upper=src->upper;
		dest->lower=src->lower;
	}
	
	int *b_max;
	int *b_min;
	double upper;
	double lower;

	region *next;
	region *last;

};

class regions {
	public:
	int n;
	region *r_list;

	regions(){
		n=1;
		r_list = new region;
		for (int user=0;user<lines;user++) { 		// intial region = entire set
			r_list->b_max[user]=MAXBITSPERTONE;
			r_list->b_min[user]=0;
		}
		r_list->upper=1;
		r_list->lower=-1;
	}
	~regions() {
		delete r_list;
	}
};


class osb_bb {

	public:
	osb_bb(const int);
	~osb_bb();

	int run();
	int bisect_l();

	void branch_regions(regions *,int,psd_vector *);
	void per_tone_exhaustive_search(int tone,psd_vector *psd);
	
	psd_vector** _psd;
	//psd_vector* _psd1;
	
	double _p_tol;
	int _rate_tol;
	double *_w;
	double *_w_max;
	double *_w_min;
	double _min_sl;
	bool _dynamic_lambda;
	int *_rate_targets;
	double _sw;
	double *_p_budget;
	bool _branch_and_bound;
	int _num_threads;
	
	int **_b;
	double **_p;
	
	private:
	void optimise_p();
	bool powers_converged();
	bool rates_converged();
	void update_l();
	void update_w();
	void update_sl();
	void init_lines();
	bool all_powers_within_tol();
	bool all_rates_within_tol();
	double tot_pow(int);
	int tot_rate(int);
	double calc_p_distance();

	region *find_initial_max_low(int tone,psd_vector* psd);
	void isb_solution(int tone,psd_vector *psd,int *b);
	double l_k(int *,int,psd_vector *);

	//region *add_to_temp_list(region*,region*);
	region *del_from_temp_list(region*);
	void print_temp_list(region *);
	int count_temp_list(region *);
	bool check_region_size(region *);



	FILE *_log;

	double *_l;
	double *_best_l;

	int _l_evals;


	double _p_distance;
	double _prev_p_distance;
	double _min_p_distance;

	double _sl;
	double _spectral_mask;

	friend class regions;
	friend class region;

};



#endif 
