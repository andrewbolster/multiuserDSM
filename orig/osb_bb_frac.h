#ifndef OSB_BB_FRAC
#define OSB_BB_FRAC

#include "multiuser_load.h"
#include "psd_vector.h"
#include <cfloat>

//#define N_THREADS 1

class region_frac {
	
	public:

	region_frac() {
		b_max = new double[lines];
		b_min = new double[lines];
		next = NULL;
		last = NULL;
	}

	~region_frac() {
		//printf("Region deleted\n");
		//this->print_region_frac();
		delete[] b_max;
		delete[] b_min;
	}

	void print_region_frac() {
		printf("[");
		for (int user=0;user<lines;user++) {
			printf("%24.22lf-%24.22lf,",this->b_min[user],this->b_max[user]);
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

		if (upper > lines*MAXBITSPERTONE) {
			printf("SRSLY?\n");
			print_vector(p_min,"p_min");
			print_vector(b_max,"b_max");
			print_vector(w,"w");
			print_vector(l,"l");
			getchar();
		}
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
	
	double *b_max;
	double *b_min;
	double upper;
	double lower;

	region_frac *next;
	region_frac *last;

};

class regions_frac {
	public:
	int n;
	region_frac *r_list;

	regions_frac(){
		n=1;
		r_list = new region_frac;
		for (int user=0;user<lines;user++) { 		// intial region_frac = entire set
			r_list->b_max[user]=MAXBITSPERTONE;
			r_list->b_min[user]=0;
		}
		r_list->upper=1;
		r_list->lower=-1;
	}
	~regions_frac() {
		delete r_list;
	}
};


class osb_bb_frac {

	public:
	osb_bb_frac();
	~osb_bb_frac();

	int run();
	int bisect_l();

	void branch_regions(regions_frac *,int);
	
	psd_vector** _psd;
	//psd_vector* _psd1;
	
	double _p_tol;
	double *_w;
	double _min_sl;
	bool _dynamic_lambda;
	double _inc;
	
	private:
	void optimise_p();
	bool converged();
	void update_l();
	void update_sl();
	void init_lines();
	bool all_powers_within_tol();
	double tot_pow(int);
	double calc_p_distance();
	

	//region_frac *add_to_temp_list(region*,region*);
	region_frac *del_from_temp_list(region_frac *);
	void print_temp_list(region_frac *);
	int count_temp_list(region_frac *);
	bool check_region_frac_size(region_frac *);



	FILE *_log;
	double **_b;
	double **_p;

	double *_l;
	double *_best_l;

	int _l_evals;

	double *_p_budget;

	double _p_distance;
	double _prev_p_distance;
	double _min_p_distance;

	double _sl;
	double _spectral_mask;

	friend class regions_frac;
	friend class region_frac;

};



#endif 
