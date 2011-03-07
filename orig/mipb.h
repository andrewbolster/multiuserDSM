#ifndef MIPB_H
#define MIPB_H

#include <cfloat>
#include <list>
#include <iostream>
#include "psd_vector.h"
#include "multiuser_load.h"

#include <boost/threadpool.hpp>

using namespace boost::threadpool;
using namespace std;

enum {
	GRAD,
	ADAPTIVE_GRAD,
	BISECTION,
	OLD
};

struct cost_entry{
	int user;
	int tone;
	double cost;
	cost_entry()
	{
		user=-1;
		cost=DBL_MAX;
		prev=NULL;
		next=NULL;
	}
	friend bool operator<(const cost_entry &left,const cost_entry &right)
	{
		return left.cost < right.cost;
	}
	friend bool operator>(const cost_entry &left,const cost_entry &right)
	{
		return left.cost > right.cost;
	}
	friend ostream& operator <<(ostream& os,const cost_entry& c)
	{
		os << "line " << c.user << "\n";		
		os << "tone " << c.tone << "\n";		
		os << "cost " << c.cost << "\n";		
		return os;
	}

	struct cost_entry *prev;
	struct cost_entry *next;
	
};




class mipb {

	public:
	mipb(const int);
	~mipb();

	int _rate_search_type;
	
	int run();
	double *_p_budget;

	bool _spectral_mask_on;
	double _spectral_mask;
	
	bool _graph_loading;
	bool _simple_search;
	bool _search1;
	bool _search2;

	double *_wxt;
	double *_w;	// per user static weighting

	int *_rate_targets;
	int _rate_tol;

	double _sw;	
	double _min_sw;	
	double _sw1;	

	bool _old;

	int *_bit_mask;

	int _threads;
	
	int calc_delta_p(int,int,int);
	
	bool _greedy;
	bool _show_solution;
	bool _is_frac;
	double _bit_inc;
	
	int _num_threads;
	//boost::shared_ptr<fifo_pool> _tp;
	//fifo_pool *_tp;
	
	void greedy_load(void);

	double _zeta;

	private:

	pool *_tp;

	double *_wp;	// per user dynamic power update weight
	double **_cost;
	double **_no_fext_bits;
	double **_ref_bits;
	double ***_delta_p;
	double ***_mask;
	double **_b;
	int **_init_b;
	double **_p;
	double *_p_used;
	double *_p_lead;
	double _p_ave;
	double *_b_total;
	int *_natural_rate;
	int **_F;
	int _total_bits;
	bool _all_tones_full;

	bool *_behind;

	double *_best_w;

	int **_b_last;
	double **_p_last;

	double *_init_p_budget;

	double *_ref_psd;
	int *_active_channels;
	
	double *_w_max;
	double *_w_min;
	
	double _last_bit_cost;

	psd_vector **_psd;

	int _num_greedys;

	struct cost_entry **_sorted_costs;	
	std::list<struct cost_entry> _L;

	struct cost_entry* _cost_list_head;
	struct cost_entry* _cost_list_tail;

	void load(void);

	void init_cost_matrix(int *,int *);
	void calc_mask(int,int);
	void init_no_fext_bits(void);
	void init_ref_bits(void);
	void init_power(void);
	void reset_data(void);
	bool total_power_constraint_broken(int,int);
	bool spectral_mask_constraint_broken(int,int);
	int min_cost(int*,int*);
	int update_power(int,int);
	int update_wp();
	double cf(int,int);
	double cf2(int,int);
	double var_cf(int,int);
	void init_lines(void);
	void write_current_stats(int,int);
	void write_stats_totals();
	
	void bisect_w(int,int);
	void bisect_w(void);

	void old_rate_algo(void);
	void update_w();
	void update_sw();
	double calc_r_distance();
	void update_w_single(int);

	void bisect_p_budget(int);

	void calc_lost_bits(int,int);

	void dump_matrix(double **,const char *);
	void dump_int_matrix(int **,const char *);

	void copy_b_and_p_to_last(void);

	double rate_spring(int);

	bool rates_converged();

	void recalc_and_sort_costs(int);
	void recalc_costs(int);
	void find_min_sorted_per_tone_costs(int&,int&);
	void get_min_in_list(int&,int&);
	void insert_new_min_cost(int);
	void recalc_all_costs_and_find_min(int&,int&,pool&);

	void threaded_recalc_and_find_min(struct cost_entry* cost,int,int);

	void add_to_cost_list(struct cost_entry *c);
	void add_to_front(struct cost_entry *c);
	void print_cost_list(void);
	void remove_cost_list_head(void);
	void add_to_cost_list_after(struct cost_entry *c,struct cost_entry *current);
	void add_to_cost_list_before(struct cost_entry *c,struct cost_entry *current);
	void cost_list_del(void);

	void dispatch_delta_p_threads(int);

	friend class greedy_nsection;
};


#endif
