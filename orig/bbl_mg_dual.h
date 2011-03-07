#ifndef BBL_H
#define BBL_H

#include "multiuser_greedy.h"

struct bbl_entry {
        int line_id;
	int tone;
        double metric;

	bbl_entry()
	{
		line_id=-1;
		tone=-1;
		metric=0;
	}
	friend bool operator<(const bbl_entry &left,const bbl_entry &right)
	{
		return left.metric < right.metric;
	}

	friend bool operator>(const bbl_entry &left,const bbl_entry &right)
	{
		return left.metric > right.metric;
	}
};


class bbl_mg_dual
{
	
	public:
	bbl_mg_dual();
	int run();
	int *rate_b_target;
	static const int e=5;
	
	private:
	static const double alpha=1.0;
	static const double beta=1.0;

	struct bbl_entry **bbl;
	struct bbl_entry **bbl_a;
	struct bbl_entry **bbl_b;

	int *current_rate_b;
	
	struct bbl_entry bbl_metric(int,int);
	void init_bad_boys_league(void);
	void print_bbl(struct bbl_entry **);
	void migrate_services();
	void set_services();
	void set_loading_algo();
	bool rate_b_targets_good();
	struct bbl_entry *find_bbl_b_min(int);
	struct bbl_entry *find_bbl_a_max(int);
	void move_to_bbl_a(struct bbl_entry *);
	void move_to_bbl_b(struct bbl_entry *);
	void null_entry(struct bbl_entry *);

	//void greedy_multiuser_load(void);
	bbl_multiuser_greedy* mg;

};

#endif
