#ifndef GREEDY_NSECTION
#define GREEDY_NSECTION

#include "multiuser_load.h"
#include "mipb.h"

#include <boost/threadpool.hpp>
#include <boost/bind.hpp>
#include <boost/mem_fn.hpp>
#include <boost/bind/mem_fn.hpp>

using namespace boost::threadpool;


class greedy_nsection {

	public:
	greedy_nsection(const int);
	
	void nsection();
	
	int _num_threads;
	int *_rate_targets;

	int _rate_tol;

	bool _nsection;
	bool _grad_search;

	double *_p_budget;

	private:

	pool *_tp;

	double *_w;
	mipb **_greedy_workers;

	bool *_first_w_run;

};




#endif
