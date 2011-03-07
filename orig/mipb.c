#include "mipb.h"
#include <algorithm>
#include <list>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cfloat>
#include <cassert>
#include "multiuser_load.h"

#include <boost/threadpool.hpp>
#include <boost/bind.hpp>
#include <boost/mem_fn.hpp>
#include <boost/bind/mem_fn.hpp>

using namespace boost::threadpool;

struct mipb_delta_p_thread_args {
	mipb *m;
	int tone;
	int lower;
	int upper;
	int thread_id;
};

void *delta_p_thread(void *p);



mipb::mipb(const int t)
{

	_num_threads=t;

	_delta_p = new double**[DMTCHANNELS];
	_mask = new double**[DMTCHANNELS];
	_cost = new double*[DMTCHANNELS];
	_sorted_costs = new struct cost_entry*[DMTCHANNELS];
	_no_fext_bits = new double*[DMTCHANNELS];
	_ref_bits = new double*[DMTCHANNELS];
	_b = new double*[DMTCHANNELS];
	_b_last = new int*[DMTCHANNELS];
	_init_b = new int*[DMTCHANNELS];
	_p = new double*[DMTCHANNELS];
	_p_last = new double*[DMTCHANNELS];
	_F = new int*[DMTCHANNELS];

	_wp = new double[lines];
	_w = new double[lines];
	_best_w = new double[lines];
	_w_max = new double[lines];
	_w_min = new double[lines];
	_wxt = new double[lines];
	_p_used=new double[lines];
	_p_lead=new double[lines];
	_p_budget=new double[lines];
	_b_total=new double[lines];
	_natural_rate=new int[lines];
	_rate_targets=new int[lines];
	_ref_psd = new double[lines];
	_init_p_budget = new double[lines];
	_active_channels = new int[lines];

	_behind = new bool[lines];

	_bit_mask = new int[lines];

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		_delta_p[tone] = new double*[lines];
		_mask[tone] = new double*[lines];
		_cost[tone] = new double[lines];
		_sorted_costs[tone] = new struct cost_entry[lines];
		_no_fext_bits[tone] = new double[lines];
		_ref_bits[tone] = new double[lines];
		_b[tone] = new double[lines];
		_b_last[tone] = new int[lines];
		_init_b[tone] = new int[lines];
		_p[tone] = new double[lines];
		_p_last[tone] = new double[lines];
		_F[tone] = new int[lines];
	}

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			_delta_p[tone][user] = new double[lines];
			_mask[tone][user] = new double[lines];
			for (int user1=0;user1<lines;user1++) {
				_delta_p[tone][user][user1]=0;
				_mask[tone][user][user1]=0;
			}
			_cost[tone][user]=0;
			_F[tone][user]=0;
			_b[tone][user]=0;
			_p[tone][user]=0;
		}
	}


	for (int user=0;user<lines;user++) {

		_wp[user]=1;
		_w[user]=1;
		_w_min[user]=0;
		_w_max[user]=100;
		_wxt[user]=1;
		_p_budget[user]=0.110;
		_behind[user]=false;
		_rate_targets[user]=NOT_SET;
	}
	//_threads=1;	// default

	/*	
	_psd = new psd_vector*[_num_threads];

	for (int i=0;i<_num_threads;i++) {
		_psd[i] = new psd_vector;
	}
	*/
	
	_p_ave=0.0;

	_spectral_mask_on=true;
	_spectral_mask=dbmhz_to_watts(-30);

	_graph_loading=true;

	_rate_tol=10;

	_min_sw=1e-4;
	_sw=_min_sw;
	_sw1=1e-5;

	_bit_mask[0]=0;		
	_bit_mask[1]=0;		
	_bit_mask[2]=0;		
	
	/*_w[0]=1;
	_w[1]=1;
	_w[2]=1e-6;
	*/

	for (int user=0;user<lines;user++) {
		_ref_psd[user] = line_array[user]->waterfill_level;
		_init_p_budget[user] = 0;
		_p_lead[user]=0;
	}
	/*
	_init_p_budget[0]=0;
	_init_p_budget[1]=0;
	_init_p_budget[2]=0.02818;
	*/
	_greedy=true;
	
	/*
	_w[0]=1;
	_w[1]=1;
	_w[2]=1;
	*/
	//_ref_psd[0]=dbmhz_to_watts(-140);
	//_ref_psd[1]=dbmhz_to_watts(-140);
	//_ref_psd[2]=dbmhz_to_watts(-140);
	//
	//
	_simple_search=false;
	_search1=false;
	_search2=true;

	_cost_list_head=NULL;
	_cost_list_tail=NULL;

	_old=false;

	_show_solution=false;
	
	_is_frac=false;

	_bit_inc=1;

	_p_lead[1]=-0.011;

	_num_greedys=0;

	//_num_threads=1;

	//_tp->size_controller().resize(_num_threads);
	//
	_rate_search_type=ADAPTIVE_GRAD;

	_tp = new pool(_num_threads);

	_zeta=1.1;
}

mipb::~mipb()
{
	/*
	for (int i=0;i<_num_threads;i++) {
		delete _psd[i];
	}
	delete _psd;
	*/
}

void mipb::reset_data()
{

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			for (int user1=0;user1<lines;user1++) {
				_delta_p[tone][user][user1]=0;
				_mask[tone][user][user1]=0;
			}
			_cost[tone][user]=0;
			_F[tone][user]=0;
			_b[tone][user]=0;
			_p[tone][user]=0;
		}
	}

	for (int user=0;user<lines;user++) {
		_b_total[user]=0;
		_p_used[user]=0;
		_wp[user]=1;
		_wxt[user]=1;
	}	
	_total_bits=0;
	//p_ave=0;
	//_L.clear();
	
	if (_search2) {
		cost_list_del();
		_cost_list_head=_cost_list_tail=NULL;
	}
}

double mipb::rate_spring(int user)
{
	double x = (double)(_rate_targets[user]-_b_total[user])/(double)(_rate_targets[user]); // percentage distance

	if(_rate_targets[user] == NOT_SET) 
		return 1;

	double spring = exp(-500*x)+1;

	//if (spring > 10)
	//	printf("Rate spring returned %6.4g on line %d\n",spring,user);
	
	return spring;
}


int mipb::run()
{

	int min_tone,min_user;
	int last;
	
	if (_is_frac) {
		_bit_inc=0.1;
	}	

	if (!_greedy) {
		_simple_search=true;
		_search1=false;
		_search2=false;
	}

	bool target=false;
	for (int user=0;user<lines;user++) {
		if (_rate_targets[user] != NOT_SET) {
			target=true;
			_graph_loading=false;
		}
	}

	if (target) {
		double *b_total_last = new double[lines];
		double *w_last = new double[lines];
		bool first_run=true;	
		bool reset_w=false;
		int hold_counter=0;
		int _osc_tol=_rate_tol/2;
		do {
			//bisect_w();
			print_vector(_w,"w before loading");
	
			switch(_rate_search_type) {
			case BISECTION:
				bisect_w();
				break;
			case GRAD:
				load();
				update_w();
				break;
			case OLD:
				load();
				old_rate_algo();
				break;
			case ADAPTIVE_GRAD:
				load();

				if (first_run) {
					first_run=false;
					memcpy(b_total_last,_b_total,sizeof(double)*lines);
					memcpy(w_last,_w,sizeof(double)*lines);
				}
				else {
					int osc=0;
					for (int user=0;user<lines;user++) {
						if (_rate_targets[user] == NOT_SET)
							continue;
						if (((b_total_last[user] - _rate_targets[user]) > _osc_tol) && ((_rate_targets[user] - _b_total[user]) > _osc_tol)) {
							osc++;
						}
						//if (b_total_last[user] < _rate_targets[user] && _b_total[user] > _rate_targets[user]) {
						if (((_b_total[user] - _rate_targets[user]) > _osc_tol) && ((_rate_targets[user] - b_total_last[user]) > _osc_tol)) {
							osc++;
						}
					}
					if (osc >= 2) {
						_sw/=2;	
						printf("Oscillation detected\n");
						printf("Decreasing step size to %lf\n",_sw);
						memcpy(_w,w_last,sizeof(double)*lines);
						memcpy(_b_total,b_total_last,sizeof(double)*lines);
						printf("Changing w back to :\n");
						print_vector(w_last,"w_last");
						print_vector(b_total_last,"b_total_last");	
						hold_counter=0;
						//continue;
					}
					else {
						if (hold_counter > 0) {
							printf("Holding\n");
							hold_counter--;
						}
						else {
							_sw*=2;
							printf("Increasing step size to %lf\n",_sw);
						}
						memcpy(b_total_last,_b_total,sizeof(double)*lines);
						memcpy(w_last,_w,sizeof(double)*lines);
					}
				}
			
				update_w();
				break;
	
			}		

			//load();
			//old_rate_algo();
			
			int diff[lines];
			for (int user=0;user<lines;user++) {
				diff[user] = _b_total[user] - _rate_targets[user];
			}
			print_vector(_b_total,"b_total");
			print_vector(_w,"w");
			print_vector(diff,"rate delta");
			//getchar();	
			//update_sw();

			//update_w();	

				/*_p_budget[user]=0.02818;
				double max = _p_budget[user];
				double min = 0;
				_init_p_budget[user]=min;	// initial p_budget for preloading is set to zero
				double last;
				while(1) {
					reset_data();
					load();
					//copy_b_and_p_to_last();
					//memcpy(_b_last,_b,sizeof(int)*lines*DMTCHANNELS);
					//memcpy(_p_last,_p,sizeof(double)*lines*DMTCHANNELS);
					//getchar();
					if (_b_total[user] < _rate_targets[user]) {
						min = _init_p_budget[user];
					}
					else if (_b_total[user] > _rate_targets[user]) {
						max = _init_p_budget[user];
						if (_init_p_budget[user] == 0) {
							printf("Line %d is strong, reaches target with no preloading, bisecting\n",user);
							bisect_p_budget(user);
							copy_b_and_p_to_last();
							break;
						}
					}
					if (abs(_b_total[user]-_rate_targets[user]) < _rate_tol) {
						printf("Sweet\n");
						copy_b_and_p_to_last();
						break;
					}
					
					_init_p_budget[user]=(max+min)/2;
					
				}*/

				
		}while (!rates_converged());

	}
	else
		load();


	if (_show_solution) {		
		dump_matrix(_b,"b");
		system("../scripts/graph-dump-matrix.sh /tmp/mipb-b.txt");
		
		//dump_matrix(_p,"p");
		//system("../scripts/graph-dump-matrix.sh /tmp/mipb-p.txt");
		getchar();
		system("killall gnuplot_x11");
	}	

	//dump_matrix(_ref_bits,"ref_bits");
	//dump_matrix(_no_fext_bits,"no_fext_bits");
	//dump_int_matrix(_init_b,"init_b");
	
	init_lines();
	
	calculate_snr();

	if (_graph_loading) {
		write_stats_totals();
	}

	printf("number of greedy evalulations = %d\n",_num_greedys);


}

void mipb::dispatch_delta_p_threads(int tone)
{
	pthread_t th[_threads];
	struct mipb_delta_p_thread_args args[_threads];

	for (int t=0;t<_threads;t++) {
		args[t].m = this;
		args[t].tone = tone;
		args[t].lower = t * (lines/_threads);
		args[t].upper = (t+1)* (lines/_threads);
		if (t==_threads-1) {
			if (args[t].upper != (lines-1)) {
				args[t].upper=lines-1;
			}
		}
		args[t].thread_id=t;
		/*cout << "Thread id" << t << endl;
		cout << "Lower " << args[t].lower << endl;
		cout << "Upper " << args[t].upper << endl;
		getchar();
		*/
		pthread_create(&(th[t]),NULL,delta_p_thread,(void *)&(args[t]));
	}

	for (int t=0;t<_threads;t++) {
		pthread_join(th[t],NULL);
	}
}

void *delta_p_thread(void *p)
{
	struct mipb_delta_p_thread_args *th = static_cast<mipb_delta_p_thread_args *>(p);
	mipb *m = static_cast<mipb *>(th->m);

	//printf("Hello from thread %d!\n",th->thread_id);

	for (int user=th->lower;user<=th->upper;user++) {
		m->calc_delta_p(th->tone,user,th->thread_id);	
	}
	return (void *)NULL;
}


void mipb::update_sw(void)
{
	static double prev_r_distance=DBL_MAX;
	static double min_r_distance=DBL_MAX;
	double r_distance;

	r_distance = calc_r_distance();
	
	if (r_distance <= prev_r_distance) {
		if (r_distance <= min_r_distance) {
			min_r_distance = r_distance;
			for (int user=0;user<lines;user++) {
				_best_w[user] = _w[user];
			}
		}
		_sw*=2;
		printf("Step size increased to %lf\n",_sw);
	}
	else {
		printf("Starting a new traj\n");
		for (int user=0;user<lines;user++) {
                	_w[user] = _best_w[user];
                }

		_sw=_min_sw;
		prev_r_distance=DBL_MAX;
	}

	prev_r_distance = r_distance;
	return;

}

double mipb::calc_r_distance(void)
{
	double sum=0.0;
		
	for (int user=0;user<lines;user++) {	
		if (_rate_targets[user] == NOT_SET)
			continue;
		sum+= pow((_rate_targets[user]-_b_total[user]),2);
	}
	
	return sqrt(sum);
	
	/*
	for (int user=0;user<lines;user++) {
		sum+=abs(_rate_targets[user]-_b_total[user]);
	}

	return sum;
	*/
	
}

void mipb::old_rate_algo(void)
{
	for (int user=0;user<lines;user++) {
                if (_rate_targets[user] == NOT_SET)
                                continue;
		//if (abs(_b_total[user]-_rate_targets[user]) < _rate_tol) {
			//printf("Not changing w\n");
		//	continue;
		//}
		if (_b_total[user] < _rate_targets[user]) {
			//printf("Decreasing w\n");
			_w[user]/=_zeta;
		}
		else if (_b_total[user] > _rate_targets[user]) {
			//printf("Increasing w\n");
			_w[user]*=_zeta;	
		}
		
	}


}


void mipb::bisect_w(void)
{
	for (int user=0;user<lines;user++) {
		if (_rate_targets[user] == NOT_SET)
				continue;
			printf("Searching on line %d\n",user);
					
			double max;
			double min;
			bool converged_early=false;


			load();
			print_vector(_b_total,"b_total");
			if (_b_total[user] < _rate_targets[user]) {		// found max, now find min
				printf("Initial value of w[%d] = %lf was max, now finding min\n",user,_w[user]);
				max = _w[user];
				_w[user]/=2;
				while(1) {
					load();
					print_vector(_b_total,"b_total");
					if (_b_total[user] < _rate_targets[user]) {
						printf("new value of w[%d] = %lf was new max\n",user,_w[user]);
						max = _w[user];
						_w[user]/=2;
					}
					else if (_b_total[user] > _rate_targets[user]) {
						printf("min now found on w[%d] = %lf\n",user,_w[user]);
						min = _w[user];
						break;	
					}
					if (abs(_b_total[user]-_rate_targets[user]) < _rate_tol) {
						printf("converged early\n");
						converged_early=true;
						break;
					}
				}
				
			}
			else if (_b_total[user] > _rate_targets[user]) { 	// found min, now find max
				printf("Initial value of w[%d] = %lf was min, now finding max\n",user,_w[user]);
				min = _w[user];	
				_w[user]*=2;
				while(1) {
					load();
					print_vector(_b_total,"b_total");
					if (_b_total[user] < _rate_targets[user]) {
						printf("max now found on w[%d] = %lf\n",user,_w[user]);
						max = _w[user];
						break;
					}
					else if (_b_total[user] > _rate_targets[user]) {
						printf("new value of w[%d] = %lf was new min\n",user,_w[user]);
						min = _w[user];
						_w[user]*=2;
					}
					if (abs(_b_total[user]-_rate_targets[user]) < _rate_tol) {
						//printf("Sweet\n");
						printf("converged early\n");
						converged_early=true;
						break;
					}
				}
			}
			if (abs(_b_total[user]-_rate_targets[user]) < _rate_tol) {
				//printf("Sweet\n");
				printf("converged early\n");
				converged_early=true;
				continue;
			}
		
			if (converged_early)
				continue;
	
			printf("About to start bisection with w[%d]max = %lf w[%d]min = %lf\n",user,max,user,min);

			_w[user] = (max+min)/2;
			while(1) {
				//reset_data();
				print_vector(_w,"w");
				load();
				print_vector(_b_total,"b_total");
				//getchar();
				if (_b_total[user] < _rate_targets[user]) {
					max = _w[user];
				}
				else if (_b_total[user] > _rate_targets[user]) {
					min = _w[user];	
				}
				if (abs(_b_total[user]-_rate_targets[user]) < _rate_tol) {
					//printf("Sweet\n");
					break;
				}
				_w[user] = (max+min)/2;

				if (_w[user]<1e-20) {
					printf("Cannot achieve rate on line %d\n",user);
					exit(0);
				}
				if (fabs(max-min) < 1e-15) {
					printf("Cannot achieve rate accuracy on line %d\n",user);
					exit(0);
				}
			}

			if (rates_converged())
				return;

	}
}

void mipb::bisect_p_budget(int user)
{
	double max=_p_budget[user];
	double min=0;

	while(1) {
		reset_data();
		load();
		if (_b_total[user] < _rate_targets[user]) {
			min = _p_budget[user];
		}
		else if (_b_total[user] > _rate_targets[user]) {
			max = _p_budget[user];
		}
		if (abs(_b_total[user]-_rate_targets[user]) < _rate_tol) {
			printf("Jeanst\n");
			break;
		}
		
		_p_budget[user]=(max+min)/2;
	}
}

void mipb::copy_b_and_p_to_last(void)
{
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			_b_last[tone][user] = _b[tone][user];
			_p_last[tone][user] = _p[tone][user];
		}
	}

}

void mipb::bisect_w(int user,int rate_user)
{
	_w_min[user]=0;
	_w_max[user]=1e-4;
	//_b_total[user]=0;
	int last=-1;
	while(1) {
		reset_data();
		//init_cost_matrix();
		load();				

		print_vector(_w,"_w");
		print_vector(_b_total,"_b_total");
		//getchar();
		if (user == rate_user) {
			if (_b_total[user] > _rate_targets[user]) {
				//printf("Setting _w_max[%d] = %lf\n",user,_w[user]);
				_w_max[user] = _w[user];
			}
			else if (_b_total[user] < _rate_targets[user]) {
				//printf("Setting _w_min[%d] = %lf\n",user,_w[user]);
				_w_min[user] = _w[user];
			}
			else {
				//printf("Fuck you.\n");
				//printf("_b_total[%d] = %d\n",user,_b_total[user]);
				//printf("_rate_targets[%d] = %d\n",user,_rate_targets[user]);
			}

			if (abs(_b_total[user] - _rate_targets[user]) < _rate_tol) {
				printf("Converged on line %d\n",user);
				break;
			}
		}
		else {
			if (_b_total[rate_user] > _rate_targets[rate_user]) {
				//printf("Setting _w_max[%d] = %lf\n",user,_w[user]);
				_w_min[user] = _w[user];
			}
			else if (_b_total[rate_user] < _rate_targets[rate_user]) {
				//printf("Setting _w_min[%d] = %lf\n",user,_w[user]);
				_w_max[user] = _w[user];
			}
			else {
				//printf("Fuck you.\n");
				//printf("_b_total[%d] = %d\n",user,_b_total[user]);
				//printf("_rate_targets[%d] = %d\n",user,_rate_targets[user]);
			}

			if (abs(_b_total[rate_user] - _rate_targets[rate_user]) < _rate_tol) {
				printf("Converged on line %d\n",user);
				break;
			}
		}
	
		/*if (_w[user] == last) {
			printf("Converged on line %d but rate target not met\n",user);
			break;
			//getchar();
		}*/

		last=_w[user];
		_w[user] = (_w_max[user]+_w_min[user])/2;
		printf("changing weight vector to\n");	
		print_vector(_w,"_w");
		print_vector(_b_total,"_b_total");

	}



}

void mipb::update_w_single(int user)
{

	if (_rate_targets[user] == NOT_SET)
		return;

	if (_w[user] == 0 && _b_total[user] > _rate_targets[user]) {
		printf("Need to give bits to other lines\n");
		for (int other_user=0;other_user<lines;other_user++) {
			if (_rate_targets[other_user] == NOT_SET || (_b_total[other_user] < _rate_targets[other_user])) {
				_w[other_user]+=_sw;		
			}
		}
	}
	else if (_b_total[user] < _rate_targets[user]) {
		_w[user] = _w[user] + _sw*(_rate_targets[user]-_b_total[user]);
	}
	else if( _b_total[user] > _rate_targets[user]) {
        	_w[user] = _w[user] - _sw*(_rate_targets[user]-_b_total[user]);
                 if (_w[user] < 0)
                 	_w[user] = 0;
	}
}


void mipb::update_w()
{
	for (int user=0;user<lines;user++) {
		if (_rate_targets[user] == NOT_SET)
			continue;
		if (abs(_b_total[user] - _rate_targets[user]) < _rate_tol) {
			_w[user]=_w[user];
			continue;
		}
		printf("Before w[%d] = %lf\n",user,_w[user]);
		if (_b_total[user] < _rate_targets[user]) {
			printf("User %d rate is less than target\n",user);
			double update = _sw*(_rate_targets[user]-_b_total[user]);
			if (_w[user] - update < 0) 
				_w[user]*=0.9;
			else 
				_w[user] = _w[user] - update;
		}
		else if (_b_total[user] > _rate_targets[user]) {
			printf("User %d rate is more than target\n",user);
			//_w[user] = _w[user] - _sw*(_rate_targets[user]-_b_total[user]);
			_w[user] = _w[user] + _sw*(_b_total[user]-_rate_targets[user]);
		}
		printf("After w[%d] = %lf\n",user,_w[user]);
	}
}

void mipb::load()
{
	int min_tone,min_user;
	int min_tone2,min_user2;
	//init_ref_bits();
	reset_data();

	//pool tp(_num_threads);
	_num_greedys++;	
	
	//init_no_fext_bits();
	/*
	for (int tone=0;tone<DMTCHANNELS;tone++) {	// preeoading based on ref bits done here
		for (int user=0;user<lines;user++) {
			//_b_total[user]+=_init_b[tone][user]=_b[tone][user] = MIN(rint((_ref_bits[tone][user])),MAXBITSPERTONE);
			_b_total[user]+=_init_b[tone][user]=_b[tone][user] = MIN((int)(_ref_bits[tone][user]),MAXBITSPERTONE);
			if (_b[tone][user] == MAXBITSPERTONE)
				_F[tone][user]=1;
		}
	}
	
	init_power();	
	*/
	
	init_cost_matrix(&min_tone,&min_user);
	_all_tones_full=false;
	while(1) {
	
		if (_greedy) {	
			if (_simple_search)
				min_cost(&min_tone,&min_user);	
			else if(_search1)
				find_min_sorted_per_tone_costs(min_tone,min_user);	
			else if (_search2) 
				get_min_in_list(min_tone,min_user);
		}
			
		//printf("min tone = %d\tmin user =%d\n",min_tone,min_user);
		//printf("min tone2 = %d\tmin user2 =%d\n",min_tone2,min_user2);

		if (_all_tones_full) {
			//printf("All tones full\n");
			if (_graph_loading) {
				write_current_stats(min_tone,min_user);
			}	
			break;
		}
		_b[min_tone][min_user]+=_bit_inc;        // add bit to min cos tone
		_b_total[min_user]+=_bit_inc;
		_total_bits+=_bit_inc;

		/*if (_b[min_tone][min_user] > MAXBITSPERTONE) {
			printf("b = %6.4lf\n",_b[min_tone][min_user]);
		}*/
		//assert(_b[min_tone][min_user] <= MAXBITSPERTONE);

		update_power(min_tone,min_user);

		if (_b[min_tone][min_user] >= MAXBITSPERTONE) {
			_F[min_tone][min_user] = 1;
		}

		if (!_greedy)
			update_wp();

		if (_total_bits % 50 == 0) {
			if (_graph_loading) {
				write_current_stats(min_tone,min_user);
			}
		}
		

		for (int user=0;user<lines;user++) {
			//tp.schedule(boost::bind<int>(boost::mem_fn(&mipb::calc_delta_p),_3,min_tone,user,0));
			//tp.schedule(boost::bind(&mipb::calc_delta_p,_3,min_tone,user,0));
			_tp->schedule(boost::bind(boost::mem_fn(&mipb::calc_delta_p),this,min_tone,user,0));
			//tp.wait();
			//calc_delta_p(min_tone,user,0);
		}
		_tp->wait();
		//dispatch_delta_p_threads(min_tone);

		
		if (!_greedy) {
			//recalc_all_costs_and_find_min(min_tone,min_user,tp);
		}
		else {
			if (_simple_search)
				recalc_costs(min_tone);
			else if(_search1)
				recalc_and_sort_costs(min_tone);
			else if(_search2) {
				recalc_and_sort_costs(min_tone);
				insert_new_min_cost(min_tone);
			}
				
		}
		//printf("Min tone2 = %d\tMin user = %d\n",min_tone2,min_user2);

	}
	/*
	for (int victim=0;victim<lines;victim++) {
		for (int xtalker=0;xtalker<lines;xtalker++) {
			if (victim == xtalker) 
				continue;
			calc_lost_bits(victim,xtalker);
		}
	}
	*/

	if (_show_solution) {		
		dump_matrix(_b,"b");
		system("../scripts/graph-dump-matrix.sh /tmp/mipb-b.txt");
		
		//dump_matrix(_p,"p");
		//system("../scripts/graph-dump-matrix.sh /tmp/mipb-p.txt");
		getchar();
		system("killall gnuplot_x11");
	}
	
	//print_vector(_b_total,"b_total");
	//print_vector(_init_p_budget,"init_p_budget");

	/*	
	dump_matrix(_ref_bits,"ref_bits");
	dump_int_matrix(_init_b,"init_b");
	dump_int_matrix(_b,"b");

	system("../scripts/graph-dump-matrix.sh /tmp/mipb-ref_bits.txt");
	system("../scripts/graph-dump-matrix.sh /tmp/mipb-init_b.txt");
	system("../scripts/graph-dump-matrix.sh /tmp/mipb-b.txt");
	system("../scripts/graph-dump-vector.sh /tmp/vector-ref_noise_line_0.txt");
	system("../scripts/graph-dump-vector.sh /tmp/vector-ref_noise_line_1.txt");
	system("../scripts/graph-dump-vector.sh /tmp/vector-ref_noise_line_2.txt");

	getchar();

	system("killall gnuplot_x11");
	*/
}


void mipb::greedy_load()
{
	int min_tone,min_user;
	int min_tone2,min_user2;
	reset_data();

	assert(_greedy==true);
	
	_all_tones_full=false;
	init_cost_matrix(&min_tone,&min_user);
	while(1) {
	
		if (_simple_search)
			min_cost(&min_tone,&min_user);	
		else if(_search1)
			find_min_sorted_per_tone_costs(min_tone,min_user);	
		else if (_search2) 
			get_min_in_list(min_tone,min_user);
			

		//printf("min tone = %d min_user = %d\n",min_tone,min_user);

		if (_all_tones_full) {
			//printf("All tones full\n");
			if (_graph_loading) {
				write_current_stats(min_tone,min_user);
			}	
			break;
		}
		_b[min_tone][min_user]+=_bit_inc;        // add bit to min cos tone
		_b_total[min_user]+=_bit_inc;
		_total_bits+=_bit_inc;

		update_power(min_tone,min_user);

		if (_b[min_tone][min_user] >= MAXBITSPERTONE) {
			_F[min_tone][min_user] = 1;
		}

		for (int user=0;user<lines;user++) {
			calc_delta_p(min_tone,user,0);
		}

		if (_simple_search) {
			recalc_costs(min_tone);
		}
		else if(_search1) {
			recalc_and_sort_costs(min_tone);
		}
		else if(_search2) {
			recalc_and_sort_costs(min_tone);
			insert_new_min_cost(min_tone);
		}
				

	}

	if (_show_solution) {		
		dump_matrix(_b,"b");
		system("../scripts/graph-dump-matrix.sh /tmp/mipb-b.txt");
		
		//dump_matrix(_p,"p");
		//system("../scripts/graph-dump-matrix.sh /tmp/mipb-p.txt");
		getchar();
		system("killall gnuplot_x11");
	}
	
	//print_vector(_b_total,"b_total");
	//print_vector(_w,"w");
	//print_vector(_init_p_budget,"init_p_budget");

	/*	
	dump_matrix(_ref_bits,"ref_bits");
	dump_int_matrix(_init_b,"init_b");
	dump_int_matrix(_b,"b");

	system("../scripts/graph-dump-matrix.sh /tmp/mipb-ref_bits.txt");
	system("../scripts/graph-dump-matrix.sh /tmp/mipb-init_b.txt");
	system("../scripts/graph-dump-matrix.sh /tmp/mipb-b.txt");
	system("../scripts/graph-dump-vector.sh /tmp/vector-ref_noise_line_0.txt");
	system("../scripts/graph-dump-vector.sh /tmp/vector-ref_noise_line_1.txt");
	system("../scripts/graph-dump-vector.sh /tmp/vector-ref_noise_line_2.txt");

	getchar();

	system("killall gnuplot_x11");
	*/
}





void mipb::recalc_all_costs_and_find_min(int &min_tone,int &min_user,pool &tp)
{

		
	double min=DBL_MAX;
	_all_tones_full=true;

	/*	
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			_cost[tone][user] = cf(tone,user);
			if (_cost[tone][user] < min) {
				//printf("New min found at tone %d user %d = %6.4g\n",tone,user,_cost[tone][user]);
				min=_cost[tone][user];
				min_tone=tone;
				min_user=user;
				_all_tones_full=false;
			}
			//if (_F[tone][user] == 0) {
			//	user_full[user]=false;
			//}
		}
	}
	
	printf("min found at tone %d user %d = %6.4g\n",min_tone,min_user,_cost[min_tone][min_user]);
	*/
	
	//_tp->resize(_num_threads);
	//_tp = fifo_pool::create_pool(_num_threads);
	struct cost_entry per_thread_costs[_num_threads];

	for (int t=0;t<_num_threads;t++) {
		int lower = t * (DMTCHANNELS/_num_threads);
		int upper = (t+1) * (DMTCHANNELS/_num_threads);
		tp.schedule(boost::bind(boost::mem_fn(&mipb::threaded_recalc_and_find_min),this,&per_thread_costs[t],lower,upper));
	}

	tp.wait();

	min=DBL_MAX;

	for (int t=0;t<_num_threads;t++) {
		//cout << per_thread_costs[t];
		if (per_thread_costs[t].cost < min) {
			min = per_thread_costs[t].cost;
			min_tone = per_thread_costs[t].tone;
			min_user = per_thread_costs[t].user;
			_all_tones_full=false;
		}
	}

	
	//printf("min found at tone %d user %d = %6.4g\n",min_tone,min_user,_cost[min_tone][min_user]);
	
	//getchar();
}


void mipb::threaded_recalc_and_find_min(struct cost_entry* cost, int lower_tone,int upper_tone)
{
	//printf("Scheduled lower = %d upper = %d\n",lower_tone,upper_tone);
	double min=DBL_MAX;
	
	for (int tone=lower_tone;tone<upper_tone;tone++) {
		for (int user=0;user<lines;user++) {
			_cost[tone][user] = cf(tone,user);
			if (_cost[tone][user] < min) {
				//printf("New min found at tone %d user %d = %6.4g\n",tone,user,_cost[tone][user]);
				min=_cost[tone][user];
				cost->cost=min;
				cost->tone=tone;
				cost->user=user;
			}
		}
	}

	//cout << "min cost found \n" << *cost;

}


void mipb::get_min_in_list(int& min_tone,int& min_user)
{
	
	while (1) {
		/*
		std::list<struct cost_entry>::iterator i;
		i=_L.begin();

		if ((*i).cost ==DBL_MAX) {
			_all_tones_full=true;
			return;
		}
	
		min_tone=(*i).tone;
		min_user=(*i).user;
		*/

		if (_cost_list_head->cost == DBL_MAX) {
			_all_tones_full=true;
			return;
		}

		min_tone=_cost_list_head->tone;
		min_user=_cost_list_head->user;

		if (total_power_constraint_broken(min_tone,min_user) || spectral_mask_constraint_broken(min_tone,min_user)) {	// this cost is no longer valid
			//cout << "entry = " << *_cost_list_head;
			//cout << "Does not work" << endl;
			//_L.erase(i);	// erase this value as it is not valid	
			remove_cost_list_head();
			recalc_and_sort_costs(min_tone);	// recalc for this tone
			insert_new_min_cost(min_tone);		// reinsert in list
			continue;				// start again
		}
		else {
			//cout << "min_entry =" << *_cost_list_head << endl;
			//_L.erase(i);
			remove_cost_list_head();
			break;
		}

	}
}

void mipb::insert_new_min_cost(int tone)
{
	struct cost_entry c;
	c = _sorted_costs[tone][0];

	//cout << "cost to insert " << c << endl;
	/*
	//std::list<struct cost_entry>::iterator i;
	std::list<struct cost_entry>::reverse_iterator r;
	int j=_L.size();
	//for(i=_L.end();i>_L.begin();i++,j--) {
	for(r=_L.rbegin();r!=_L.rend();r++,j--) {
		//cout << "list position " << j << endl;	
		//cout << *r;
		if (c.cost > (*r).cost) {
			//cout << "Inserted at position " << j;
			_L.insert(r.base(),c);
			return;
		}	
	}	
	*/



	if (c.cost == DBL_MAX) {
		//cout << "Adding at tail" << endl;
		add_to_cost_list_after(&c,_cost_list_tail);
		return;
	}

		
	struct cost_entry *current=_cost_list_tail;
	int i=DMTCHANNELS-1;
	while (current!=NULL) {
		if (c.cost > current->cost) {
			//if (i==0) {
			//	cout << c << endl;
			//	cout << *current << endl;
			//}
			//cout << "Adding at position " << i << endl;
			//cout << c << endl;
			add_to_cost_list_after(&c,current);
			//print_cost_list();
			//getchar();
			return;
		}
		current=current->prev;
		i--;
	}
	
	//cout << "Adding to front" << endl;	
	add_to_front(&c);
	return;

	/*	
	struct cost_entry *current=_cost_list_head;
	int i=0;
	while (current!=NULL) {
		if (c.cost < current->cost) {
			//if (i==0) {
			//	cout << c << endl;
			//	cout << *current << endl;
			//}
			//cout << "Adding at position " << i << endl;
			//cout << c << endl;
			//add_to_cost_list_after(&c,current);
			add_to_cost_list_before(&c,current);
			//print_cost_list();
			//getchar();
			return;
		}
		current=current->next;
		i++;
	}

	//cout << "Adding to end" << endl;
	add_to_cost_list_after(&c,_cost_list_tail);
	return;
	*/

}

void mipb::recalc_and_sort_costs(int tone)
{

	double min_cost=DBL_MAX;
	int min_user;
	for (int user=0;user<lines;user++) {
		_sorted_costs[tone][user].cost=_cost[tone][user] = cf(tone,user);
                _sorted_costs[tone][user].user=user;
                _sorted_costs[tone][user].tone=tone;
		if (_cost[tone][user] < min_cost) {
			min_cost=_cost[tone][user];
			min_user=user;
		}
	}
	_sorted_costs[tone][0].cost=min_cost;
	_sorted_costs[tone][0].user=min_user;
	//std::sort(_sorted_costs[tone],_sorted_costs[tone]+lines);	// sort so minimum cost is in element zero

}

void mipb::recalc_costs(int tone)
{
	for (int user=0;user<lines;user++) {
                _cost[tone][user] = cf(tone,user);
        }
}

void mipb::find_min_sorted_per_tone_costs(int& min_tone,int& min_user)
{ 

	_all_tones_full=false;	
	double min=DBL_MAX;
	for (int tone=0;tone<DMTCHANNELS;tone++) {	// search all tones, min cost is the first element
		if (_sorted_costs[tone][0].cost < min) {
			int user=_sorted_costs[tone][0].user;
			if (total_power_constraint_broken(tone,user) || spectral_mask_constraint_broken(tone,user)) {	// this cost is no longer valid
				recalc_and_sort_costs(tone);
				_F[tone][user] = 1;
				tone--;		// try again on this tone after recalculating
				continue;
			}	
			min_user=user;
			min_tone=tone;
			min=_sorted_costs[tone][0].cost;
		}
	}
		
	if (min == DBL_MAX)
		_all_tones_full=true;

}


void mipb::calc_lost_bits(int victim, int xtalker)
{
	double p[lines];
	double b[lines];
	double total_lost_bits=0.0;	

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		memcpy(p,_p[tone],sizeof(double)*lines);
		//print_vector(p,"p_before");
		//print_vector(_b[tone],"b_beore");
		p[xtalker]=0;
		calculate_b_vector_from_psd(p,12.95,tone,b);
		//print_vector(p,"p");
		//print_vector(b,"b");
		total_lost_bits+=b[victim]-_b[tone][victim];
	}

	printf("Total bits lost on line %d due to line %d = %lf\n",victim,xtalker,total_lost_bits);
}

void mipb::init_cost_matrix(int *min_tone,int *min_user)
{

	cost_entry temp_c[DMTCHANNELS];

	double min=DBL_MAX;
        for (int tone=0;tone<DMTCHANNELS;tone++) {
                for(int user=0;user<lines;user++) {
                        calc_delta_p(tone,user,0);
                        _sorted_costs[tone][user].cost=_cost[tone][user] = cf(tone,user);
			//printf("cf returned %g\n",_cost[tone][user]);
			_sorted_costs[tone][user].user=user;
			if (_cost[tone][user] < min) {
				min=_cost[tone][user];
				*min_tone=tone;
				*min_user=user;
			}
			//printf("Cost on user %d tone %d = %g\n",user,tone,_cost[tone][user]);
			//calc_mask(tone,user);
                }
		//std::sort(_sorted_costs[tone],_sorted_costs[tone]+lines);
		recalc_and_sort_costs(tone);
		if (_search2) {
			temp_c[tone] = _sorted_costs[tone][0];
			//_L.push_back(_sorted_costs[tone][0]);

		}
        }

	if(_search2) {
		//_L.sort();
		std::sort(temp_c,temp_c+DMTCHANNELS);
		//cout << _L.size() << endl;
		//getchar();
	}
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		add_to_cost_list(&temp_c[tone]);
	}
	//print_cost_list();
	//exit(0);
}

void mipb::calc_mask(int tone,int user)
{

	double ref_psd = minus_36_5_dbmhz_watts;

	for (int victim=0;victim<lines;victim++) {
		if (victim == user)
			;
		else {
			double b[lines];
			double p[lines];
			memset(p,0,sizeof(double)*lines);
			memset(b,0,sizeof(double)*lines);
			double bits = no_fext_bits(tone,victim,ref_psd);
			if (bits < 1) {
				_mask[tone][user][victim]=0;
				continue;
			}
			p[user]=ref_psd;
			p[victim]=ref_psd;
			calculate_b_vector_from_psd(p,12.95,tone,b);	// FIXME
			double lost_bits = bits - b[victim];
			_mask[tone][user][victim]=_w[victim]*lost_bits;
			if (_mask[tone][user][victim] < 0) {
				printf("Whoops\n");
				printf("Victim = %d\n",victim);
				printf("User = %d\n",user);
				printf("Bits = %lf\n",bits);
				print_vector(p,"p");
				print_vector(b,"b");
				getchar();
			}
		}
	}

}


bool mipb::rates_converged()
{
	for (int user=0;user<lines;user++) {
		if (_rate_targets[user] == NOT_SET)
			continue;
		if (abs(_rate_targets[user] - _b_total[user]) > _rate_tol) {
			return false;
		}
	}

	return true;
}


int mipb::min_cost(int *channel,int *line)
{

	int tone,user;
	int min_tone,min_user;
	double min=DBL_MAX;
	bool *user_full = new bool[lines];	

	for (int user=0;user<lines;user++)
		user_full[user] = true;

	_all_tones_full=true;		// if i put this variable below this line, it segfaults.... compiler bug or weird bug in my code??

	for (tone=0;tone<DMTCHANNELS;tone++) {
		for (user=0;user<lines;user++) {
			if (_F[tone][user]==1) {
				continue;
			}
			user_full[user]=false;		// if we have reached here, the current current is not full!
			assert(_cost[tone][user] > 0);
			//assert(_cost[tone][user] < DBL_MAX);
			if (_cost[tone][user] < min) {
				if (total_power_constraint_broken(tone,user) || spectral_mask_constraint_broken(tone,user)) {
					_F[tone][user] = 1;
					continue;
				}
				_all_tones_full=false;
				min=_cost[tone][user];
				min_tone=tone;
				min_user=user;
					
			}
		}
	}
	
	*channel=min_tone;
	*line=min_user;

	return 0;
}

int mipb::update_wp()
{

	double a=1e3;
	double b=99;
	double c=9190.24;
	

	double p_ave=0.0;
	double div=lines;
	for (int user=0;user<lines;user++) {
		if (!_behind[user])
			p_ave+=_p_used[user];
		else
			div--;
	}
	p_ave/=div;
	
	for (int user=0;user<lines;user++) {
		//double x=(_p_used[user]*100/_p_budget[user])-p_ave;
		//double x = _p_used[user] - (_p_ave+_p_lead[user]);
		//double x = _p_used[user] - _p_ave;
		double x = _p_used[user] - p_ave;
		//_wp[user] = UNdB(2000/(100*(p_ave+0.00001))*(x));
		//_wp[user] = UNdB(10*(_p_used[user]*100/_p_budget[user] - _p_ave));
		//_wp[user] = pow(exp((_p_used[user] - _p_ave)),100);
		
		//linear plus sigmoid
		//_wp[user] = a*x+1.001-(a*x+1)*exp(-c*x)/(exp(-c*x) + b);
	

		// linear	
		/*		
		double grad = 1/(1.5*fabs(x));
		_wp[user] = grad*x+1;
		
		if (_behind[user])
			_wp[user]=1;
		*/
		//_wp[user]=1e2*x+1;
		//_wxt[user] = grad*x/10 + 1;
		//
		if (_behind[user]) {
			_wp[user]=0;
			continue;
		}
		
		if (x < -100e-3) {
			//printf("Looks like used %d is behind\n",user);
			_behind[user]=true;
			/*
			printf("p_ave = %g\n",p_ave);
			for (int user=0;user<lines;user++) {
				printf("p_used[%d] =%g\n",user,_p_used[user]);
			}
			*/
			for (int tone=0;tone<DMTCHANNELS;tone++) {
				_F[tone][user]=1;
			}
			continue;
		}

		/*		
		_wp[user] = 1e3*(x/_last_bit_cost)+1;
	
		if (_wp[user] < 0)	
			_wp[user]=10*DBL_MIN;
		*/
	
			
		if (x > 0) {
			//_behind[user]=0;
			//_wp[user] = (1e5/_last_bit_cost)*x+1;
			assert(_last_bit_cost > 0);
			double z = MIN(x/_last_bit_cost,709);
			_wp[user] = exp(0.25*z);
			
			//double z = x/_last_bit_cost;
			//_wp[user] = 1e6*x+1;
			
			//double z = x;
			//_wp[user] = exp(100000*z);
			if (_wp[user] >= DBL_MAX) {
				printf("x = %6.4g\n",x);
				printf("z = %lf\n",z);
				printf("last_bit_cost = %6.4g\n",_last_bit_cost);
				getchar();
			}
			//_wp[user] = 1e5*x*_w[user]+1;
			assert(_wp[user] > 0);
		}	
		else {
			_wp[user] = 1;
			//_behind[user]++;
		}
		

		/*
		double z = MIN(x/_last_bit_cost,709);
		_wp[user] = exp(0.5*z);

		if (_behind[user])
			_wp[user]=1;
		*/
	
		/*	
		if (_behind[user] >= 100 && _b_total[user] > 0) {
			_behind[user]=0;
			printf("User %d seems to be behind\n",user);
		}
		*/
	
		// original
		//
		//_wp[user] = pow(exp(x),100);

		// step
	
		/*		
		if (x > 0) {
			_wp[user]=1e15;
		}
		else if (x < 0) {
			_wp[user]=1e-15;
		}
		else {
			_wp[user]=1;
		}*/

		//_wp[user] = UNdB(1e5*x);
		//assert(_wp[user] < DBL_MAX);
		assert(_wp[user] > 0);
	}

}

void mipb::init_no_fext_bits()
{
	double ref_psd = minus_36_5_dbmhz_watts;
	
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			if (tone <= _active_channels[user]-1)
				_no_fext_bits[tone][user] = no_fext_bits(tone,user,line_array[user]->waterfill_level);
			else
				_no_fext_bits[tone][user] = 0;
			//printf("No fext bits on line %d tone %d = %lf\n",user,tone,_no_fext_bits[tone][user]);
		}
	}

}

void mipb::init_ref_bits()
{
	double ref_psd = dbmhz_to_watts(-60);

	double power_check[lines];
	memset(power_check,0,sizeof(double)*lines);

	/*
	for (int user=0;user<lines;user++) {
		_active_channels[user] = MIN(DMTCHANNELS,_p_budget[user]/_ref_psd[user]);
		printf("Line %d active channels = %d\n",user,_active_channels[user]);
	}
	*/

	double p_mat[lines][DMTCHANNELS];	// do not be alarmed, this is the other way around than usual on purpose
	for (int user=0;user<lines;user++) {
		double noise[DMTCHANNELS];
		memset(noise,0,sizeof(double)*DMTCHANNELS);
		for (int xtalker=0;xtalker<lines;xtalker++) {
		//	if (xtalker == user || _rate_targets[xtalker] == NOT_SET)
		//		continue;
			for (int tone=0;tone<DMTCHANNELS;tone++) {
				noise[tone]+=_ref_psd[user]*get_xtalk_gain(xtalker,user,tone);
			}
		}
		char tag[100];
		sprintf(tag,"ref_noise_line_%d",user);
		dump_vector(noise,tag);
		int ret=waterfill(user,_init_p_budget[user],NULL,p_mat[user],NULL);	// waterfill on all lines
		printf("Waterfilling rate on line %d = %d\n",user,ret);
	}

	
	for (int tone=0;tone<DMTCHANNELS;tone++) {	
		double p[lines];
		double b[lines];
		for (int user=0;user<lines;user++) {
			/*if (tone <= _active_channels[user]-1)
				p[user] = nearest_integer_bit_p(tone,user,_ref_psd[user]);
			else
				p[user] = 0;
			*/
			p[user] = p_mat[user][tone];

			power_check[user]+=p[user];
		}
		calculate_b_vector_from_psd(p,12.95,tone,b);	
		
		for (int user=0;user<lines;user++) {
			_ref_bits[tone][user] = b[user];
		}
	}

	for (int user=0;user<lines;user++) {
		printf("power_check[%d] = %lf\n",user,power_check[user]);
	}	
}

double mipb::cf(int tone,int line_id)
{

	double cost=0.0;
	double a=0;
	double b=0;
	double c[lines];

	if (total_power_constraint_broken(tone,line_id) || spectral_mask_constraint_broken(tone,line_id)) {
		//printf("Setting tone %d line %d full because of spectral mask or total power broken\n",tone,line_id);
		_F[tone][line_id] = 1;
		//printf("Returning max dbl as tone says power constraint broken\n");	
		return DBL_MAX;
	}

	if (_F[tone][line_id] == 1) {
		//printf("Returning max dbl as tone says it is full\n");	
		return DBL_MAX;
	}

	if (_b[tone][line_id] >= MAXBITSPERTONE) {
		_F[tone][line_id]=1;
		//printf("Returning max dbl as tone says it is past MAX\n");	
		return DBL_MAX;
	}

	for (int user=0;user<lines;user++) {
		if (user == line_id) {
			if (_old) {
				cost+=_delta_p[tone][line_id][user];		// power for the actual bit loaded
			}
			else if (!_greedy) {
				a+=cost+=_wp[user]*_delta_p[tone][line_id][user];
			}
			else {
				c[line_id]=cost+=_w[line_id]*_delta_p[tone][line_id][user];		// power for the actual bit loaded
			}
		}
		else {
			/*if (_b[tone][user] == 0) {	// victim line has no bits yet
				//printf("Victim has no bits yet\n");
				if (_no_fext_bits[tone][user] > 1) { 	// victim can theoretically use this tone
					int b[lines];
					double p[lines];
					double p_old[lines];
					double g[lines];
					//printf("Line looking for cost is %d\n",line_id);
					//printf("Victim is %d\n",user);
					for (int user1=0;user1<lines;user1++) {
						g[user1] = line_array[user1]->gamma[line_array[user1]->service[tone]];
					}
					memset(p,0,sizeof(double)*lines);
					memset(p_old,0,sizeof(double)*lines);
					memcpy(b,&_b[tone][0],sizeof(int)*lines);	// get current b vector
					//b[user]=(int)_ref_bits[tone][user];					// pretend victim has one bit
					//b[user]=MIN((int)_no_fext_bits[tone][user],_w[user]*_ref_bits[tone][user]);
					b[user]=1;
					_psd->calc(b,g,tone,p_old);				// get psd vector
					b[line_id]++;		
					_psd->calc(b,g,tone,p);				// get psd vector

					cost+=_wp[user]*(p[user]-p_old[user]);
		
					if (_wp[user]*(p[user]-p_old[user]) < 0) {
						printf("Line_id = %d\n",line_id);
						printf("victim = %d\n",user);
						printf("tone = %d\n",tone);
						print_vector(_b[tone],"b_current");
						print_vector(p_old,"p_old");
						print_vector(b,"b");
						print_vector(p,"p");
						print_vector(_wp,"wp");
						getchar();
					}
					//printf("Extra cost incurred is %6.4g\n",p[user]-p_old[user]);
					
				}
			}*/
			//cost+=_wxt[user]*_delta_p[tone][line_id][user];		// power for the extra crosstalk incurred
			//cost+=_wp[user]*_delta_p[tone][line_id][user]+_mask[tone][line_id][user];		// power for the extra crosstalk incurred
			if(_old) {
				cost+=_delta_p[tone][line_id][user]; 	// power for the extra crosstalk incurred
			}
			else if (!_greedy) {
				b+=cost+=_wp[user]*_delta_p[tone][line_id][user];
			}
			else {
				if ((1-_w[user]) > 0)
					cost+=c[user]=(1-_w[user])*_delta_p[tone][line_id][user]; 	// power for the extra crosstalk incurred
				//cost+=c[user]=_delta_p[tone][line_id][user]; 	// power for the extra crosstalk incurred
				//cost+=c[user]=(exp(-0.01*_w[user]))*_wp[user]*_delta_p[tone][line_id][user]; 	// power for the extra crosstalk incurred
				//cost+=_w[line_id]*_wp[user]*_delta_p[tone][line_id][user]; 	// power for the extra crosstalk incurred
			}
		}	
	
	}

	/*	
	if (b > a && _wp[line_id] < 0.9) {
		printf("Cost for bit on line %d tone %d is dominated by crosstalk terms even though wp is less than 1\n",line_id,tone);
		printf("_wp[%d] = %g\ta = %g\tb = %g\n",line_id,_wp[line_id],a,b);
		print_vector(_p_used,"_p_used");
		printf("Total bits = %d\n",_total_bits);
		getchar();
	}*/

	//printf("Total cost is %6.4g\n",cost);
	//getchar();
	//

	/*	
	if (b > 2*a) {
		printf("Crosstalk twice price of bit on line %d tone %d\n",line_id,tone);
		print_vector(_wp,"wp");
		print_vector(_b_total,"b_total");
	}
	*/

	if (cost == 0) {
		printf("WTF\n");
		printf("Line_id = %d\n",line_id);
		printf("tone = %d\n",tone);
		//print_vector(_b[tone],"b_current");
		//print_vector(p_old,"p_old");
		//print_vector(b,"b");
		//print_vector(p,"p");
		print_vector(_w,"w");
		getchar();
	}

	/*
	for (int user=0;user<lines;user++) {
		if (c[user] < 0) {
			printf("Line_id = %d\n",line_id);
			print_vector(c,"c");
			print_vector(_w,"w");
			print_vector(_wp,"wp");
			print_vector(_delta_p[tone][line_id],"delta_p");
			print_vector(_b[tone],"b");
			double ratio = c[line_id]/(-1*c[user]);
			printf("Ratio = %6.4g\n",ratio);
			if (ratio < 10)
				getchar();
			break;
		}
	}
	*/

	if (_old)
		return _w[line_id]*cost;
	else 
		return cost;
}

double mipb::var_cf(int tone,int line_id)
{

	double new_p_ave=0.0;
	double new_var=0.0;

	for (int user=0;user<lines;user++ ) {
		new_p_ave+=_p_used[user] + _delta_p[tone][line_id][user];
	}

	new_p_ave/=lines;

	// calculate new _p_ave

	// find min new variance in _p_used
	for (int user=0;user<lines;user++) {
		//new_var+=pow((_p_used[user] + _delta_p[tone][line_id][user]) - new_p_ave),2);
	}

	new_var/=lines;
	
	return new_var;
	
}

bool mipb::total_power_constraint_broken(int tone,int user)
{
        for (int i=0;i<lines;i++) {
                //if (p_used[i] + _delta_p[tone][user][i] > _p_budget[i]) {
                if (_p_used[i] + _delta_p[tone][user][i] > _p_budget[i]) {
                	//printf("adding a bit on line %d tone %d will break the power constraint on line %d\n",user,tone,i);
			//printf("p_budget[%d] = %lf\n",i,_p_budget[i]);
                        return true;
		}
	}

	return false;
}

void mipb::init_power()
{

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		double g[lines];
		int b[lines];
		for (int user=0;user<lines;user++) {
			g[user] = line_array[user]->gamma[line_array[user]->service[tone]];
			b[user] = (int)_b[tone][user];		
		}
		//_psd[0]->calc(&_b[tone][0],g,tone,&_p[tone][0]);
		if (_is_frac)
			calculate_psd_vector(&_b[tone][0],g,tone,&_p[tone][0]);
		else 
			_psd[0]->calc(b,g,tone,&_p[tone][0]);

		for (int user=0;user<lines;user++) {
			assert(_p[tone][user] >= 0);
			_p_used[user]+=_p[tone][user];
		}
	}
	
	double p_sum=0.0;
	for (int user=0;user<lines;user++) {
		p_sum+=_p_used[user];
		printf("p_used[%d] = %6.4g\n",user,_p_used[user]);
		assert(_p_used[user] - _p_budget[user] < 0.001);
	}
	_p_ave=p_sum/lines;

	dump_matrix(_p,"init_power");

	printf("p_ave = %6.4g\n",_p_ave);


}
              
int mipb::update_power(int tone, int user)
{

	_last_bit_cost=0.0;
        for (int i=0;i<lines;i++) {
                _last_bit_cost+=_p[tone][i]+=_delta_p[tone][user][i];
                _p_used[i]+=_delta_p[tone][user][i];
		_p_ave+=_delta_p[tone][user][i]/lines;
        }
	//printf("Last bit cost = %6.4g\n",_last_bit_cost);
	//printf("bit cost = %6.4g\n",_delta_p[tone][user][user]);
	return 0;
}


int mipb::calc_delta_p(int tone,int line_id,int thread_id)
{
	
	double old_p_tot=0.0,new_p_tot=0.0;

	int b[lines];
	double p[lines];
	double old_p[lines];
	double g[lines];


	for (int user=0;user<lines;user++) {
		g[user] = line_array[user]->gamma[line_array[user]->service[tone]];
		//old_p_tot+=_p[tone][user];
		old_p[user]=_p[tone][user];
		b[user] = (int)_b[tone][user];
		if (b[user] > MAXBITSPERTONE) {
			printf("Boom b[%d] = %d\n",user,b[user]);
		}
		assert(b[user] >= 0);
		assert(b[user] <= MAXBITSPERTONE);
	}

	b[line_id]++;

	if (_is_frac) {
		_b[tone][line_id]+=_bit_inc;
		calculate_psd_vector(&_b[tone][0],g,tone,p);
		_b[tone][line_id]-=_bit_inc;
	}
	else {
		//_psd[thread_id]->calc(b,g,tone,p);
		calculate_psd_vector(b,g,tone,p,cache);
	}

	for (int user=0;user<lines;user++) {
		if (p[user] < 0) {
			//printf("No solution on tone %d\n",tone);
			_F[tone][line_id]=1;
			for (int user=0;user<lines;user++) {
				_delta_p[tone][line_id][user] = DBL_MAX;
			}
			return 1;
		}
		_delta_p[tone][line_id][user] = p[user] - old_p[user];
	}	
	
	return 0;
}




bool mipb::spectral_mask_constraint_broken(int tone,int user)
{
        if(_spectral_mask_on)
                if (_p[tone][user] + _delta_p[tone][user][user] > _spectral_mask)
                        return true;
        return false;
}




void mipb::init_lines()
{

        int tone,user;
        struct line* current;

        for (user=0;user<lines;user++) {
                current=get_line(user);
		if (_is_frac) {
			strcpy(current->loading_algo,"MIPB_FRAC");
			current->is_frac=true;
		}
		else {
			strcpy(current->loading_algo,"MIPB");
		}
                for (tone=0;tone<DMTCHANNELS;tone++) {
			if (_is_frac)
				current->_b[tone]=_b[tone][user];
                        current->b[tone]=_b[tone][user];
                        current->psd[tone]=watts_to_dbmhz(_p[tone][user]);
                }
        }


}

void mipb::write_current_stats(int tone,int user)
{

	char fn[80];
	char fn1[80];
	char fn2[80];
	FILE *fp=NULL;
	FILE *fp1=NULL;
	FILE *fp2=NULL;

	sprintf(fn,"%s/data/mipb/stats/%d.txt",ROOT_DIR,_total_bits);
	sprintf(fn1,"%s/data/mipb/stats/b_and_p_stats.txt",ROOT_DIR);
	sprintf(fn2,"%s/data/mipb/stats/cost_%d.txt",ROOT_DIR,_total_bits);

	fp = fopen(fn,"w");

	if (fp==NULL)
		exit(2);

	fp1 = fopen(fn1,"a");

	if (fp1==NULL)
		exit(2);

	fp2 = fopen(fn2,"w");

	if (fp2==NULL)
		exit(2);

	
	fprintf(fp,"#%dlines\n",lines);
	fprintf(fp,"#bit added was on tone %d line %d\n",tone,user);

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		//fprintf(fp,"%d\t",tone);
#ifdef VDSL_UPSTREAM
                if (tone == LOWERCHANNELS) {
                        for (int times=0;times<2;times++) {
                                if (times == 0)
                                        fprintf(fp,"%6.4lf\t",freq_matrix[tone-1]+CHANNEL_BANDWIDTH);
                                else
                                        fprintf(fp,"%6.4lf\t",freq_matrix[tone]-CHANNEL_BANDWIDTH);
                                for (int user=0;user<lines;user++) {
                                        fprintf(fp,"%d ",0);
                                }
                                for (int user=0;user<lines;user++) {
                                        fprintf(fp,"%s ","-inf");
                                }
                                fprintf(fp,"\n");
                        }
                }
#endif
		fprintf(fp,"%6.4lf\t",freq_matrix[tone]);
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%4.2lf ",_b[tone][user]);
		}
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%4.2lf ",watts_to_dbmhz(_p[tone][user]));
		}
		/*for (int user=0;user<lines;user++) {
			fprintf(fp,"%4.2e\t",cost[tone][user]);
		}
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%d\t",F[tone][user]);
		}
		*/
		if (tone < lines) {
			fprintf(fp,"%6.4lg\t",10*log10(_p_used[tone]*1e3));
			fprintf(fp,"%6.4lg\t",_wp[tone]);
			fprintf(fp,"%d\t",_b_total[tone]);
		}

		if (tone == 0) {
			//double p_ave=0.0;
			//for (int user=0;user<lines;user++) {
			//	p_ave+=p_used[user];
			//}
			//p_ave/=lines;
			fprintf(fp,"%6.4g",_p_ave*1e3);

		}
		
		fprintf(fp,"\n");
	}

	fprintf(fp1,"%d\t",_total_bits);

	for (int user=0;user<lines;user++) {
		fprintf(fp1,"%d   ",_b_total[user]);
	}

	for (int user=0;user<lines;user++) {
		fprintf(fp1,"%6.4lg   ",_p_used[user]);
	}

	fprintf(fp1,"\n");


	for (int tone=0;tone<DMTCHANNELS;tone++) {
		fprintf(fp2,"%d ",tone);
		for (int user=0;user<lines;user++) {
			fprintf(fp2,"%6.4g ",_cost[tone][user]);
		}
		fprintf(fp2,"\n");
	}



	fclose(fp);
	fclose(fp1);
	fclose(fp2);

}

void mipb::write_stats_totals()
{

        char fn[80];
        FILE *fp;

        sprintf(fn,"%s/data/mipb/stats/stats.txt",ROOT_DIR);

        fp = fopen(fn,"w");

        if (fp==NULL)
                exit(2);

        fprintf(fp,"#%dlines\n",lines);
        fprintf(fp,"#%dbits\n",_total_bits);
        fprintf(fp,"#%dtones\n",DMTCHANNELS);

        fclose(fp);

}

void mipb::dump_matrix(double **matrix,const char *tag)
{
	char fn[300];
	FILE *fp;

	sprintf(fn,"/tmp/mipb-%s.txt",tag);

	fp = fopen(fn,"w");

	if (fp==NULL)
		exit(2);

	for (int tone=0;tone<DMTCHANNELS;tone++) {
#ifdef VDSL_UPSTREAM
                if (tone == LOWERCHANNELS) {
                        for (int times=0;times<2;times++) {
                                if (times == 0) {
                                        fprintf(fp,"%6.4lf\t",freq_matrix[tone-1]+CHANNEL_BANDWIDTH);
				}
                                else {
                                        fprintf(fp,"%6.4lf\t",freq_matrix[tone]-CHANNEL_BANDWIDTH);
				}
                                
				for (int user=0;user<lines;user++) {
                                        fprintf(fp,"%d ",0);
                                }
				/*
                                for (int user=0;user<lines;user++) {
                                        fprintf(fp,"%s ","-inf");
                                }*/
                                fprintf(fp,"\n");
                        }
                }
#endif
		fprintf(fp,"%6.4lf\t",freq_matrix[tone]);
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%6.4g ",matrix[tone][user]);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

}

void mipb::dump_int_matrix(int **matrix,const char *tag)
{
	char fn[300];
	FILE *fp;

	sprintf(fn,"/tmp/mipb-%s.txt",tag);

	fp = fopen(fn,"w");

	if (fp==NULL)
		exit(2);

	for (int tone=0;tone<DMTCHANNELS;tone++) {
#ifdef VDSL_UPSTREAM
                if (tone == LOWERCHANNELS) {
                        for (int times=0;times<2;times++) {
                                if (times == 0)
                                        fprintf(fp,"%6.4lf\t",freq_matrix[tone-1]+CHANNEL_BANDWIDTH);
                                else
                                        fprintf(fp,"%6.4lf\t",freq_matrix[tone]-CHANNEL_BANDWIDTH);
                                for (int user=0;user<lines;user++) {
                                        fprintf(fp,"%d ",0);
                                }
				/*
                                for (int user=0;user<lines;user++) {
                                        fprintf(fp,"%s ","-inf");
                                }*/
                                fprintf(fp,"\n");
                        }
                }
#endif
		fprintf(fp,"%6.4lf\t",freq_matrix[tone]);
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%d ",matrix[tone][user]);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

}

/*
			for(int user=0;user<lines;user++) {
				if (_rate_targets[user] == NOT_SET)
					continue;

				_w_min[user]=0;
				_w_max[user]=1e-4;
				//_b_total[user]=0;
				last=-1;
				while(1) {
					reset_data();
					init_cost_matrix();
					load();				

					print_vector(_w,"_w");
					print_vector(_b_total,"_b_total");
					//getchar();

					if (_b_total[user] > _rate_targets[user]) {
						//printf("Setting _w_max[%d] = %lf\n",user,_w[user]);
						_w_max[user] = _w[user];
					}
					else if (_b_total[user] < _rate_targets[user]) {
						//printf("Setting _w_min[%d] = %lf\n",user,_w[user]);
						_w_min[user] = _w[user];
					}
					else {
						//printf("Fuck you.\n");
						//printf("_b_total[%d] = %d\n",user,_b_total[user]);
						//printf("_rate_targets[%d] = %d\n",user,_rate_targets[user]);
					}

					if (abs(_b_total[user] - _rate_targets[user]) < _rate_tol) {
						printf("Converged on line %d\n",user);
						break;
					}
				
					if (_w[user] == last) {
						printf("Converged on line %d but rate target not met\n",user);
						break;
						//getchar();
					}

					last=_w[user];
					_w[user] = (_w_max[user]+_w_min[user])/2;
					printf("changing weight vector to\n");	
					print_vector(_w,"_w");
					print_vector(_b_total,"_b_total");

				}
			}
			*/
			/*
			load();
			for (int user=0;user<lines;user++) {
				_natural_rate[user] = _b_total[user];
				printf("Natural rate on line %d = %d\n",user,_natural_rate[user]);
			}

			bool rates_too_big=true;
			for (int user=0;user<lines;user++) {
				if (_rate_targets[user] == NOT_SET)
					continue;
				if (abs(_natural_rate[user] - _rate_targets[user]) > _rate_tol && (_natural_rate[user] < _rate_targets[user])) {
					printf("Going to bisect on user %d\n",user);
					printf("Natural rate = %d Target = %d\n",_natural_rate[user],_rate_targets[user]);
					getchar();
					bisect_w(user,user);
					rates_too_big=false;
					break;	
				}
				if (_natural_rate[user] > _rate_targets[user] && _w[user] > 0) {
					printf("Going to bisect on user %d\n",user);
					printf("Natural rate = %d Target = %d\n",_natural_rate[user],_rate_targets[user]);
					getchar();
					bisect_w(user,user);
					rates_too_big=false;
					break;	
				}
			}

			if (rates_too_big) {
				printf("Rates are all too big\n");
				for (int user=0;user<lines;user++) {
					if (_rate_targets[user] == NOT_SET) {
						printf("Rate target on line %d is not set\n",user);
						for (int rateuser=0;rateuser<lines;rateuser++ ) {
							if (_rate_targets[rateuser] == NOT_SET)
								continue;
							if (_natural_rate[rateuser] > _rate_targets[rateuser]) {
								printf("Selected users weight to bisect is %d\n",rateuser);
								bisect_w(user,rateuser);
								break;
							}
						}
						break;
					}
				}
			}
			*/
			/*
			reset_data();
			init_cost_matrix();
			load();
			print_vector(_w,"_w");
			print_vector(_b_total,"_b_total");
			if (abs(_b_total[user] - _rate_targets[user]) < _rate_tol) {
					printf("Converged on line %d\n",user);
					break;
			}
			*/

void mipb::add_to_cost_list(struct cost_entry *c)
{

	if (_cost_list_head==NULL) {
		_cost_list_head = new cost_entry;
		_cost_list_tail = _cost_list_head;
	
		_cost_list_head->user=c->user;
		_cost_list_head->tone=c->tone;
		_cost_list_head->cost=c->cost;

		return;
	}

	cost_entry *current = new cost_entry;

	current->user=c->user;
	current->tone=c->tone;
	current->cost=c->cost;
	
	current->prev=_cost_list_tail;
	current->prev->next=current;
	_cost_list_tail=current;
	return;

}

void mipb::add_to_cost_list_after(struct cost_entry *c,struct cost_entry *current) 
{
	struct cost_entry *n = new cost_entry;	// new entry
	n->user=c->user;
	n->tone=c->tone;
	n->cost=c->cost;

	if (current==_cost_list_tail) {		// adding right at the end
		current->next=n;
		n->prev=current;
		_cost_list_tail=n;
		//cout << "Added one to the tail" << endl;
		return;
	}

	n->next=current->next;	// link in new entry
	n->prev=current;

	current->next->prev=n;	// fix next entry in the list prev porinter
	current->next=n;	// fix current next pointer to now point to n

	assert(n->next !=NULL);
	return;

}

void mipb::add_to_cost_list_before(struct cost_entry *c,struct cost_entry *current) 
{
	struct cost_entry *n = new cost_entry;	// new entry
	n->user=c->user;
	n->tone=c->tone;
	n->cost=c->cost;

	if (current==_cost_list_head) {		// adding right at the front
		current->prev=n;
		n->next=current;
		_cost_list_head=n;
		//cout << "Added one to the tail" << endl;
		return;
	}

	n->next=current;	// link in new entry
	n->prev=current->prev;

	current->prev->next=n;	// fix next entry in the list prev porinter
	current->prev=n;	// fix current next pointer to now point to n

	//assert(n->next !=NULL);
	return;

}


void mipb::add_to_front(struct cost_entry *c)
{
	struct cost_entry *n = new cost_entry;	// new entry
	n->user=c->user;
	n->tone=c->tone;
	n->cost=c->cost;

	_cost_list_head->prev=n;	// link old head to new entry
	n->next=_cost_list_head;	// set up new entries next to point to old head
	
	_cost_list_head=n;
}

void mipb::print_cost_list()
{
	struct cost_entry *current = _cost_list_head;
	int entries=0;

	while (current!=NULL) {
		entries++;
		//cout << *current << endl;
		if (current->next != NULL) {
			if (current->cost > current->next->cost) {
				cout << "fuck up" << endl;
				getchar();
			}
		}
		current=current->next;
	}
	cout << "Entries = " << entries << endl;
}

void mipb::remove_cost_list_head()
{
	struct cost_entry *current=_cost_list_head;

	if (current->next == NULL) {
		print_cost_list();
	}

	_cost_list_head=current->next;
	_cost_list_head->prev=NULL;
	delete current;

	return;

}

void mipb::cost_list_del()
{
	struct cost_entry *current=_cost_list_head;
	struct cost_entry *temp;

	while (current!=NULL){
		temp=current;
		delete temp;
		current=current->next;
	}

}

double mipb::cf2(int tone,int line_id)
{

	double cost=0.0;
	double a=0;
	double b=0;
	double c[lines];

	if (total_power_constraint_broken(tone,line_id) || spectral_mask_constraint_broken(tone,line_id)) {
		_F[tone][line_id] = 1;
		return DBL_MAX;
	}

	if (_F[tone][line_id] == 1)
		return DBL_MAX;

	if (_b[tone][line_id] >= MAXBITSPERTONE) {
		_F[tone][line_id]=1;
		return DBL_MAX;
	}

	for (int user=0;user<lines;user++) {
		if (user == line_id) {
			if (_old) {
				cost+=_delta_p[tone][line_id][user];		// power for the actual bit loaded
			}
			else if (!_greedy) {
				a+=cost+=_wp[user]*_delta_p[tone][line_id][user];
			}
			else {
				c[line_id]=cost+=_w[line_id]*_delta_p[tone][line_id][user];		// power for the actual bit loaded
			}
		}
		else {
			/*if (_b[tone][user] == 0) {	// victim line has no bits yet
				//printf("Victim has no bits yet\n");
				if (_no_fext_bits[tone][user] > 1) { 	// victim can theoretically use this tone
					int b[lines];
					double p[lines];
					double p_old[lines];
					double g[lines];
					//printf("Line looking for cost is %d\n",line_id);
					//printf("Victim is %d\n",user);
					for (int user1=0;user1<lines;user1++) {
						g[user1] = line_array[user1]->gamma[line_array[user1]->service[tone]];
					}
					memset(p,0,sizeof(double)*lines);
					memset(p_old,0,sizeof(double)*lines);
					memcpy(b,&_b[tone][0],sizeof(int)*lines);	// get current b vector
					//b[user]=(int)_ref_bits[tone][user];					// pretend victim has one bit
					//b[user]=MIN((int)_no_fext_bits[tone][user],_w[user]*_ref_bits[tone][user]);
					b[user]=1;
					_psd->calc(b,g,tone,p_old);				// get psd vector
					b[line_id]++;		
					_psd->calc(b,g,tone,p);				// get psd vector

					cost+=_wp[user]*(p[user]-p_old[user]);
		
					if (_wp[user]*(p[user]-p_old[user]) < 0) {
						printf("Line_id = %d\n",line_id);
						printf("victim = %d\n",user);
						printf("tone = %d\n",tone);
						print_vector(_b[tone],"b_current");
						print_vector(p_old,"p_old");
						print_vector(b,"b");
						print_vector(p,"p");
						print_vector(_wp,"wp");
						getchar();
					}
					//printf("Extra cost incurred is %6.4g\n",p[user]-p_old[user]);
					
				}
			}*/
			//cost+=_wxt[user]*_delta_p[tone][line_id][user];		// power for the extra crosstalk incurred
			//cost+=_wp[user]*_delta_p[tone][line_id][user]+_mask[tone][line_id][user];		// power for the extra crosstalk incurred
			if(_old) {
				cost+=_delta_p[tone][line_id][user]; 	// power for the extra crosstalk incurred
			}
			else if (!_greedy) {
				b+=cost+=_wp[user]*_delta_p[tone][line_id][user];
			}
			else {
				//if ((1-_w[user]) < 0)
				cost+=c[user]=(1-_w[user])*_delta_p[tone][line_id][user]; 	// power for the extra crosstalk incurred
				//cost+=c[user]=_delta_p[tone][line_id][user]; 	// power for the extra crosstalk incurred
				//cost+=c[user]=(exp(-0.01*_w[user]))*_wp[user]*_delta_p[tone][line_id][user]; 	// power for the extra crosstalk incurred
				//cost+=_w[line_id]*_wp[user]*_delta_p[tone][line_id][user]; 	// power for the extra crosstalk incurred
			}
		}	
	
	}

	/*	
	if (b > a && _wp[line_id] < 0.9) {
		printf("Cost for bit on line %d tone %d is dominated by crosstalk terms even though wp is less than 1\n",line_id,tone);
		printf("_wp[%d] = %g\ta = %g\tb = %g\n",line_id,_wp[line_id],a,b);
		print_vector(_p_used,"_p_used");
		printf("Total bits = %d\n",_total_bits);
		getchar();
	}*/

	//printf("Total cost is %6.4g\n",cost);
	//getchar();
	//

	if (b > 2*a) {
		printf("Crosstalk twice price of bit on line %d tone %d\n",line_id,tone);
		print_vector(_wp,"wp");
		print_vector(_b_total,"b_total");
	}

	if (cost == 0) {
		printf("WTF\n");
		printf("Line_id = %d\n",line_id);
		printf("tone = %d\n",tone);
		//print_vector(_b[tone],"b_current");
		//print_vector(p_old,"p_old");
		//print_vector(b,"b");
		//print_vector(p,"p");
		print_vector(_w,"w");
		getchar();
	}

	/*
	for (int user=0;user<lines;user++) {
		if (c[user] < 0) {
			printf("Line_id = %d\n",line_id);
			print_vector(c,"c");
			print_vector(_w,"w");
			print_vector(_wp,"wp");
			print_vector(_delta_p[tone][line_id],"delta_p");
			print_vector(_b[tone],"b");
			double ratio = c[line_id]/(-1*c[user]);
			printf("Ratio = %6.4g\n",ratio);
			if (ratio < 10)
				getchar();
			break;
		}
	}
	*/

	if (_old)
		return _w[line_id]*cost;
	else 
		return cost;
}
