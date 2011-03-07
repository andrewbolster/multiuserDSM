#include "multiuser_load.h"
#include "greedy_nsection.h"
#include "mipb.h"

#include <boost/threadpool.hpp>
#include <boost/bind.hpp>
#include <boost/mem_fn.hpp>
#include <boost/bind/mem_fn.hpp>

using namespace boost::threadpool;

greedy_nsection::greedy_nsection(const int t)
{
	_num_threads=t;

	_greedy_workers = new mipb*[_num_threads];

	_rate_targets = new int[lines];
	_w = new double[lines];
	_first_w_run = new bool[lines];
	_p_budget = new double[lines];
	
	_tp = new pool(_num_threads);

	for (int user=0;user<lines;user++) {
		_w[user] = 1;
		_first_w_run[user]=true;
		_p_budget[user] = 0.111;
	}

	_rate_tol=80;

	_nsection=true;
	_grad_search=false;
}


void greedy_nsection::nsection()
{

	
	for (int thread_id=0;thread_id<_num_threads;thread_id++) {
		_greedy_workers[thread_id] = new mipb(1);
		_greedy_workers[thread_id]->_rate_tol=_rate_tol;
		_greedy_workers[thread_id]->_graph_loading=false;
		_greedy_workers[thread_id]->_greedy=true;
		_greedy_workers[thread_id]->_old=true;
		_greedy_workers[thread_id]->_simple_search=false;
		_greedy_workers[thread_id]->_search1=true;
		_greedy_workers[thread_id]->_search2=false;

		for (int user=0;user<lines;user++) {
			_greedy_workers[thread_id]->_rate_targets[user]=_rate_targets[user];
			_greedy_workers[thread_id]->_p_budget[user]=_p_budget[user];
		}
	}
	int ngreedys=0;
	bool done=false;
	int winning_thread;
	
	if (_nsection) {
		assert(_num_threads >= 2);
		printf("Hi!\n");
		while(1) {

			for (int user=0;user<lines;user++) {
				if (_rate_targets[user] == NOT_SET)
					continue;
				
				double max=DBL_MAX;
				int max_rate=-INT_MAX;
				double min=0;
				int min_rate=INT_MAX;
				int max_thread;
				int min_thread;

				//double w_range_bottom=0.001;
				double w_range_bottom;
				//double w_range_top=10;
				double w_range_top;
				//_w[user] = 0.001;
				//
				

				if (_first_w_run[0]) {
					//_first_w_run[0]=false;
					w_range_bottom = 1/(pow((double)2,(double)((double)_num_threads/4)));
					//w_range_top = 1*(pow((double)2,(double)((double)_num_threads/4)));
					//w_range_top = _num_threads*0.5;
					w_range_top = pow(_num_threads,0.5);
				}
				else {
					
					print_vector(_w,"w");
					print_vector(_greedy_workers[winning_thread]->_w,"winning thread w");
					print_vector(_greedy_workers[winning_thread]->_b_total,"winning thread b_total");
					if (_greedy_workers[winning_thread]->_b_total[user] > _rate_targets[user]) {
						printf("Winning thread b[%d] = %lf which is greater than rate_target\n",user,_greedy_workers[winning_thread]->_b_total[user]);
						printf("Setting bottom to %lf\n",_w[user]);
						printf("Setting top to %lf\n",_w[user]*_num_threads);
						w_range_bottom=_w[user];
						//w_range_top=_w[user]*_num_threads+1;
						//w_range_top=_w[user]*_num_threads;
						w_range_top=_w[user]*2;
					}
					else if (_greedy_workers[winning_thread]->_b_total[user] < _rate_targets[user])	{
						printf("Winning thread b[%d] = %lf which is less than rate_target\n",user,_greedy_workers[winning_thread]->_b_total[user]);
						printf("Setting top to %lf\n",_w[user]);
						printf("Setting bottom to %lf\n",_w[user]/_num_threads);
						w_range_top=_w[user];
						//w_range_bottom=_w[user]/_num_threads-1;
						//w_range_bottom=_w[user]/_num_threads;
						w_range_bottom=_w[user]/2;
					}
					if(abs(_greedy_workers[winning_thread]->_b_total[user]-_rate_targets[user]) < _rate_tol){
						printf("Wow w already converged\n");
						_w[user]=_greedy_workers[winning_thread]->_w[user];
						continue;
					}
				}
		
		
				while (1) {
					printf("Threads = %d\ttop = %lf\tbottom = %lf\n",_num_threads,w_range_top,w_range_bottom);
					for (int thread_id=0;thread_id<_num_threads;thread_id++) {
						//_greedy_workers[thread_id]->_w[user] = _w[user]*(thread_id+1);
						memcpy(_greedy_workers[thread_id]->_w,_w,sizeof(double)*lines);
						/*
						if (_num_threads==2) {
							_greedy_workers[thread_id]->_w[user] = w_range_bottom+((w_range_top-w_range_bottom)/(_num_threads+1))*(thread_id);
						}
						else {
							_greedy_workers[thread_id]->_w[user] = w_range_bottom+((w_range_top-w_range_bottom)/(_num_threads-1))*(thread_id);
						}
						*/
						if (_first_w_run[0] || _num_threads==2)
							_greedy_workers[thread_id]->_w[user] = w_range_bottom+((w_range_top-w_range_bottom)/(_num_threads+1))*(thread_id+1);
						else 
							_greedy_workers[thread_id]->_w[user] = w_range_bottom+((w_range_top-w_range_bottom)/(_num_threads-1))*(thread_id);
							
						printf("Thread %d w = %lf\n",thread_id,_greedy_workers[thread_id]->_w[user]);
						//printf("%d %lf\n",_num_threads,_greedy_workers[thread_id]->_w[user]);
						_tp->schedule(boost::bind(boost::mem_fn(&mipb::greedy_load),_greedy_workers[thread_id]));	
					}
					//getchar();
					ngreedys++;
					_tp->wait();
					
					for (int thread_id=0;thread_id<_num_threads;thread_id++) {
						print_vector(_greedy_workers[thread_id]->_b_total,"b_total");	
						print_vector(_greedy_workers[thread_id]->_w,"w");	
						if (_greedy_workers[thread_id]->_b_total[user] > _rate_targets[user]) {
							if(_greedy_workers[thread_id]->_b_total[user] < min_rate) {
								min=_greedy_workers[thread_id]->_w[user];
								min_rate=_greedy_workers[thread_id]->_b_total[user];
								min_thread=thread_id;
							} 	
						}
						else if (_greedy_workers[thread_id]->_b_total[user] < _rate_targets[user]) {
							if (_greedy_workers[thread_id]->_b_total[user] > max_rate) {
								max=_greedy_workers[thread_id]->_w[user];
								max_rate=_greedy_workers[thread_id]->_b_total[user];
								max_thread=thread_id;
							}
						}
					}


					printf("max on line %d = %lf rate = %d\n",user,max,max_rate);
					printf("min on line %d = %lf rate = %d\n",user,min,min_rate);

					//getchar();
					if (max == DBL_MAX) {
						if (abs(rate_targets[user]-min_rate) < _rate_tol) {
							_w[user]=min;
							winning_thread=min_thread;
							break;
						}
						w_range_top*=2;
						w_range_bottom=min;
						printf("Restarting because we havent found a max w on line %d\n",user);
						continue;
					}

					if (min == 0) {
						if (abs(rate_targets[user]-max_rate) < _rate_tol) {
							_w[user]=max;
							winning_thread=max_thread;
							break;
						}
						w_range_top=max;
						w_range_bottom/=2;
						//w_range_top=max;
						printf("Restarting because we havent found a min w on line %d\n",user);
						continue;
					}	
			
					if ( (_rate_targets[user]-max_rate) < (min_rate-_rate_targets[user]) ) {	// max is closer than min
						if (abs(rate_targets[user]-max_rate) < _rate_tol) {
							_w[user]=max;
							winning_thread=max_thread;
							break;
						}
					}
					else {
						if (abs(rate_targets[user]-min_rate) < _rate_tol) {
							_w[user]=min;
							winning_thread=min_thread;
							break;
						}		
					}
					

					w_range_bottom=min;
					w_range_top=max;

				}

				if (_first_w_run[0])
					_first_w_run[0]=false;

				printf("bisection of w complete on line %d = %lf\n",user,_w[user]);

				if (_greedy_workers[winning_thread]->rates_converged()) {
					printf("Sweet!\n");
					_greedy_workers[winning_thread]->init_lines();
					done=true;
					break;
				}

			}
			if (done)
				break;

		}
	}

	else if (_grad_search) {

		int increments=0;

		for (int thread_id=0;thread_id<_num_threads;thread_id++) {
			increments+=(thread_id+1)*(thread_id+1);
		}

		
		for (int thread_id=0;thread_id<_num_threads;thread_id++) {	// each worker expect the first starts with randomised w vector
			if (thread_id==0)
				_tp->schedule(boost::bind(boost::mem_fn(&mipb::greedy_load),_greedy_workers[thread_id]));
			/*else {
				for (int user=0;user<lines;user++) {
					if (_rate_targets[user] != NOT_SET)
						_greedy_workers[thread_id]->_w[user]=((double)rand())/RAND_MAX;
				}
				_tp->schedule(boost::bind(boost::mem_fn(&mipb::greedy_load),_greedy_workers[thread_id]));
			}*/
		}


		ngreedys++;
		_tp->wait();
		

		bool first_run=true;

		//double min_step;
		double *b_total_last = new double[lines];
		while (1) {
		
		

			double min=DBL_MAX;
			double r_distance;
			for (int thread_id=0;thread_id<_num_threads;thread_id++) {	
				printf("Thread %d results\n",thread_id);	
				print_vector(_greedy_workers[thread_id]->_b_total,"b_total");	
				print_vector(_greedy_workers[thread_id]->_w,"w");
				r_distance = _greedy_workers[thread_id]->calc_r_distance();
				printf("R distance = %lf\n",r_distance);
				
				if (r_distance < min) {
					min = r_distance;
					winning_thread=thread_id; 
				}
				
			}

			printf("Winning thread was %d\n",winning_thread);

			
			


			if (first_run) {
				first_run=false;
			}
			else {
				for (int user=0;user<lines;user++) {
					if (rate_targets[user] == NOT_SET)
						continue;
					if (_greedy_workers[winning_thread]->_b_total[user] > _rate_targets[user] &&
						b_total_last[user] < _rate_targets[user]) {
						printf("Oscillation detected\n");
					}
					
					if (_greedy_workers[winning_thread]->_b_total[user] < _rate_targets[user] &&
						b_total_last[user] > _rate_targets[user]) {
						printf("Oscillation detected\n");
					}
				}
				memcpy(b_total_last,_greedy_workers[winning_thread]->_b_total,sizeof(double)*lines);	
			}


			if (_greedy_workers[winning_thread]->rates_converged()) {
				printf("Sweet!\n");
				_greedy_workers[winning_thread]->init_lines();
				break;	
			}

			//print_vector(_greedy_workers[winning_thread]->_b_total,"winning thread b_total");	
			/*	
			double max_step=DBL_MAX;
			double step;

			for (int user=0;user<lines;user++) {
				if (_rate_targets[user] == NOT_SET)
					continue;
				if (_greedy_workers[winning_thread]->_b_total[user] < _rate_targets[user]) {
					step = (-1*_greedy_workers[winning_thread]->_w[user])/(_greedy_workers[winning_thread]->_b_total[user] - _rate_targets[user]);
					if (step < max_step)
						max_step=step;
				}
					 
			}

			if (max_step==DBL_MAX) {
				max_step=0.005;
			}
	
			printf("Absolute Max step = %lf\n",max_step);
			//max_step=max_step*(1-exp(-0.04*_num_threads));
			max_step/=10;
			printf("Max step used = %lf\n",max_step);
			
			//min_step=max_step/increments; 
			double min_step=max_step/_num_threads; 
			//min_step=max_step/128; 
			*/


			double min_step=0.0001;
			//double max_step=0.005;

			/*
			double inc;
			if (_num_threads>1)
				inc=max_step-min_step/(_num_threads-1);
			else
				inc=0;
			*/

			for (int thread_id=0;thread_id<_num_threads;thread_id++) {
				//print_vector(_greedy_workers[winning_thread]->_b_total,"winning thread b_total");	
				//printf("Thread %d\n",thread_id);
				//printf("winning thread = %d\n",winning_thread);	
				//print_vector(_greedy_workers[thread_id]->_w,"w_before");
				for (int user=0;user<lines;user++) {
					if (_rate_targets[user] != NOT_SET) {
						//_greedy_workers[thread_id]->_w[user] = _greedy_workers[winning_thread]->_w[user] + (_greedy_workers[winning_thread]->_b_total[user] - _rate_targets[user])*min_step*(thread_id+1)*(thread_id+1);
						
						//print_vector(_greedy_workers[winning_thread]->_b_total,"winning thread b_total");	
						//printf("Changing w on thread %d user %d\n",thread_id,user);
						//printf("w = %lf\n",_greedy_workers[thread_id]->_w[user]);
						//printf("winning thread b_total[%d] = %lf\n",user,_greedy_workers[winning_thread]->_b_total[user]);
						//printf("rate target[%d] = %d\n",user,_rate_targets[user]);
						_greedy_workers[thread_id]->_w[user] = _greedy_workers[winning_thread]->_w[user] + (_greedy_workers[winning_thread]->_b_total[user] - _rate_targets[user])*min_step*(thread_id+1);
						//_greedy_workers[thread_id]->_w[user] = _greedy_workers[winning_thread]->_w[user] + (_greedy_workers[winning_thread]->_b_total[user] - _rate_targets[user])*(inc*(thread_id)+min_step);
						//printf("w_new = %lf\n",_greedy_workers[thread_id]->_w[user]);
						if (_greedy_workers[thread_id]->_w[user] < 0) {
							_greedy_workers[thread_id]->_w[user] = 1e-10;
						}
					}
				}
				//print_vector(_greedy_workers[thread_id]->_w,"w_after");
				//_tp->schedule(boost::bind(boost::mem_fn(&mipb::greedy_load),_greedy_workers[thread_id]));	
			}
			ngreedys++;

			for (int thread_id=0;thread_id<_num_threads;thread_id++) {
				_tp->schedule(boost::bind(boost::mem_fn(&mipb::greedy_load),_greedy_workers[thread_id]));	
			}

			_tp->wait();	

		}

	}

	printf("Number of threaded greedy evals = %d\n",ngreedys);
 
}
