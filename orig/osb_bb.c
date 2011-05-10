#include "multiuser_load.h"
#include "osb_bb.h"
#include <cassert>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include "util.h"
static double per_tone_opt_time=0;
static double total_opt_time=0;
static double n_opt=0;
static bool sol_debug=false;
static bool list_debug=false;
static void *opt_p_top_half(void *p);
static void *opt_p_bottom_half(void *p);
static void *opt_p_thread(void * p);
struct thread_args{
	osb_bb *o;
	int lower;
	int upper;	
	int thread_id;
};
osb_bb::osb_bb(const int t)
{
	_log = fopen("osb_bb.log","w");
	if (_log == NULL) {
		printf("Cannot open log file\n");
		exit(2);
	}
	_b = new int*[DMTCHANNELS];
	_p = new double*[DMTCHANNELS];
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		_b[tone] = new int[lines];
		_p[tone] = new double[lines];
	}
	_l = new double[lines];
	_best_l = new double[lines];
	_w = new double[lines];
	_w_max = new double[lines];
	_w_min = new double[lines];
	_rate_targets = new int[lines];
	_p_distance=DBL_MAX;
	_prev_p_distance=DBL_MAX;
	_min_p_distance=DBL_MAX;
	_l_evals=0;
	_p_budget=new double[lines];
	for (int user=0;user<lines;user++) {
		_l[user]=0;
		_w[user]=1/(double)lines;
		_w_max[user]=100;
		_w_min[user]=0;
		_p_budget[user]=0.110;
		_rate_targets[user]=NOT_SET;
	}
	_num_threads=t;
	_psd = new psd_vector*[_num_threads];
	for (int n=0;n<_num_threads;n++) {
		_psd[n] = new psd_vector;
	}
	_min_sl=500;
	_p_tol=0.015;
	_dynamic_lambda=false;
	_spectral_mask=dbmhz_to_watts(100000);
	_sw=1e-5;
	_branch_and_bound=true;
	// initialisation of _b and _p to shut up valgrind
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			_b[tone][user] = 0;
			_p[tone][user] = 0;
		}
	}
}
osb_bb::~osb_bb()
{
	fclose(_log);
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		delete[] _b[tone];
		delete[] _p[tone];
	}
	delete[] _b;
	delete[] _p;
	delete[] _l;
	delete[] _w;
	delete[] _p_budget;
	for (int n=0;n<_num_threads;n++) {
		delete _psd[n];
	}
	delete _psd;
}
int osb_bb::run()
{
	bool _started_dynamic=false;
	if (_dynamic_lambda) {
		_started_dynamic=true;
		_min_sl=1;
	}
	_sl=_min_sl;		// set step size in here in case we changed the default value
	bool target=false;
	for (int user=0;user<lines;user++) {
		if (_rate_targets[user] != NOT_SET)
			target=true;
	}
	if (target) {
		do {
			for(int user=0;user<lines;user++) {
				if (_rate_targets[user] == NOT_SET)
					continue;
				print_vector(_w,"current weights before bisecting rate on line");
				printf("%d\n",user);
				bisect_l();
				if (abs(tot_rate(user) - _rate_targets[user]) < _rate_tol) {
					printf("Converged very early on line %d\n",user);
					continue;
				}
				bool early=false;
				if (tot_rate(user) > _rate_targets[user]) {	// w is max
					printf("Intial w[%d] = %lf was max\n",user,_w[user]);
					_w_max[user]=_w[user];
					_w[user]/=2;					// now find min
					while (1) {
						bisect_l();
						if (tot_rate(user) > _rate_targets[user]) {
							_w_max[user]=_w[user];
							_w[user]/=2;
						}
						if (tot_rate(user) < _rate_targets[user]) {
							_w_min[user]=_w[user];
							printf("Found the min on line %d = %lf\n",user,_w[user]);
							break;
						}
						if (abs(tot_rate(user) - _rate_targets[user]) < _rate_tol) {
							printf("Converged early on line %d\n",user);
							early=true;
							break;
						}
					}
				}
				else if (tot_rate(user) < _rate_targets[user]) {	// w is min
					printf("Intial w[%d] = %lf was min\n",user,_w[user]);
					_w_min[user]=_w[user];
					_w[user]*=2;					// now find max
					while (1) {
						bisect_l();
						if (tot_rate(user) < _rate_targets[user]) {
							_w_min[user]=_w[user];
							_w[user]*=2;
						}
						if (tot_rate(user) > _rate_targets[user]) {
							printf("Found the max on line %d = %lf\n",user,_w[user]);
							_w_max[user]=_w[user];
							break;
						}
						if (abs(tot_rate(user) - _rate_targets[user]) < _rate_tol) {
							printf("Converged early on line %d\n",user);
							early=true;
							break;
						}
					}
				}	
				if (early)
					continue;
				_w[user] = (_w_max[user]+_w_min[user])/2;
				printf("About to start rate bisection on line %d starting with w[%d] = %lf\n",user,user,_w[user]);	
				while(1) {
					bisect_l();
					printf("Number of lambda evaluations = %d\n",_l_evals);
					print_vector(_w,"w");
					//update_w();
					print_vector(_w,"w");
					if (tot_rate(user) > _rate_targets[user]) {
						_w_max[user] = _w[user];
					}
					else if (tot_rate(user) < _rate_targets[user]) {
						_w_min[user] = _w[user];
					}
					if (abs(tot_rate(user) - _rate_targets[user]) < _rate_tol) {
						printf("Rate bisection converged on line %d\n",user);
						//getchar();
						break;
					}
					_w[user] = (_w_max[user]+_w_min[user])/2;
					//getchar();
				}
				if (rates_converged())
					break;
			}
		}while (!rates_converged());
	}
	else {
		bisect_l();
	}
	init_lines();
	calculate_snr();
	return 0;
}
int osb_bb::bisect_l()
{
	double *l_last = new double[lines];
	double *p_current = new double[lines];
	while (1) {	
		for (int user=0;user<lines;user++) {
			double l_min;
			double l_max;
			double pow,last;
			_l[user]=0;
			while(1) {
				//print_vector(_l,"_l");
				optimise_p();
				if (tot_pow(user) > _p_budget[user]) {
					if (_l[user] ==0) {
						_l[user]++;
					}
					else {
						_l[user]*=2;
					}
				}
				else 
					break;
			} 
			//printf("l_max[%d] = %lf\n",user,_l[user]);
			l_max=_l[user];
			l_min=_l[user]/2;
			while(1) {
				_l[user] = (l_max + l_min)/2;
				//print_vector(_l,"_l");
				optimise_p();
				pow = tot_pow(user);
				if (pow > _p_budget[user]) {
					l_min = _l[user];
				}
				else {
					l_max = _l[user];
				}
				if (pow==last || (fabs(_p_budget[user]-pow) < _p_tol)) {
					//printf("we're done on user %d\n",user);
					break;
				}
				last=pow;
			}
		}
		//print_vector(_l,"_l");
		//print_vector(l_last,"l_last");
		for (int user=0;user<lines;user++) {
			p_current[user] = tot_pow(user);
		}
		//print_vector(p_current,"p_current");
		bool done=true;
		for (int user=0;user<lines;user++) {
			if (l_last[user] != _l[user]) {
				done=false;
			}
		}
		if (done) {
			//printf("The whole lot has converged\n");
			break;
		}
		else {
			//printf("Still going\n");
		}
		memcpy(l_last,_l,sizeof(double)*lines);
		if (all_powers_within_tol()) {
			break;
		}
	}
	//init_lines();
	//calculate_snr();
	return 0;
}
void osb_bb::update_sl()
{
	static double min_p_distance_last,p_distance_last_stuck;
	static bool stuck;
	if (_l_evals==0) {
		stuck=false;
		min_p_distance_last=DBL_MAX;
		p_distance_last_stuck=DBL_MAX;
	}
	if (stuck) {
		_sl=_min_sl;
		stuck=false;
		return;
	}
	_p_distance = calc_p_distance();
	if (_p_distance <= _prev_p_distance) {		// new distance is better than last one
		printf("Last distance was %lf, current distance is %lf\n",_prev_p_distance,_p_distance);
		if (_p_distance <= _min_p_distance) {	// new distance is best so far
			printf("Min distance was %lf, new min distance is %lf\n",_min_p_distance,_p_distance);
			_min_p_distance = _p_distance;
			for (int user=0;user<lines;user++) {
				_best_l[user] = _l[user];	// record best l vector
			}
		}
		_sl*=2;
		printf("Step size increased to %lf\n",_sl);
	}
	else {
		stuck=false;
		if (_min_p_distance == min_p_distance_last) {
                        printf("Looks like we are stuck\n");
                        //getchar();
                        _dynamic_lambda=false;
			_sl=100;
                        /*if (min_p_distance == p_distance_last_stuck) {
				stuck++;
			}*/
			/*
			if (_min_sl < 0.01) {
                        	printf("Really stuck\n");
				getchar();
			}*/
			stuck=true;
                }
                min_p_distance_last = _min_p_distance;
		printf("Starting a new trajectory\n");
		for (int user=0;user<lines;user++) {
			_l[user] = _best_l[user];
		}
		if (stuck) 
			_sl = rand()%1000;
		else
			_sl = _min_sl;
		_prev_p_distance = DBL_MAX;
		return;
	}	
	_prev_p_distance = _p_distance;
	return;
}
void osb_bb::optimise_p()
{
	/*
	static bool first_run=true;
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		// create our intial region for each tone, contains all values
		regions *rs = new regions;
		if (tone==0 && first_run) {
			int *b_init = new int[lines];
			first_run=false;
			for (int user=0;user<lines;user++) {
				b_init[user]=8;
			}
			branch_regions(rs,tone,b_init);
		}
		else if (tone ==0 && !first_run) {
			branch_regions(rs,tone,&_b[0][0]);	// this is the value from before we update _l
		}
		else {
			// use the last solution as a starting point to speed up convergance
			branch_regions(rs,tone,&_b[tone-1][0]);
		}
		delete rs;
		//printf("This tone ave opt time was %lf seconds\n",per_tone_opt_time/DMTCHANNELS);
		printf("Total average opt time is %lf seconds\n",total_opt_time/n_opt);
		per_tone_opt_time=0;
	}
	*/
	/*
	pthread_t bottom_half,top_half;
	struct thread_args th1,th2;
	th1.o = this;
	th2.o = this;	
	pthread_create(&bottom_half,NULL,opt_p_bottom_half,(void *)&th1);
	pthread_create(&top_half,NULL,opt_p_top_half,(void *)this);
	pthread_join(bottom_half,NULL);
	pthread_join(top_half,NULL);
	*/
	struct timeval time_a,time_b;
	gettimeofday(&time_a,NULL);
	pthread_t th[_num_threads];
	struct thread_args args[_num_threads];
	for (int t=0;t<_num_threads;t++) {
		args[t].o = this;
		args[t].lower = t *(DMTCHANNELS/_num_threads);
		args[t].upper = (t+1) *(DMTCHANNELS/_num_threads);
		args[t].thread_id=t;
		pthread_create(&(th[t]),NULL,opt_p_thread,(void *)&(args[t]));
	}
	for (int t=0;t<_num_threads;t++) {
		pthread_join(th[t],NULL);
	}
	gettimeofday(&time_b,NULL);
	double t = ((time_b.tv_sec - time_a.tv_sec)*1e6 + (time_b.tv_usec - time_a.tv_usec))/1e6;
	if (_branch_and_bound) {
		//printf("average time taken for branch and bound = %lf seconds\n",t/DMTCHANNELS);
	}
	else {
		//printf("average time taken for per_tone_opt = %lf seconds\n",t/DMTCHANNELS);
	}
	//exit(0);
}
void *opt_p_thread(void * p)
{
	thread_args *th = static_cast<thread_args *>(p);
	osb_bb *o = static_cast<osb_bb *>(th->o);	
	//printf("Hello from thread %d!\n",th->thread_id);	
	//printf("My lower bound is %d\n",th->lower);	
	//printf("My upper bound is %d\n",th->upper);
	//printf("My psd pointer is %p\n",o->_psd[th->thread_id]);
	//getchar();	
	for (int tone=th->lower;tone<th->upper;tone++) {
	//for (int tone=134;tone<136;tone++) {
		if (o->_branch_and_bound) {
			regions *rs = new regions;
			//printf("tone = %d\n",tone);
			//int *b = new int[lines];
			//o->per_tone_exhaustive_search(tone,o->_psd[th->thread_id]);
			//memcpy(b,o->_b[tone],sizeof(int)*lines);
			//print_vector(o->_b[tone],"b");
			//print_vector(o->_p[tone],"p");
			o->branch_regions(rs,tone,o->_psd[th->thread_id]);
			//if (memcmp(b,o->_b[tone],sizeof(int)*lines) != 0) {
			//	print_vector(o->_b[tone],"bb");
			//	print_vector(b,"osb");
			//	getchar();
			//}
			//print_vector(o->_b[tone],"b");
			//print_vector(o->_p[tone],"p");
			//getchar();
			//printf("Done on tone %d by thread %d\n",tone,th->thread_id);
			//delete[] b;
			delete rs;
		}
		else {
			o->per_tone_exhaustive_search(tone,o->_psd[th->thread_id]);
		}
	}
	return (void *)NULL;	// silence warning
}
/*
void *opt_p_bottom_half(void * p)
{
	printf("Hello from bottom half thread!\n");	
	for (int tone=0;tone<DMTCHANNELS/2;tone++) {
		regions *rs = new regions;
		osb_bb *o = static_cast<osb_bb *>(p);	
		o->branch_regions(rs,tone,o->_psd1);
		delete rs;
	}
}
void *opt_p_top_half(void *p)
{
	printf("Hello from top half thread!\n");	
	for (int tone=DMTCHANNELS/2;tone<DMTCHANNELS;tone++) {
		regions *rs = new regions;
		osb_bb *o = static_cast<osb_bb *>(p);
		o->branch_regions(rs,tone,o->_psd2);
		delete rs;
	}
}
*/
region *osb_bb::find_initial_max_low(int tone,psd_vector* psd)
{
	double *p_max = new double[lines];
	double *p_min = new double[lines];
	double max_low;	
	double *g = new double[lines];
	region *ret= new region;
	for (int user=0;user<lines;user++) {	// get gamma values on this tone
		g[user] = line_array[user]->gamma[0];
	}
	region *isb_sol = new region;
	isb_solution(tone,psd,isb_sol->b_min);
	memcpy(isb_sol->b_max,isb_sol->b_min,sizeof(int)*lines);
	psd->calc(isb_sol->b_min,g,tone,p_min);
	psd->calc(isb_sol->b_max,g,tone,p_max);
	isb_sol->calc_lower_bound(_w,_l,_spectral_mask,p_max);
	isb_sol->calc_upper_bound(_w,_l,_spectral_mask,p_min);	
	max_low=isb_sol->lower;
	region::copy(ret,isb_sol);
	for (int i=0;i<3;i++) {
		int init_tone=tone+(i-1);
		if (init_tone < 0 || init_tone >=DMTCHANNELS)
			continue;
		region *init = new region;
		memcpy(init->b_min,_b[init_tone],sizeof(int)*lines);
		memcpy(init->b_max,_b[init_tone],sizeof(int)*lines);
		psd->calc(init->b_min,g,tone,p_min);
	        psd->calc(init->b_max,g,tone,p_max);
        	init->calc_lower_bound(_w,_l,_spectral_mask,p_max);
	        init->calc_upper_bound(_w,_l,_spectral_mask,p_min);	
		if (init->lower > max_low) {
			max_low = init->lower;
			region::copy(ret,init);
		}
		delete init;
	}
	delete isb_sol;
	delete[] p_max;
	delete[] p_min;
	delete[] g;
	//ret->print_region();
	//getchar();
	return ret;
}
void osb_bb::branch_regions(regions *rs,int tone,psd_vector *psd)
{
	//struct timeval time_a,time_b;
	//gettimeofday(&time_a,NULL);
	region *r = rs->r_list;
	double max_low=-DBL_MAX;
	bool max_low_changed=false;
	double *p_max = new double[lines];
	double *p_min = new double[lines];
	double *g = new double[lines];
	for (int user=0;user<lines;user++) {	// get gamma values on this tone
		g[user] = line_array[user]->gamma[0];
	}
	// find initial max_low here
	region *init=find_initial_max_low(tone,psd);	
	max_low=init->lower;
	//init->print_region();
	//getchar();
	/*	
	region *init = new region;
	if ((tone==0 || ((tone % (DMTCHANNELS/N_THREADS)) == 0)) && _l_evals==0) {
		isb_solution(tone,psd,init->b_min);
		//printf("isb solution - ");
		//print_vector(init->b_min,"init->b_min");
		//memcpy(init->b_max,init->b_min,sizeof(int)*lines);
		//psd->calc(init->b_min,g,tone,p_min);
		//psd->calc(init->b_max,g,tone,p_max);
		//init->calc_lower_bound(_w,_l,_spectral_mask,p_max);
		//init->calc_upper_bound(_w,_l,_spectral_mask,p_min);
		//max_low=init->lower;
		//printf("Max low = %lf\n",max_low);
		//printf("First run on tone %d\n",tone);
		//calculate_b_vector_flat_psd(dbmhz_to_watts(-40),g,tone,init->b_min);
		memcpy(init->b_max,init->b_min,sizeof(int)*lines);
	}
	else {
		int init_tone;
		if (tone == 0) {
			init_tone=tone;
		}
		else {
			init_tone=tone-1;
		}
		memcpy(init->b_min,_b[init_tone],sizeof(int)*lines);
		memcpy(init->b_max,_b[init_tone],sizeof(int)*lines);
	}
	if (sol_debug || list_debug) {
		printf("Tone %d\n",tone);
	}
	psd->calc(init->b_min,g,tone,p_min);
	psd->calc(init->b_max,g,tone,p_max);
	init->calc_lower_bound(_w,_l,_spectral_mask,p_max);
	init->calc_upper_bound(_w,_l,_spectral_mask,p_min);
	//printf("Initial region\n");
	//init->print_region();
	//printf("Press a key\n");
	//getchar();	
	max_low=init->lower;
	while (max_low==-DBL_MAX) {	// if first init resulted in spectral mask breakage or cocktail party
		//printf("That b_init didnt work\n");
		//print_vector_psd(p_max,"p_max");
		for (int user=0;user<lines;user++) {
			if (p_max[user] > _spectral_mask || p_max[user] < 0) {
				//printf("Reducing b[%d]\n",user);
				init->b_min[user]--;
				init->b_max[user]--;
				break;
			}
		}
		psd->calc(init->b_max,g,tone,p_max);
		psd->calc(init->b_min,g,tone,p_min);
		init->calc_lower_bound(_w,_l,_spectral_mask,p_max);
		init->calc_upper_bound(_w,_l,_spectral_mask,p_min);
		max_low=init->lower;
		//printf("Press a key\n");
		//getchar();
	}
	*/
	//if (max_low > MAXBITSPERTONE * lines) {
	//	printf("WTF??\n");
	//	getchar();
	//}
	//init->print_region();
	//print_vector_psd(p_max,"p");
	//printf("init lower bound = %.50lf\n",max_low);
	//getchar();
	while(1) {
		if (check_region_size(r)) {
			//if (count_temp_list(r) > 1) {
				//printf("All regions are unit size but there is more than one solution....\n");
				//print_temp_list(r);
			//}
			rs->r_list=r;
			memcpy(&_b[tone][0],r->b_max,sizeof(int)*lines);
			psd->calc(r->b_max,g,tone,p_max);
			memcpy(&_p[tone][0],p_max,sizeof(double)*lines);
			if (sol_debug) {
				printf("Solution is ");
				r->print_region();
				print_vector_psd(p_max,"p");
				printf("Press a key\n");
				getchar();
			}
			goto out;
		}
		// temp list pointers which are linked back to rs->r_list on every branching
		region *list_head=NULL;
		region *list_tail=NULL;
		int initial_discards=0;
		int recheck_discards=0;
		int new_regions=0;		
		//printf("Current list is \n");
		//getchar();
		//print_temp_list(r);	
		//count_temp_list(r);
		while (r != NULL) {	//iterate over all current regions
			unsigned int n = (int)pow(2,lines);	// number of new regions 2^N
			//printf("Region to split = \n");
			//r->print_region();
			//getchar();
			for (unsigned int j=0;j<n;j++) {		// for all new regions
				region *temp = new region;
				for (int user=0;user<lines;user++) {  // for each user in the new region
					if (j&(1<<user)) {	// top half of old region
						temp->b_max[user] = r->b_max[user];
						temp->b_min[user] = (r->b_max[user]+r->b_min[user])/2 + 1;
					}
					else {			// bottom half
						temp->b_max[user] = (r->b_max[user]+r->b_min[user])/2;
						temp->b_min[user] = r->b_min[user];
					}	
					if (temp->b_max[user] < temp->b_min[user]) {
						printf("WTF??");
						temp->print_region();
						printf("Region split was \n");
						r->print_region();
						getchar();
					}
				}
				// calc upper bound of new region
				psd->calc(temp->b_min,g,tone,p_min);
				temp->calc_upper_bound(_w,_l,_spectral_mask,p_min);
				if (temp->upper < max_low) {	// discard the region straight away
					if (region::region_is_within_region(init,temp) && (max_low_changed==false)) {
						printf("Fuck up occured\n");						
						printf("Tone = %d\n",tone);
						getchar();
						print_vector(temp->b_min,"current b_min");
						psd->calc(temp->b_min,g,tone,p_min);
						print_vector(p_min,"current p_min");
						temp->calc_upper_bound(_w,_l,_spectral_mask,p_min);
						printf("current region\n");
						temp->print_region();
						printf("init region\n");
						init->print_region();	
					}
					//printf("Discarding region\n");
					//temp->print_region();
					initial_discards++;
					delete temp;
					continue;	// next region please
				}
				new_regions++;
				// calc lower bound of new region
				psd->calc(temp->b_max,g,tone,p_max);
				temp->calc_lower_bound(_w,_l,_spectral_mask,p_max);
				if (temp->lower >= max_low) {	// ah-ha! a new max_low!
					//printf("New max low found\n");
					//temp->print_region();
					if(!max_low_changed) {
						max_low_changed=true;
					}
					max_low=temp->lower;
				}
				if (list_head == NULL) {
					list_head=temp;
					list_tail=temp;
				}
				else {
					temp->last=list_tail;
					temp->last->next=temp;
					list_tail=temp;
				}
			}
			//getchar();
			//print_temp_list(list_head);
			//getchar();	
			// move to the next region to split and delete the old one
			region *t = r;
			r=r->next;			
			delete t;
		}	
		// after iterating over the entire region list, recheck all the old entries incase a new max low was found
		region *current=list_head;
		while(current!=NULL) {
			region *next=current->next;
			if (current->upper < max_low) {
				recheck_discards++;
				new_regions--;
				if (current==list_head) {
					list_head=del_from_temp_list(current);
					//printf("Deleted a region from the list head\n");
				}
				else if (current==list_tail) {
					list_tail=del_from_temp_list(current);
					//printf("Deleted a region from the list tail\n");
				}
				else {
					del_from_temp_list(current);
				}
			}
			current=next;
		}
		if (list_debug) {
			printf("%d regions intially discarded\n",initial_discards);
			printf("%d regions discarded on recheck\n",recheck_discards);
			printf("New list has %d regions\n",new_regions);
			count_temp_list(list_head);
		}
		// link the temp_list back into r to start a new branching level
		if (new_regions==0) {
			printf("Ahhhh! no regions in list!\n");
			if (!max_low_changed) {
				printf("The max low never changed after init\n");
			}
			printf("Max low = %.50lf\n",max_low);
			printf("Tone = %d\n",tone);
			printf("Region =\n");
			init->print_region();
			for (int user=0;user<lines;user++) {
				printf("l[%d] = %.60lf\n",user,_l[user]);
			}
			getchar();
		}
		r=list_head;
	}
	out:
	delete init;
	delete[] p_max;
	delete[] p_min;
	delete[] g;
	//gettimeofday(&time_b,NULL);
	//double t = ((time_b.tv_sec - time_a.tv_sec)*1e6 + (time_b.tv_usec - time_a.tv_usec))/1e6;
	//printf("Time taken for branch and bound on tone %d is %lf seconds\n",tone,t);
}
void osb_bb::update_l()
{
	for (int user=0;user<lines;user++) {
                double pow = tot_pow(user);
                double update = _sl*(pow-_p_budget[user]);
                if (_l[user] + update < 0)
                        _l[user] = _l[user]*0.1;
		else
			_l[user] = _l[user] + update;
	}
}
void osb_bb::update_w()
{
	for (int user=0;user<lines;user++) {
		if (_rate_targets[user] == NOT_SET)
			continue;
		int rate = tot_rate(user);
                double update = _sw*(rate-_rate_targets[user]);
                if (_w[user] - update < 0)
                        _w[user] = _w[user]*0.1;
		else
			_w[user] = _w[user] - update;
	}
}
bool osb_bb::powers_converged()
{
	if (all_powers_within_tol())
                return true;
        else
                return false;
}
bool osb_bb::rates_converged()
{
	printf("Checking rates\n");
	for (int user=0;user<lines;user++)
		tot_rate(user);		
	if (all_rates_within_tol())
		return true;
	else
		return false;
}
void osb_bb::init_lines()
{
        struct line *current;
        for (int user=0;user<lines;user++) {
                current=get_line(user);
                current->is_dual=0;
                strcpy(current->loading_algo,"OSB BB");
                for (int tone=0;tone<DMTCHANNELS;tone++) {
                        current->b[tone] = _b[tone][user];
                        current->psd[tone] = watts_to_dbmhz(_p[tone][user]);
                }
        }
}
bool osb_bb::all_powers_within_tol()
{
 	for (int user=0;user<lines;user++) {
		if (_l[user] < 1e-20) {	// if l equals 0 it is a non active power constraint
			continue;
		}
        	if (fabs(tot_pow(user) - _p_budget[user]) > _p_tol) {
			return false;
		}
	}
	return true;
}
bool osb_bb::all_rates_within_tol()
{
	for (int user=0;user<lines;user++) {
		if (_rate_targets[user] == NOT_SET)
			continue;
		if (abs(tot_rate(user) - _rate_targets[user]) > _rate_tol) {
			return false;
		}
	}
	return true;
}
int osb_bb::tot_rate(int user)
{
	int sum=0;
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		sum+=_b[tone][user];
	}
	printf("Rate of user %d = %d\n",user,sum);
	return sum;
}
double osb_bb::tot_pow(int user)
{
	double sum=0;
        for (int tone=0;tone<DMTCHANNELS;tone++) {
                sum+=_p[tone][user];
        }
        //printf("power used by user %d = %16.14lf\n",user,sum);
        fprintf(_log,"power used by user %d = %16.14lf\n",user,sum);
        return sum;
}
/*
region *osb_bb::add_to_temp_list(region * tail,region *r)
{
}
*/
region *osb_bb::del_from_temp_list(region *current)
{
	if (current->last == NULL) {	// first entry in list
		region *ret=current->next;	// this will be the new list_head
		current->next->last=NULL;
		delete current;
		return ret;
	}
	if (current->next == NULL) {	// last entry in list
		region *ret=current->last;	// this will be new list tail
		current->last->next=NULL;
		delete current;
		return ret;
	}
	current->last->next=current->next;
	current->next->last=current->last;
	delete current;
	return NULL;
}
void osb_bb::print_temp_list(region *head)
{
	region *current=head;
	while(current!=NULL) {
		current->print_region();
		current=current->next;
	}
}
int osb_bb::count_temp_list(region *head)
{
	region *current=head;
	int n=0;
	while(current!=NULL) {
		n++;
		current=current->next;
	}
	printf("List has %d members\n",n);
	return n;
}
bool osb_bb::check_region_size(region * head)
{
	region *current=head;
	//assert(current!=NULL);
	if (current == NULL) {
		printf("Whoops!! the list is empty\n");
		printf("Please investigate\n");
		getchar();
	}
	while (current!=NULL) {
		if (current->lower != current->upper) {
			return false;
		}
		current=current->next;
	}
	return true;
}
double osb_bb::calc_p_distance()
{
        double sum=0.0;
        for (int user=0;user<lines;user++) {
                sum += pow((_p_budget[user] - tot_pow(user)),2);
        }
	printf("Current p_distance is %lf\n",sqrt(sum));
        return sqrt(sum);
}
void osb_bb::isb_solution(int tone,psd_vector *psd,int *b)
{
	double lk,lk_max;
	bool converged;
	int b_max,b_max_2;
	int *b_last = new int [lines];
	for (int user=0;user<lines;user++) {
		b[user]=0;
	}
	while(1) {
		for (int user=0;user<lines;user++) {
			b_last[user]=b[user];
		}	
		for (int user=0;user<lines;user++) {
			for (b[user]=0,lk=0.0,lk_max=-DBL_MAX;b[user]<=MAXBITSPERTONE;b[user]++) {
				lk = l_k(b,tone,psd);
				if (lk>lk_max) {
					lk_max=lk;
					b_max=b[user];
				}
				else {
					b_max=b[user]-1;
					break;
				}
			}
			b[user]=b_max;
		}
		converged=true;
		for (int user=0;user<lines;user++) {
			if (b_last[user]!=b[user]) {
				converged=false;
				break;
			}
		}
		if (converged) {
			//print_vector(b,"isb b");
			//getchar();
			delete[] b_last;
			break;
		}
	}
}
void osb_bb::per_tone_exhaustive_search(int tone, psd_vector *psd)
{
	double lk=0.0;
	double lk_max=-DBL_MAX;
	int *b_max = new int[lines];
	double *g = new double[lines];
	double *p = new double[lines];
	for (int user=0;user<lines;user++) {    // get gamma values on this tone
                g[user] = line_array[user]->gamma[0];
        }
	vector_next(_b[tone],lines,MAXBITSPERTONE,true);
	do {
		lk = l_k(_b[tone],tone,psd);
		//if (tone == 62) {
		//	printf("lk = %lf\n",lk);
		//	print_vector(_b[tone],"b");
		//}
		if (lk >= lk_max) {
			lk_max=lk;
			//if (tone == 62) 
			//	printf("New lk_max = %lf\n",lk_max);
			for (int user=0;user<lines;user++) {
                 		b_max[user] = _b[tone][user];
                        }
		}
		//getchar();
	}while(vector_next(_b[tone],lines,MAXBITSPERTONE,false));
	//if (tone == 62) {
	//	getchar();
	//}
	psd->calc(b_max,g,tone,_p[tone]);
	for (int user=0;user<lines;user++) {
		_b[tone][user] = b_max[user];
	}
	delete[] b_max;
	delete[] g;
	delete[] p;
}
double osb_bb::l_k(int *b,int tone,psd_vector *psd)
{
        double b_sum=0,p_sum=0;
        double *p = new double [lines];
        double *g = new double [lines];
        for (int user=0;user<lines;user++) {
                b_sum+=b[user]*_w[user];
		g[user]=line_array[user]->gamma[0];
        }
        /*calculate_psd_vector(_b,g,tone,p);*/
        psd->calc(b,g,tone,p);
        for (int user=0;user<lines;user++) {
                if ((p[user] < 0) || (p[user] > _spectral_mask)) {
                      /*
                      p[user]=DBL_MAX;
                      if (p[user] > _spectral_mask) {
	                      printf("Spectral mask broken\n");
                              getchar();
                      }                                                                                                                                                                 */
                	delete[] p;
			delete[] g;
                        return -DBL_MAX;
                }
                p_sum+=_l[user]*p[user];
                /*last_psd[user]=p[user];*/
        }
        delete[] p;
        delete[] g;
        return b_sum-p_sum;
}
