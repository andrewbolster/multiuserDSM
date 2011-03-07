#include "multiuser_load.h"
#include "osb_bb_frac.h"
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
static bool init_debug=false;
static bool new_region_debug=false;

static void *opt_p_top_half(void *p);
static void *opt_p_bottom_half(void *p);
static void *opt_p_thread(void * p);


struct thread_args{
	osb_bb_frac *o;
	int lower;
	int upper;	
	int thread_id;
	struct tone_lock *locks;
};

struct tone_lock{
	tone_lock()
	{
		done=false;
		pthread_mutex_init(&mutex,NULL);
	}
	bool done;
	pthread_mutex_t mutex;
};


osb_bb_frac::osb_bb_frac()
{

	_log = fopen("osb_bb_frac.log","w");

	if (_log == NULL) {
		printf("Cannot open log file\n");
		exit(2);
	}

	_b = new double*[DMTCHANNELS];
	_p = new double*[DMTCHANNELS];

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		_b[tone] = new double[lines];
		_p[tone] = new double[lines];
	}

	_l = new double[lines];
	_best_l = new double[lines];
	_w = new double[lines];

	_p_distance=DBL_MAX;
	_prev_p_distance=DBL_MAX;
	_min_p_distance=DBL_MAX;

	_l_evals=0;

	_p_budget=new double[lines];

	for (int user=0;user<lines;user++) {
		_l[user]=0;
		_w[user]=1;
		_p_budget[user]=0.110;
	}

	//_psd = new psd_vector*[N_THREADS];

	//for (int n=0;n<N_THREADS;n++) {
	//	_psd[n] = new psd_vector;
	//}

	_min_sl=500;
	_p_tol=0.015;

	_dynamic_lambda=false;

	_spectral_mask=dbmhz_to_watts(-30);

	// initialisation of _b and _p to shut up valgrind
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			_b[tone][user] = 0;
			_p[tone][user] = 0;
		}
	}

	_inc=0.1;

}

osb_bb_frac::~osb_bb_frac()
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

	/*
	for (int n=0;n<N_THREADS;n++) {
		delete _psd[n];
	}
	delete[] _psd;
	*/

}

int osb_bb_frac::run()
{

	if (_dynamic_lambda)
		_min_sl=1;

	_sl=_min_sl;		// set step size in here in case we changed the default value

	while(!converged() || _l_evals == 0) {
		print_vector(_l,"l");
		//print_vector_to_file(_l,"l",_log);
		optimise_p();
		
		if (_dynamic_lambda)
			update_sl();
		_l_evals++;
		//fprintf(_log,"Number of lambda evaluations = %d\n",_l_evals);
		update_l();
		//getchar();
		printf("Current distance = %.15lf\n",calc_p_distance());
	}

	printf("Number of lambda evaluations = %d\n",_l_evals);

	init_lines();

	calculate_snr();

	return 0;
}

int osb_bb_frac::bisect_l()
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

				if (pow==last) {
					//printf("we're done on user %d\n",user);
					break;
				}

				last=pow;
			}

		}

		print_vector(_l,"_l");
		print_vector(l_last,"l_last");
		
		for (int user=0;user<lines;user++) {
			p_current[user] = tot_pow(user);
		}

		print_vector(p_current,"p_current");

		bool done=true;
		for (int user=0;user<lines;user++) {
			if (l_last[user] != _l[user]) {
				done=false;
			}
		}
		if (done) {
			printf("The whole lot has converged\n");
			break;
		}
		else {
			printf("Still going\n");
		}
		memcpy(l_last,_l,sizeof(double)*lines);

		if (all_powers_within_tol()) {
			break;
		}

	}
	init_lines();

	calculate_snr();

	return 0;

}


void osb_bb_frac::update_sl()
{
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
		printf("Starting a new trajectory\n");
		for (int user=0;user<lines;user++) {
			_l[user] = _best_l[user];
		}
		_sl = _min_sl;
		_prev_p_distance = DBL_MAX;
		return;
	}	
	
	_prev_p_distance = _p_distance;
	return;
}

void osb_bb_frac::optimise_p()
{

	/*
	static bool first_run=true;

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		// create our intial region_frac for each tone, contains all values
		region_fracs_frac *rs = new regions;

		if (tone==0 && first_run) {
			int *b_init = new int[lines];
			first_run=false;
			for (int user=0;user<lines;user++) {
				b_init[user]=8;
			}
			
			branch_region_fracs_frac(rs,tone,b_init);
		}
		else if (tone ==0 && !first_run) {
			branch_region_fracs_frac(rs,tone,&_b[0][0]);	// this is the value from before we update _l
		}
		else {
			// use the last solution as a starting point to speed up convergance
			branch_region_fracs_frac(rs,tone,&_b[tone-1][0]);
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

	struct tone_lock lock_array[DMTCHANNELS];
	
	pthread_t th[N_THREADS];
	struct thread_args args[N_THREADS];

	for (int t=0;t<N_THREADS;t++) {
		args[t].o = this;
		args[t].lower = t *(DMTCHANNELS/N_THREADS);
		args[t].upper = (t+1) *(DMTCHANNELS/N_THREADS);
		args[t].thread_id=t;
		args[t].locks=lock_array;

		pthread_create(&(th[t]),NULL,opt_p_thread,(void *)&(args[t]));
	}

	for (int t=0;t<N_THREADS;t++) {
		pthread_join(th[t],NULL);
	}


	gettimeofday(&time_b,NULL);

	double t = ((time_b.tv_sec - time_a.tv_sec)*1e6 + (time_b.tv_usec - time_a.tv_usec))/1e6;

	printf("average time taken for branch and bound = %lf seconds\n",t/DMTCHANNELS);

}

void *opt_p_thread(void * p)
{

	thread_args *th = static_cast<thread_args *>(p);
	osb_bb_frac *o = static_cast<osb_bb_frac *>(th->o);	

	printf("Hello from thread %d!\n",th->thread_id);	
	//printf("My lower bound is %d\n",th->lower);	
	//printf("My upper bound is %d\n",th->upper);
	//printf("My psd pointer is %p\n",o->_psd[th->thread_id]);
	//getchar();	
	/*
	while(1) {
		bool all_done=true;
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			if (th->locks[tone].done == false) {
				all_done=false;
				if (pthread_mutex_trylock(&(th->locks[tone].mutex)) == 0) {
					regions_frac *rs = new regions_frac;
					printf("Thread %d tone %d\n",th->thread_id,tone);
					o->branch_regions(rs,tone);
					delete rs;
					th->locks[tone].done=true;
					pthread_mutex_unlock(&(th->locks[tone].mutex));
				}
				else {
					continue;
				}
			}
		}
		if (all_done) {
			printf("All done!\n");
			break;
		}
	}
	*/
	
	for (int tone=th->lower;tone<th->upper;tone++) {
	//for (int tone=134;tone<136;tone++) {
		regions_frac *rs = new regions_frac;
		//o->branch_regions(rs,tone,o->_psd[th->thread_id]);
		o->branch_regions(rs,tone);
		//printf("Done on tone %d\n",tone);
		delete rs;
	}
	printf("Thread %d done\n",th->thread_id);

	return (void *)NULL;	// silence warning
}

/*
void *opt_p_bottom_half(void * p)
{
	printf("Hello from bottom half thread!\n");	

	for (int tone=0;tone<DMTCHANNELS/2;tone++) {
		region_fracs_frac *rs = new regions;
		osb_bb *o = static_cast<osb_bb *>(p);	
		o->branch_region_fracs_frac(rs,tone,o->_psd1);
		
		delete rs;
	}


}

void *opt_p_top_half(void *p)
{

	printf("Hello from top half thread!\n");	
	
	for (int tone=DMTCHANNELS/2;tone<DMTCHANNELS;tone++) {
		region_fracs_frac *rs = new regions;
		osb_bb *o = static_cast<osb_bb *>(p);
		o->branch_region_fracs_frac(rs,tone,o->_psd2);
		delete rs;

	}
}
*/


void osb_bb_frac::branch_regions(regions_frac *rs,int tone)
{

	
	//struct timeval time_a,time_b;

	//gettimeofday(&time_a,NULL);

	region_frac *r = rs->r_list;

	double max_low=-DBL_MAX;
	bool max_low_changed=false;
	double *p_max = new double[lines];
	double *p_min = new double[lines];
	double *g = new double[lines];

	for (int user=0;user<lines;user++) {	// get gamma values on this tone
		g[user] = line_array[user]->gamma[0];
	}

	// find initial max_low here
	
	/*if (tone == 51) {
		list_debug=true;
		sol_debug=true;
		init_debug=true;
		new_region_debug=true;
	}*/

	region_frac *init = new region_frac;
	
	if ((tone==0 || ((tone % (DMTCHANNELS/N_THREADS)) == 0)) && _l_evals==0) {
		//printf("First run on tone %d\n",tone);
		calculate_b_vector_flat_psd(dbmhz_to_watts(-40),g,tone,init->b_min);
		memcpy(init->b_max,init->b_min,sizeof(double)*lines);
	}
	else {
		int init_tone;
		if (tone == 0) {
			init_tone=tone;
		}
		else {
			init_tone=tone-1;
		}
	
		memcpy(init->b_min,_b[init_tone],sizeof(double)*lines);
		memcpy(init->b_max,_b[init_tone],sizeof(double)*lines);

		if (init_debug) {
			printf("Init tone = %d\n",init_tone);
		}

	}

	if (sol_debug || list_debug) {
		printf("Tone %d\n",tone);
	}

	calculate_psd_vector(init->b_min,g,tone,p_min);
	calculate_psd_vector(init->b_max,g,tone,p_max);
	init->calc_lower_bound(_w,_l,_spectral_mask,p_max);
	init->calc_upper_bound(_w,_l,_spectral_mask,p_min);
	//printf("Initial region_frac\n");
	//init->print_region_frac();
	//printf("Press a key\n");
	//getchar();	
	//
		
	if (init_debug) {
		printf("Tone = %d\n",tone);
		printf("Init region\n");
		init->print_region_frac();
	}
	
	max_low=init->lower;
	while (max_low==-DBL_MAX) {	// if first init resulted in spectral mask breakage or cocktail party
		if (any_neg(init->b_min)) {
			printf("Oh shit\n");
			getchar();
		}
		if (init_debug) {
			printf("That b_init didnt work\n");
			print_vector_psd(p_max,"p_max");
		}
		for (int user=0;user<lines;user++) {
			//if (p_max[user] > _spectral_mask || p_max[user] < 0) {
			if (init_debug) {
				printf("Reducing b[%d]\n",user);
			}
			if (init->b_min[user] >= _inc) {
				if (init_debug) {
					printf("init->b_min[%d] = %50.48lf\n",user,init->b_min[user]);
					printf("%d\n",(int)(init->b_min[user]/_inc));
					printf("%50.48lf\n",(double)((int)(init->b_min[user]/_inc))*_inc);
					printf("%50.48lf\n",(double)((int)(init->b_min[user]/_inc))*_inc-_inc);
				}		
				init->b_min[user] = (double)((int)(init->b_min[user]/_inc))*_inc-_inc;
			}
			if (init->b_max[user] >= _inc) {
				if (init_debug) {
					printf("init->b_max[%d] = %50.48lf\n",user,init->b_max[user]);
                                        printf("%d\n",(int)(init->b_max[user]/_inc));
                                        printf("%50.48lf\n",(double)((int)(init->b_max[user]/_inc))*_inc);
                                        printf("%50.48lf\n",(double)((int)(init->b_max[user]/_inc))*_inc-_inc);
				}
				init->b_max[user] = (double)((int)(init->b_max[user]/_inc))*_inc-_inc;
			}
		}
		if (init_debug) {
			printf("About to re try to calculate the bounds of the initial region on tone %d with reduced b\n",tone);
			print_vector_hires(init->b_min,"b_min");
			print_vector_hires(init->b_max,"b_max");
		}
		calculate_psd_vector(init->b_max,g,tone,p_max);
		calculate_psd_vector(init->b_min,g,tone,p_min);
		if (init_debug) {
			print_vector(p_min,"p_min");	
			print_vector(p_max,"p_max");
			print_vector_psd(p_max,"p");
			print_vector(init->b_min,"b_min");
			print_vector(init->b_max,"b_max");
		}
		init->calc_lower_bound(_w,_l,_spectral_mask,p_max);
		init->calc_upper_bound(_w,_l,_spectral_mask,p_min);
		if (init_debug) {
			init->print_region_frac();
			getchar();
		}
		max_low=init->lower;
		if (init_debug) {
			printf("Press a key\n");
			getchar();
		}
	}

	

	//if (max_low > MAXBITSPERTONE * lines) {
	//	printf("WTF??\n");
	//	getchar();
	//}
	if (init_debug) {
		init->print_region_frac();
		print_vector_psd(p_max,"p");
		printf("init lower bound = %.50lf\n",max_low);
		getchar();
	}

	while(1) {

		if (check_region_frac_size(r)) {
			//if (count_temp_list(r) > 1) {
				//printf("All region_fracs_frac are unit size but there is more than one solution....\n");
				//print_temp_list(r);
			//}
			rs->r_list=r;
			memcpy(&_b[tone][0],r->b_max,sizeof(double)*lines);
			calculate_psd_vector(r->b_max,g,tone,p_max);
			memcpy(&_p[tone][0],p_max,sizeof(double)*lines);
			if (sol_debug) {
				printf("Solution is ");
				r->print_region_frac();
				print_vector_psd(p_max,"p");
				printf("Press a key\n");
				getchar();
			}
			goto out;
		}
		
		// temp list pointers which are linked back to rs->r_list on every branching
		region_frac *list_head=NULL;
		region_frac *list_tail=NULL;

		int initial_discards=0;
		int recheck_discards=0;
		int new_regions=0;		

		//printf("Current list is \n");
		//getchar();
		//print_temp_list(r);	
		//count_temp_list(r);

		while (r != NULL) {	//iterate over all current region_fracs_frac
			unsigned int n = (int)pow(2,lines);	// number of new region_fracs_frac 2^N

			//printf("Region to split = \n");
			//r->print_region_frac();
			//getchar();
			
			for (unsigned int j=0;j<n;j++) {		// for all new region_fracs_frac
				region_frac *temp = new region_frac;
				
				for (int user=0;user<lines;user++) {  // for each user in the new region_frac
					if (r->b_max[user] == r->b_min[user]) { // dont split this dimension
						temp->b_max[user] = r->b_max[user];
						temp->b_min[user] = r->b_min[user];
						continue;
					}
					if (j&(1<<user)) {	// top half of old region_frac
						//double t;
						//printf("Top half of old region\n");
						temp->b_max[user] = r->b_max[user];
						//temp->b_min[user] = (r->b_max[user]+r->b_min[user]) + 0.1;
						temp->b_min[user] = (double)((int)(((r->b_max[user]+r->b_min[user])/2)*(1/_inc)))/(1/_inc) + _inc;
					}
					else {			// bottom half
						//printf("Bottom half of old region\n");
						temp->b_max[user] = (double)((int)(((r->b_max[user]+r->b_min[user])/2)*(1/_inc)))/(1/_inc);
						temp->b_min[user] = r->b_min[user];
					}	
	
					if (temp->b_max[user] < temp->b_min[user]) {
						printf("WTF??");
						temp->print_region_frac();
						printf("Region split was \n");
						r->print_region_frac();
						getchar();
					}
				}

				// calc upper bound of new region_frac
				
				/*if ((tone == 51) && (temp->b_max[0] > 8.7) && (temp->b_min[1] < 6.2)) {
					temp->print_region_frac();
					printf("Is this the right region?\n");
					getchar();
				}*/
				calculate_psd_vector(temp->b_min,g,tone,p_min);
				temp->calc_upper_bound(_w,_l,_spectral_mask,p_min);

				//printf("New region_frac\n");
				//temp->print_region_frac();
				//getchar();

				if (temp->upper < max_low) {	// discard the region_frac straight away
					if (new_region_debug) {
						printf("Discarding region_frac\n");
						temp->print_region_frac();
					}
					initial_discards++;
					delete temp;
					continue;	// next region_frac please
				}

				new_regions++;

				// calc lower bound of new region_frac
				calculate_psd_vector(temp->b_max,g,tone,p_max);
				temp->calc_lower_bound(_w,_l,_spectral_mask,p_max);

				if (new_region_debug) {
					printf("New region\n");
					temp->print_region_frac();
					getchar();
				}

				if (temp->lower > max_low) {	// ah-ha! a new max_low!
					//printf("New max low found\n");
					//temp->print_region_frac();
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
			// move to the next region_frac to split and delete the old one
			region_frac *t = r;
			r=r->next;			
			delete t;
		}	


		// after iterating over the entire region_frac list, recheck all the old entries incase a new max low was found
		region_frac *current=list_head;
		while(current!=NULL) {
			region_frac *next=current->next;
			if (current->upper < max_low) {
				recheck_discards++;
				new_regions--;
				if (current==list_head) {
					list_head=del_from_temp_list(current);
					//printf("Deleted a region_frac from the list head\n");
				}
				else if (current==list_tail) {
					list_tail=del_from_temp_list(current);
					//printf("Deleted a region_frac from the list tail\n");
				}
				else {
					del_from_temp_list(current);
				}
			}
			current=next;
		}

		if (list_debug) {
			printf("%d regions initially discarded\n",initial_discards);
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
			init->print_region_frac();
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

void osb_bb_frac::update_l()
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

bool osb_bb_frac::converged()
{
	if (all_powers_within_tol())
                return true;
        else
                return false;
}

void osb_bb_frac::init_lines()
{

        struct line *current;

        for (int user=0;user<lines;user++) {
                current=get_line(user);
                current->is_dual=0;
		current->is_frac=true;
                strcpy(current->loading_algo,"OSB BB FRAC");

                for (int tone=0;tone<DMTCHANNELS;tone++) {
                        current->b[tone] = _b[tone][user];
                        current->_b[tone] = _b[tone][user];
                        current->psd[tone] = watts_to_dbmhz(_p[tone][user]);
                }
        }


}

bool osb_bb_frac::all_powers_within_tol()
{

 	for (int user=0;user<lines;user++) {
		if (_l[user] < 1e-100) {	// if l equals 0 it is a non active power constraint
			continue;
		}
        	if (fabs(tot_pow(user) - _p_budget[user]) > _p_tol) {
			return false;
		}
	}

	return true;

}

double osb_bb_frac::tot_pow(int user)
{

	double sum=0;

        for (int tone=0;tone<DMTCHANNELS;tone++) {
                sum+=_p[tone][user];
        }

        printf("power used by user %d = %16.14lf\n",user,sum);
        fprintf(_log,"power used by user %d = %16.14lf\n",user,sum);

        return sum;


}
/*
region_frac *osb_bb::add_to_temp_list(region * tail,region *r)
{
}
*/

region_frac *osb_bb_frac::del_from_temp_list(region_frac *current)
{

	if (current->last == NULL) {	// first entry in list
		region_frac *ret=current->next;	// this will be the new list_head
		current->next->last=NULL;
		delete current;
		return ret;
	}
	if (current->next == NULL) {	// last entry in list
		region_frac *ret=current->last;	// this will be new list tail
		current->last->next=NULL;
		delete current;
		return ret;
	}

	current->last->next=current->next;
	current->next->last=current->last;

	delete current;

	return NULL;
}

void osb_bb_frac::print_temp_list(region_frac *head)
{
	region_frac *current=head;

	while(current!=NULL) {
		current->print_region_frac();
		current=current->next;
	}

}

int osb_bb_frac::count_temp_list(region_frac *head)
{
	region_frac *current=head;
	int n=0;
	while(current!=NULL) {
		n++;
		current=current->next;
	}

	printf("List has %d members\n",n);

	return n;

}


bool osb_bb_frac::check_region_frac_size(region_frac * head)
{

	region_frac *current=head;

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

double osb_bb_frac::calc_p_distance()
{

        double sum=0.0;

        for (int user=0;user<lines;user++) {
                sum += pow((_p_budget[user] - tot_pow(user)),2);
        }

	printf("Current p_distance is %lf\n",sqrt(sum));

        return sqrt(sum);

}


