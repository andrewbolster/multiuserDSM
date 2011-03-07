#include "multiuser_weighted.h"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cfloat>

//static double normalised_channel_gain(int line_id,int channel);
//static double get_max_gain(void);

//static double f(double n_delta_xtalk_power,int line_id,int user,int channel);

//static double relative_intraline_goodness(int line_id,int channel);

//static double relative_interline_goodness(int line1,int line2);

static bool DESC(const xtalker_metric &left,const xtalker_metric &right)
{
	return left.lost_bits > right.lost_bits;
}

static double p_tol=0.001;


multiuser_weighted::multiuser_weighted()
{

	flanagans_constant=new double*[DMTCHANNELS];

	_lost_bits = new double**[lines];
	_actual_lost_bits = new double*[lines];
	_xtalker_table = new struct xtalker_metric*[lines];

	w = new double[lines];
	w_min = new double[lines];
	w_max = new double[lines];
	a = new double [lines];
	linear_grad = new double [lines];
	linear_int = new double [lines];
	user_re_init = new bool [lines];

	ave_bits = new double [lines];
	ave_gain = new double [lines];
	
	bits = new double *[lines];

	ave_xtalk = new double *[lines];

	for (int user=0;user<lines;user++) {
		bits[user] = new double [DMTCHANNELS];
		ave_xtalk[user] = new double [lines];
		_lost_bits[user] = new double *[DMTCHANNELS];
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			_lost_bits[user][tone] = new double [lines];
		}
		_actual_lost_bits[user] = new double [lines];
		_xtalker_table[user] = new struct xtalker_metric [lines];
	}

	max_met = new double [lines];
	iter_count=0;
	e=1;

/*
	linear_xt_coeffs = new double **[DMTCHANNELS];
	
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		linear_xt_coeffs[tone] = new double *[lines];
		for (int user=0;user<lines;user++) {
			linear_xt_coeffs[tone][user]=new double [lines];
		}
	}
*/
	

	w[0]=100;
	//a[1]=10.75;
	w[1]=9;
	w[2]=9;

	linear_grad[0]=0.5;
	linear_grad[1]=0.5;

	linear_int[0]=200*(linear_grad[0]);
	linear_int[1]=200*(linear_grad[1]);


	sprintf(cf_string0,"1 + exp(x + %4.2lf)",a[0]);
	sprintf(cf_string1,"1 + exp(x + %4.2lf)",a[1]);

	psd->_caching_on=true;

	for (int user=0;user<lines;user++) {
		w[user] = 10000;
		w_min[user] = NOT_SET;
		w_max[user] = NOT_SET;
		user_re_init[user]=false;
		max_met[user]=-1;	
	}
	
	this->calc_delta_p_debug=false;
	this->cost_function_debug=false;
	this->min_cost_debug=false;
	this->total_bits=0;
	//this->graph_loading=true;
	this->graph_loading=false;

	this->rate_targets_off=false;
	this->w_updates_off=false;

	for (int tone=0;tone<DMTCHANNELS;tone++) {
			this->flanagans_constant[tone] = new double[lines];
			this->flanagans_constant[tone][0]=1.25e8;	// 1.25e8 works well for original near-far
			this->flanagans_constant[tone][1]=1.25e8;
	}
	
	//w[0]=100;
	//w[1]=100;
	//w[2]=100;

}

multiuser_weighted::~multiuser_weighted()
{

//	delete psd;

}



int multiuser_weighted::run_a()
{

        int tone,user,t_line,runs=0;
	FILE* fp[lines];
	bool last_run=false;
	
        reset_data();
	
	calc_ave_bits();

	calc_ave_xtalk();

	for (int user=0;user<lines;user++) {
		double max=0;
		int line=0;
		for (int xtalker=0;xtalker<lines;xtalker++) {
			if (xtalker==user)
				continue;
			if (ave_xtalk[user][xtalker] > max) {
				max = ave_xtalk[user][xtalker];
				line = xtalker;
			}
		}
		printf("Worst xtalker into line %d = line %d\n",user,line);
	}

	//exit(0);

	for (int user=0;user<lines;user++) {
		printf("rate target on line %d is %d\n",user,b_target[user]);
	}

	do {
		iter_count++;
		reset_data();
	
		init_cost_matrix();

		for (int user=0;user<lines;user++) {
			printf("w[%d] = %4.2lf\n",user,w[user]);
		}
	

		while(1) {
			min_cost(&tone,&user);
		
			if(all_tones_full) {
				printf("All tones are full on iteration number %d\n",iter_count);
				break;
			}
			

			b[tone][user]++;        // add bit to min cos tone
			b_total[user]++;
			total_bits++;

			update_power(tone,user);        //update power on tone

			if (graph_loading)
				write_current_stats(tone,user);

			if (b[tone][user] == MAXBITSPERTONE) {
				F[tone][user]=1;
			}

			if (!rate_targets_off) {
				if (b_total[user] == b_target[user]) {
					freeze_user(user);
					printf("just froze user %d cos he got his target\n",user);
					printf("Power used by user %d when frozen = %lf\n",user,p_used[user]);
					init_cost_matrix();
					user_re_init[user]=true;
				}		
			}	

			/*for (int user=0;user<lines;user++) {
				if (user_frozen[user] == true && user_re_init[user] == false) {
					init_cost_matrix();	
					user_re_init[user]=true;
				}	
		
			}*/

			for (int user=0;user<lines;user++)
				calc_delta_p_cost(tone,user);   // update delta_p on all lines for this tone
		}

		for (int user=0;user<lines;user++) {
			printf("rate on line %d is %d\n",user,b_total[user]);
		}

		//getchar();


		//calc_lost_bits();

		//sort_lost_bits();

		//print_lost_bits();

		if (1) {
			char tag[100];
			sprintf(tag,"iteration%d",iter_count);
			init_lines();
			calculate_snr();
			for (int user=0;user<lines;user++) {
				write_mw_stats(user,tag);
			}	
		}


		update_w();	

		//getchar();
	
	} while(!finished());

	if (all_rates_good()) {
		printf("All lines met their rate target\n");
	}
	else {
		double sum=0;
		for (int user=0;user<lines;user++) {
			sum+=pow(b_total[user]-b_target[user],2);
		}
		printf("Total distance from target = %lf\n",sqrt(sum));
	}
	
        init_lines();

        calculate_snr();

        set_loading_algo();

	return 0;
}


bool multiuser_weighted::finished()
{
	//static int b_total_last;
	bool ret;	

	if (w_updates_off)
		return true;

	if (all_rates_good())
		ret = true;
	else
		ret = false;

	//b_total_last = b_total[0];

	return ret;

}


bool multiuser_weighted::all_rates_good()
{

	bool ret=true;

	for (int user=0;user<lines;user++) {
		if (b_target[user] == NOT_SET) {
			//ret = false;
			continue;
		}
		if (b_target[user] - b_total[user] > e)
			return false;
	}
	
/*
	if (ret == false) {
		printf("All rates with targets are good\n");
		printf("Setting all other targets to current rate\n");
		for (int user=0;user<lines;user++) {
			if (b_target[user] == NOT_SET) {
				b_target[user] = b_total[user]+80;
			}
		}
	}
*/
	return ret;

}

void multiuser_weighted::update_w()
{

	double s = 10000;
	bool stop=false;

	if (w_updates_off)
		return;
/*
	for (int user=0;user<lines;user++) {
		if (b_target[user] == NOT_SET) {
			w[user] = w[user];
			continue;
		}
		if (b_total[user] != b_target[user]) {
			double update = s*(b_total[user] - b_target[user]);
			if (w[user] + update < 0)
				w[user] = 0.9*w[user];
			else
				w[user] = w[user] + update;
		}
		else
			w[user] = w[user];
	}
*/

/*
	for (int user=0;user<lines;user++) {
		int bit_diff = b_target[user] - b_total[user];
		double percent_bit_diff = (double)bit_diff/(double)b_target[user];
		if (b_total[user] < b_target[user]) {
			printf("Line %d has not met its rate target\n",user);
			printf("b_target[%d] - b_total[%d] = %d = %lf %\n",user,user,bit_diff,percent_bit_diff);
			if (w[user] > 0) {
				printf("It already has a positive weight so its getting reduced\n");
				printf("Old w[%d] = %lf\n",user,w[user]);
				w[user] = w[user]-w[user]*percent_bit_diff;
				printf("New w[%d] = %lf\n",user,w[user]);
				//w[user] = 0;
				//continue;
			}
			for (int xtalker=0;xtalker<lines;xtalker++) {
			//for (int xtalker=0;xtalker<3;xtalker++) {
				int xt = _xtalker_table[user][xtalker].line_id;
				printf("Number %d xtalker into this line is %d\n",xtalker,xt);
				if (b_total[xt] >= b_target[xt] && xt != user) {
					printf("This line has met its rate target so we will increase its weight\n");
					//w[xt] = w[xt] + 0.1*(b_target[user] - b_total[user])*_xtalker_table[user][xtalker].lost_bits;
					//w[xt] = w[xt] + 0.1*_xtalker_table[user][xtalker].lost_bits;
					printf("This line causes %lf lost bits on line %d\n",_xtalker_table[user][xtalker].lost_bits,user);
					printf("Old w[%d] = %lf\n",xt,w[xt]);
					//w[xt] = w[xt] + s*(b_target[user]-b_total[user])*_xtalker_table[user][xtalker].lost_bits;
					w[xt] = w[xt] + s*percent_bit_diff*_xtalker_table[user][xtalker].percentage;
					printf("New w[%d] = %lf\n",xt,w[xt]);
				}
				else if (b_total[xt] < b_target[xt] && xt != user) {
					printf("This xtalker has not met its rate target so we're going to leave it alone\n");
				}
				if (stop)
					getchar();
			}	
		}
		else { 
			
			//printf("Line %d has met its rate target, so we will adjust its own weight\n",user);
			//printf("Old w[%d] = %lf\n",user,w[user]);
			//w[user] = w[user] - 10*w[user]*percent_bit_diff;
			//printf("New w[%d] = %lf\n",user,w[user]);
			
		}

		if (stop)
			getchar();
	}
	*/

	for (int user=0;user<lines;user++) {
		if (b_total[user] == b_target[user]) {
			if (w[user] == 0)
				w[user]=1;
			else
				w[user]*=1.5;
		}
		else {
			if (w[user] > 0)
				w[user]/=1.5;
		}
	}




}


int multiuser_weighted::update_power(int tone, int user)
{

	for (int i=0;i<lines;i++) {
		p[tone][i]+=_delta_p[tone][user][i];
		p_used[i]+=_delta_p[tone][user][i];
		//_actual_lost_bits[user][i]+=_lost_bits[user][tone][i];
	}

	return 0;
}

int multiuser_weighted::calc_lost_bits()
{	
	struct line* current;
	double bits;
	double *b_vec = new double [lines];
	double *p_vec = new double [lines];

	for (int user=0;user<lines;user++) {
		current=get_line(user);
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			//double r_user;
			double noise=minus_140_dbmhz_watts+current->alien_xtalk[tone];
			memset(b_vec,0,sizeof(double)*lines);	
			memset(p_vec,0,sizeof(double)*lines);
			bits=log2(1+(minus_40_dbmhz_watts*get_channel_gain(user,tone))/(noise*UNdB(12.95)));
			p_vec[user] = minus_40_dbmhz_watts;
			//r_user=relative_intraline_goodness(user,tone);

			for (int xtalker=0;xtalker<lines;xtalker++) {
				if (xtalker!=user) {
					b_vec[xtalker] = (double)b[tone][xtalker];
					//printf("Current user is %d current xtalker is %d\n",user,xtalker);
					//printf("Bits on user %d with no FEXT = %lf\n",user,bits);
					//print_vector(b_vec,"b_vector_before");
					//print_vector(p_vec,"p_vector_before");
					calculate_b_vector_one_bit_set(p_vec,current->gamma[0],tone,xtalker,b_vec);
					//_actual_lost_bits[xtalker][user]+=(bits-b_vec[user])*r_user;
					_actual_lost_bits[xtalker][user]+=(bits-b_vec[user]);
					//print_vector(b_vec,"b_vector_after");
					//print_vector(p_vec,"p_vector_after");
					//printf("bits on user %d with xtalker %d bits at %lf = %lf\n",user,xtalker,b_vec[xtalker],b_vec[user]);
					//printf("bits lost on line %d from line %d on tone %d if %lf\n",user,xtalker,tone,bits-b_vec[user]);
					b_vec[xtalker]=0;
					p_vec[xtalker]=0;
					b_vec[user]=0;
				
				}
			}
			//getchar();			
		}
	}

	delete[] b_vec;
	delete[] p_vec;
	
	return 0;

}


/*void calc_xtalk()
{

	




}*/

double multiuser_weighted::calc_delta_p_cost(int channel,int line_id)
{

	double old_p=0.0,new_p=0.0;
	double *_old_p;
	double *_p,*g;
	int *_b;
	struct line* current;	

	_b=(int *)malloc(sizeof(int)*lines);
	_p=(double *)malloc(sizeof(double)*lines);
	_old_p=(double *)malloc(sizeof(double)*lines);
	g=(double *)malloc(sizeof(double)*lines);


	if(min_cost_debug)
		cost_function_debug=true;

	if(calc_delta_p_debug)
		printf("Calculating delta_p if we add a bit to line %d tone %d\n",line_id,channel);

	for (int user=0;user<lines;user++) {
		current=get_line(user);
		g[user]=current->gamma[current->service[channel]];
		old_p+=p[channel][user];		// current power on sub-channel
		_old_p[user]=p[channel][user];		// copy current power vector on channel to local vector
		_b[user]=b[channel][user];		// copy current b vector to local vector	
	}

	if (calc_delta_p_debug) {
		printf("current p vector = [");
		for (int user=0;user<lines;user++) {
			//printf("%4.2lf ",watts_to_dbmhz(p[channel][user]));
			printf("%6.4g ",p[channel][user]);
		}
		printf("]\n");
		
		printf("current b vector = [");
		for (int user=0;user<lines;user++) {
			printf("%d ",b[channel][user]);
		}
		printf("]\n");
	}

	_b[line_id]++;			// add one bit to local b vector to calculate new p

	if (calc_delta_p_debug) {
		printf("new b vector = [");
		for (int user=0;user<lines;user++) {
			printf("%d ",_b[user]);
		}
		printf("]\n");
	}

	psd->calc(_b,g,channel,_p);

	if (calc_delta_p_debug) {
		printf("new p vector = [");
		for (int user=0;user<lines;user++) {
			//printf("%4.2lf ",_p[user]);
			printf("%6.4g ",_p[user]);
		}
		printf("]\n");
	}

	for (int i=0;i<lines;i++) {
		if (_p[i] < 0) {
		/*	
			printf("a value of p was negative\n");
			printf("b vector is [");
			for (int i=0;i<lines;i++) {
				printf("%d ",_b[i]);
			}
			printf("]\n");

			printf("p vector is [");
			for (int i=0;i<lines;i++) {
                                printf("%4.2lf ",watts_to_dbmhz(_p[i]));
                        }
                        printf("]\n");
			getchar();
		*/
			F[channel][line_id]=1;			
			return 0;
		}
		_delta_p[channel][line_id][i] = _p[i]-_old_p[i];
		new_p+=_p[i];
	}

	if (calc_delta_p_debug) {
		printf("delta_p vector = [");
		for (int user=0;user<lines;user++) {
			printf("%6.4g ",_delta_p[channel][line_id][user]);
		}	
		printf("]\n");
	}

	cost[channel][line_id]=0;	

	/*for (int user=0;user<lines;user++) {
		if (user==line_id) {
			cost[channel][line_id]+=weighting(line_id)*(_p[line_id]-_old_p[line_id]);		// power required to add a bit to this tone
		}
		else {
			//if (_p[user] == 0)
			//	cost[channel][line_id]=_p[line_id]-_old_p[line_id];
			//else
				cost[channel][line_id]+=_b[line_id]*flanagans_constant*_p[line_id]*get_xtalk_gain(line_id,user,channel)*get_xtalk_gain(line_id,line_id,channel)/(get_xtalk_gain(user,user,channel));	
		}
	}*/

	//cost[channel][line_id] = cost_function1(channel,line_id,_p,_old_p,_b,g);
	cost[channel][line_id] = last_cf(channel,line_id,_p,_old_p,_b,g);
	//cost[channel][line_id] = linear_cost_function(channel,line_id,_p,_old_p,_b,g);
	//cost[channel][line_id] = cost_function(channel,line_id,_p,_old_p,_b);
	//cost[channel][line_id] = xtalk_only_cost_function(channel,line_id,_p,_old_p,_b);

	if (calc_delta_p_debug) {
		int c;
		printf("cost = %6.4g\n",cost[channel][line_id]);
		c=getchar();
		if (c=='q')
			calc_delta_p_debug=false;
	}


	free(_b);
	free(_p);
	free(_old_p);
	free(g);

	return 0;


}


int multiuser_weighted::min_cost(int *channel,int *line)
{
	int tone,user;
	int min_tone=0,min_user=0;
	double min=1e100;
	bool *user_full = new bool[lines];	

	for (int user=0;user<lines;user++)
		user_full[user] = true;

	all_tones_full=true;		// if i put this variable below this line, it segfaults.... compiler bug or weird bug in my code??

	for (tone=0;tone<DMTCHANNELS;tone++) {
		for (user=0;user<lines;user++) {
			if (F[tone][user]==1) {
				continue;
			}
			user_full[user]=false;		// if we have reached here, the current current is not full!
			if (cost[tone][user] < min) {
				if (total_power_constraint_broken(tone,user) || spectral_mask_constraint_broken(tone,user)) {
					F[tone][user] = 1;
					continue;
				}
				all_tones_full=false;
				min=cost[tone][user];
				min_tone=tone;
				min_user=user;
					
			}
		}
	}
	
	if (all_tones_full) {
		//printf("All tones now frozen out!\n");
		//getchar();
	}

	
	
	if (min_cost_debug) {
		int c;
		double dummy;
		printf("min cost was on channel %d line %d with value %e\n",min_tone,min_user,min);
		printf("new number of bits = %d\n\n",b[min_tone][min_user]+1);
		cost_function_debug=true;
		calc_delta_p_cost(min_tone,min_user);	
		c=getchar();
		if (c=='q')
			min_cost_debug=false;
	}


	//bool re_init=false;
	for (user=0;user<lines;user++) {
		if (user_full[user] == true && user_re_init[user]==false) {	
			//re_init=true;
			user_frozen[user]=true;
			user_re_init[user]=true;
			init_cost_matrix();
			printf("Re-initialising cost matrix as line %d has used all its power budget\n",user);
			printf("Iteration number = %d\n",iter_count);
		}
	}



	*channel=min_tone;
	*line=min_user;

	return 0;

}

bool multiuser_weighted::total_power_constraint_broken(int tone,int user)
{
	for (int i=0;i<lines;i++) {
		if (p_used[i] + _delta_p[tone][user][i] > p_budget+p_tol) {
			//printf("adding a bit on line %d tone %d will break the power constraint on line %d\n",user,tone,i);
			return true;
		}
	}

	return false;
}

bool multiuser_weighted::spectral_mask_constraint_broken(int tone,int user)
{
	if(spectral_mask_on)
		if (p[tone][user] + _delta_p[tone][user][user] > spectral_mask)
			return true;
	return false;
}


double multiuser_weighted::xtalk_only_cost_function(int channel,int line_id,double*p,double *old_p,int *b)
{
	double xtalk_ratio=0;
	double delta_p=0;

	for (int user=0;user<lines;user++) {
		if (user==line_id)
			delta_p=p[line_id]-old_p[line_id];
		else
			xtalk_ratio+=(p[line_id]-old_p[line_id])*get_xtalk_gain(line_id,user,channel)/get_xtalk_gain(user,user,channel);
	}

	return xtalk_ratio;
}


double multiuser_weighted::cost_function(int channel,int line_id,double*p,double *old_p,int *b)
{

	double delta_p=0;
	double ret=0;

	delta_p=p[line_id]-old_p[line_id];

	for (int user=0;user<lines;user++) {
		if (user==line_id)
			ret+=delta_p;
		else
			ret+=b[line_id]*flanagans_constant[channel][user]*p[line_id]*get_xtalk_gain(line_id,user,channel)*get_xtalk_gain(line_id,line_id,channel)/(get_xtalk_gain(user,user,channel));
		/*
		else 
			pen=get_xtalk_gain(line_id,user,channel)/get_xtalk_gain(user,user,channel);
		*/	
	}

	//if (b_target[line_id]==b_total[line_id])
	//	return 1e300;

	//return exp(0.01/(p_budget-p_used[line_id]))*exp(10/(b_target[line_id]-b_total[line_id]))*ret;

//	printf("channel %d\tline %d\tdelta_p = %e\tcost = %e\tratio = %.15lf\t",channel,line_id,delta_p,ret,ret/delta_p);
	/*
	for (int user=0;user<lines;user++) {
		FILE *fp; 
		fp = fopen("test.data","a");
		if (user!=line_id)
			fprintf(fp,"%.15lf %.15lf\n",ret/delta_p,10*log10(delta_p*get_xtalk_gain(line_id,user,channel)/get_channel_gain(user,channel)));
		fclose(fp);
	}
	*/
//      getchar();


	

	return ret;

}

double multiuser_weighted::ya_cost_function(int channel,int line_id,double*p,double *old_p,int *b,double *g)
{

	double delta_p=p[line_id]-old_p[line_id];
	double pen=1;

	for (int user=0;user<lines;user++) {
		if (user==line_id)
			;
		else {
			double xtalk_ratio=10*log10(delta_p*get_channel_gain(line_id,channel)*get_xtalk_gain(line_id,user,channel)/get_channel_gain(user,channel));
			pen += exp(xtalk_ratio+a[line_id]);
		}
	}
	
	return delta_p*(pen);

}


double multiuser_weighted::linear_cost_function(int channel,int line_id,double*p,double *old_p,int *b,double *g)
{

	double delta_p=p[line_id]-old_p[line_id];
	double pen=1;

	for (int user=0;user<lines;user++) {
		if (user==line_id)
			;
		else {
			double xtalk_ratio=10*log10(delta_p*get_xtalk_gain(line_id,user,channel)/get_channel_gain(user,channel));
			if (xtalk_ratio < -200)
				;
			else
				pen += linear_grad[line_id]*xtalk_ratio+linear_int[line_id];
		}
	}
	
	return delta_p*(pen);

}

double multiuser_weighted::cost_function1(int channel,int line_id,double*p,double *old_p,int *b,double *g)
{


	double bits_per_watt=0;
	double total_bits_lost=0;
	double delta_p;
	double xtalk_power;
	double n_xtalk_power;
	double minn_xtalk_power;
	double fn_xtalk_power;
	double n_delta_xtalk_power;
	double fn_delta_xtalk_power;
	double n_gain;
	double fn_gain;
	double dn_gain;
	double xt_gain_ratio;
	double pen=1;
	double _cost;

	delta_p=p[line_id]-old_p[line_id];

	//if (channel==108)
		//cost_function_debug=true;
	
	if (cost_function_debug) {
		printf("delta_p on line %d tone %d is %4.2e W\n",line_id,channel,delta_p);
		//printf("normalised channel gain is %4.2lf\n",normalised_channel_gain(line_id,channel));
	}

	for (int user=0;user<lines;user++) {
		if (user==line_id)
			;
		else {	
			
			if (!user_frozen[user]) {
				//dn_gain=normalised_channel_gain(line_id,channel) - normalised_channel_gain(user,channel);
				//fn_gain=500/(1+exp(-dn_gain*100000+1000));
				//xt_gain_ratio=get_xtalk_gain(line_id,user,channel)/get_channel_gain(user,channel);
				/*if (dn_gain > 0)
					fn_gain=100*dn_gain;
				else
					fn_gain=0;*/
				xtalk_power=p[line_id]*get_xtalk_gain(line_id,user,channel);
				n_xtalk_power=10*log10(xtalk_power/get_channel_gain(user,channel));
				//fn_delta_xtalk_power=linear_grad[line_id]*n_delta_xtalk_power+linear_int[line_id];	
				minn_xtalk_power=10*log10(pow(10,g[line_id]/10)*4.31e-14*get_xtalk_gain(line_id,user,channel)/(get_channel_gain(line_id,channel)*get_channel_gain(user,channel)));
				//fn_delta_xtalk_power=1+exp((n_delta_xtalk_power-minn_delta_xtalk_power)-a[line_id]);	
				// try to figure out the constant in the expression below
				// it should return a high number, i.e. no damping if the other line channel gain is below a threshold
				// it should return a lower number the better the other channel is, maybe the gradient is the adjustable?	
				// its also a function of the crosstalk
				n_delta_xtalk_power=n_xtalk_power-minn_xtalk_power;
				fn_delta_xtalk_power=f(n_delta_xtalk_power,line_id,user,channel);
				fn_xtalk_power=200/(1+exp(-10*(n_delta_xtalk_power)+fn_delta_xtalk_power));	
				//pen+=exp(n_delta_xtalk_power+a[line_id])+fn_gain;
				//pen+=pow(10,fn_gain/10)*xt_gain_ratio+pow(10,fn_delta_xtalk_power/10);
				pen=pow(10,fn_xtalk_power/10);
				if (cost_function_debug) {
					printf("n_xtalk_p = %4.2lf\tminn_xtalk_power = %4.2lf\tn_delta_xt_p = %4.2lf\n",n_xtalk_power,minn_xtalk_power,n_delta_xtalk_power);
					printf("f(n_delta_xtalk_power) = %4.2lf\n",fn_delta_xtalk_power);
					printf("f(fn_xtalk_power) = %4.2lf\n",fn_xtalk_power);
	/*				printf("delta_xtalk_power into line %d is %4.2e W\n",user,delta_xtalk_power);
					printf("normalised delta_xtalk_power into line %d is %4.2lf dB\n",user,n_delta_xtalk_power);
					printf("normalised channel gain of line %d is %.12lf\n",user,normalised_channel_gain(user,channel));
					printf("delta normalised channel gain = %.12lf\n",dn_gain);
					printf("a[%d] = %4.2lf\n",line_id,a[line_id]);
					printf("f(n_delta_xtalk_gain) = %4.2e delta_xt_penalty = %4.2lf dB\n",pow(10,fn_delta_xtalk_power/10),fn_delta_xtalk_power);
					printf("dn_gain penalty = %4.2lf dB\n",fn_gain);
					printf("xtalk/gain ratio = %4.2lf dB\n",10*log10(xt_gain_ratio));*/
					printf("penalty is %4.2e = %4.2lf dB\n",pen,10*log10(pen));
				}
			}
			else
				pen+=0;
			
		}

	}


	//printf("channel %d\tline %d\tdelta_p = %e\tcost = %e\n",channel,line_id,delta_p,_cost);

	//return total_bits_lost*delta_p;
	if (cost_function_debug) {
		printf("Cost is %4.2e\n",delta_p*(pen));
		cost_function_debug=false;
		getchar();
	}
	return delta_p*(pen);


}
/* this should return between zero and one, increasing in cost if a line is near its rate target */

double multiuser_weighted::last_cf(int channel,int line_id,double*p,double *old_p,int *b,double *g)
{

	double pen=0,bpen=0,rpen=0,ypen=0;
	//double delta_p=p[line_id]-old_p[line_id];
	double delta_p=0;
	int nbits=1;

	double *_p = new double [lines];
	double *_oldb = new double [lines];
	double *_b = new double [lines];

	for (int user=0;user<lines;user++) {
		delta_p+=p[user]-old_p[user];
	}	

	if (((channel<20) || (channel>120 && channel<150)) && line_id==1) {
		//cost_function_debug=true;
	}
	
	_oldb[line_id]=this->b[channel][line_id];

	//print_vector(_oldb,"b vector before");

	for (int user=0;user<lines;user++) {			// set p on all other lines to -40
		if (user!=line_id) {
			//_p[user]=dbmhz_to_watts(-36.5);
			_p[user]=minus_40_dbmhz_watts;
		}
		else
			_p[user]=0;
	}

//	print_vector(_p,"old p vector");

	calculate_b_vector_one_bit_set(_p,g,channel,line_id,_oldb);		// _oldb will contain the vector of bits

	bool recalc=false;
	for (int user=0;user<lines;user++) {
		if (_oldb[user] > MAXBITSPERTONE) {
			_oldb[user]=MAXBITSPERTONE;
			recalc=true;
		}
	}

	if (recalc)
		calculate_psd_vector(_oldb,g,channel,_p);		// FIXME	this one is doubles for b vector


//	print_vector(_oldb,"b after setting all other lines to -40");

	_b[line_id]=this->b[channel][line_id]+nbits;			// add one bit to the line we're calculating cost for
	
	calculate_b_vector_one_bit_set(_p,g,channel,line_id,_b);

//	print_vector(_b,"b after adding 10 bit to line id");

/*
	calculate_b_vector_from_psd(_p,g,channel,_oldb);

	_p[line_id]=dbmhz_to_watts(-40);

	calculate_b_vector_from_psd(_p,g,channel,_b);
*/
	for (int user=0;user<lines;user++) {
		if (user==line_id)
			;
		else {
			double lost_bits = _oldb[user]-_b[user];
			_lost_bits[line_id][channel][user]=lost_bits;
			if (!user_frozen[user]) {
				double r_user,r_line_id,r_inter,r_inter_channel;
				double this_line_lost_bits;	// bits lost when line user has psd of -40
				double this_line_bits;		// bits with zero crosstalk on line $line_id
				//double *b1=new double [lines];
				if (cost_function_debug) {
					printf("delta_p on line %d channel %d is %4.2e\n",line_id,channel,delta_p);
					printf("bits lost on line %d by loading %d bits on line %d\n",user,nbits,line_id);
					printf("bits lost = %4.2lf\n",lost_bits);
				}
				r_user=relative_intraline_goodness(user,channel);
				//r_line_id=relative_intraline_goodness(line_id,channel);
				//r_inter=relative_interline_goodness(line_id,user);
				//r_inter_channel=relative_interline_goodness_channel(line_id,user,channel);

				bpen+=w[line_id]*lost_bits*r_user;
				/*
				rpen+=bpen*1/(1+exp(-3000*(r_line_id-r_user)+3100));
				//rpen+=bpen*1/(1+exp(-3000*(r_line_id-r_user)+800));
			
				for(int user1=0;user1<lines;user1++) {
					_p[user1]=0;
				}	
	
				this_line_bits=log2(1+minus_40_dbmhz_watts*get_channel_gain(line_id,channel)/(UNdB(g[line_id])*minus_140_dbmhz_watts));
				_p[line_id]=minus_40_dbmhz_watts;
				_p[user]=minus_40_dbmhz_watts;

				calculate_b_vector_from_psd(_p,g,channel,b1);

				this_line_lost_bits=this_line_bits-b1[line_id];
				
				ypen+=15/(1+exp(3000*(r_line_id-r_user)+3100))*this_line_lost_bits;
				//ypen+=15/(1+exp(3000*(r_line_id-r_user)+800))*this_line_lost_bits;
				*/
				pen+=bpen;
				if (pow(10,pen/10) == DBL_MAX) {
					printf("Ouch!\n");
					getchar();
				}
				/*
				if (cost_function_debug) {
					printf("penalty due to bits lost = %4.2lf dB\n",bpen);
					printf("r_line_id-r_user = %4.2lf\n",r_line_id-r_user);
					printf("penalty reduction = %4.2lf dB\n",rpen);
					printf("bits lost on line %d due to -40 dbmhz on line %d = %4.2lf\n",line_id,user,this_line_lost_bits);
					printf("yielding penalty is %4.2lf dB\n",ypen);
					printf("total penalty is %4.2lf dB\n",pen);
					printf("relative intra-line goodness of line_id %d is %4.2lf\n",line_id,
									r_line_id);
					printf("relative intra-line goodness of user %d is %4.2lf\n",user,
									r_user);
					//printf("relative inter-line goodness is %4.2lf\n",
					//				r_inter);
					//printf("relative inter-line goodness on channel %d is %4.2lf\n",channel,
					//				r_inter_channel);


					getchar();
				}
				*/
				//delete[] b1;
			}
			else {
				pen+=0;
			}
		}
	}

//	getchar();

	//cost_function_debug=false;
	//
	delete[] _p;
	delete[] _oldb;
	delete[] _b;

	return delta_p*pow(10,pen/10);


}


void multiuser_weighted::sort_lost_bits()
{

 	for (int user=0;user<lines;user++) {
                double lost_bits_total=0;
                
		for (int xtalker=0;xtalker<lines;xtalker++) {
			_xtalker_table[user][xtalker].line_id=xtalker;
                        _xtalker_table[user][xtalker].lost_bits=_actual_lost_bits[xtalker][user];
			lost_bits_total+=_actual_lost_bits[xtalker][user];
                }
	
                std::sort(_xtalker_table[user],_xtalker_table[user]+lines,DESC);
		
		for (int xtalker=0;xtalker<lines;xtalker++) {
			_xtalker_table[user][xtalker].percentage=_xtalker_table[user][xtalker].lost_bits/lost_bits_total;			
		}

        }


}


void multiuser_weighted::print_lost_bits()
{
	for (int user=0;user<lines;user++) {
		for (int xtalker=0;xtalker<lines;xtalker++) {
			printf("Bits lost on line %d caused by line %d = %lf percentage = %lf\n",user
												,_xtalker_table[user][xtalker].line_id
												,_xtalker_table[user][xtalker].lost_bits
												,_xtalker_table[user][xtalker].percentage);
		}
		printf("\n");
	}


}

double multiuser_weighted::f(double n_delta_xtalk_power,int line_id,int user,int channel)
{
	double bit_threshold=pow(10,-90/10);
	double r1,r2,r3,ret,met;
	
	if (get_channel_gain(user,channel) < bit_threshold) {
		printf("no damping from line %d into line %d on channel %d because user %d cant load any bits here\n",line_id,user,channel,user);
		getchar();
		return 1000;		// no damping for line_id as user cant use this channel anyway
	}
	r1=relative_intraline_goodness(user,channel);		// as this goes up, damping should increase, i.e return value down
	r2=relative_intraline_goodness(line_id,channel);	// as this goes up, damping should decrease
	
	r3=relative_interline_goodness(line_id,user);		// if line_id is way better than user then damping should increase, return value down
	
	
	
	met=get_xtalk_gain(line_id,user,channel)*(r2*r3/r1);

	if (met > max_met[line_id])
		max_met[line_id] = met;

	printf("f for line %d into line %d on channel %d was %4.2e = %4.2lf\n",line_id,user,channel,met,10*log10(met));
	printf("max met on this line is %4.2lf\n",10*log10(max_met[line_id]));

	ret = ((10*log10(met)-10*log10(max_met[line_id]))*-1*a[line_id]); 

	printf("going to return %lf\n",ret);
	getchar();
	
	return ret;

}

double multiuser_weighted::relative_intraline_goodness(int line_id,int channel)
{
	
	return bits[line_id][channel]/ave_bits[line_id];
	//return MAX(15,bits[line_id][channel])/ave_bits[line_id];
	//return get_channel_gain(line_id,channel)/ave_gain[line_id];
}

double multiuser_weighted::relative_interline_goodness(int line1,int line2)
{

	return ave_bits[line1]/ave_bits[line2];

}

double multiuser_weighted::relative_interline_goodness_channel(int line1,int line2,int channel)
{

	return bits[line1][channel]/bits[line2][channel];

}


void multiuser_weighted::calc_ave_bits()
{

	for (int user=0;user<lines;user++) {
		ave_bits[user]=0;
		ave_gain[user]=0;
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			bits[user][tone] = log2(1+(get_channel_gain(user,tone)*minus_40_dbmhz_watts)/(pow(10,9.95/10)*minus_140_dbmhz_watts));
			ave_bits[user]+=bits[user][tone];
			ave_gain[user]+=get_channel_gain(user,tone);
		}
		ave_bits[user]/=DMTCHANNELS;
		ave_gain[user]/=DMTCHANNELS;
	}

}

void multiuser_weighted::calc_ave_xtalk()
{

	for (int user=0;user<lines;user++) {
		for (int xtalker=0;xtalker<lines;xtalker++) {
			ave_xtalk[user][xtalker]=0;
			if (xtalker==user) {
				;
			}
			else {	
				for (int tone=0;tone<DMTCHANNELS;tone++) {
					ave_xtalk[user][xtalker]+=get_xtalk_gain(xtalker,user,tone)/get_channel_gain(xtalker,tone);	
				}
			}
							
		}
	}
}
/*
static double normalised_channel_gain(int line_id,int channel)
{

	double gain = get_channel_gain(line_id,channel);
	double max_gain = get_max_gain();

	return gain/max_gain;

}
*/
/*
static double get_max_gain()
{
	double gain;
	double max=DBL_MIN;
	for (int user=0;user<lines;user++) {
		gain=get_channel_gain(user,0);
		if (gain>max)
			max=gain;
	}

	return max;	
}
*/


double multiuser_weighted::weighting(int line_id)
{
	double factor;
	double b_tar=(double)b_target[line_id];
	double b_tot=(double)b_total[line_id];

	return 1;

	factor = (b_tar-b_tot)/b_tar;

	if (b_target[line_id] == DMTCHANNELS*MAXBITSPERTONE)
		return 1;
	else {
		//printf("return %lf\n",1-factor);
		return 1-factor;
	}
}

void multiuser_weighted::set_loading_algo()
{
	int user;
	struct line* current;

	for (user=0;user<lines;user++) {
		current=get_line(user);
		strcpy(current->loading_algo,"multiuser_weighted");
	}

}

void multiuser_weighted::write_current_stats(int tone,int user)
{

	char fn[50];
	FILE *fp;

	sprintf(fn,"%s/data/mw_test/stats/%d.txt",ROOT_DIR,total_bits);

	fp = fopen(fn,"w");

	if (fp==NULL)
		exit(2);
	
	fprintf(fp,"#%s\n",cf_string0);
	fprintf(fp,"#%s\n",cf_string1);
	fprintf(fp,"#bit added was on tone %d line %d\n",tone,user);

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		fprintf(fp,"%d\t",tone);
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%d\t",b[tone][user]);
		}
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%4.2lf\t",watts_to_dbmhz(p[tone][user]));
		}
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%4.2e\t",cost[tone][user]);
		}
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%d\t",F[tone][user]);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

}


void multiuser_weighted::write_mw_stats(int line_id,char *tag)
{

	int i;
	struct line* c;
	char fn[255];
	FILE *fp;

	c = get_line(line_id);

	if (tag != NULL)
		sprintf(fn,"%s/data/mw_loading/line%d_%s_%s.txt",ROOT_DIR,line_id,scen,tag);
	else
		sprintf(fn,"%s/data/mw_loading/line%d_%s.txt",ROOT_DIR,line_id,scen);		

	if ((fp = fopen(fn,"w")) == NULL) {
		printf("cannot open file to write\n");
		exit(1);
	}

	if (c->is_dual == 1) {
		printf("warning this is a dual qos line but you are printing stats for single qos\n");
		fprintf(fp,"#warning, this line is dual qos\n");
	}

	fprintf(fp,"#line_id=%d\n",line_id);
	fprintf(fp,"#length = %4.2lf lt = %4.2lf nt = %4.2lf\n",c->length,c->lt,c->nt);
	fprintf(fp,"#p_total = %lf mW\n",c->p_total*1e3);
	fprintf(fp,"#b_total = %d  data_rate = %d approx\n",c->b_total,c->b_total*FRAMERATE);
	fprintf(fp,"#background noise = %4.2lf dBm/Hz\n",c->background_noise);
	fprintf(fp,"#w[%d] = %6.4lf\n",line_id,w[line_id]);
	fprintf(fp,"#tone\tbits\tpsd\t\tsnr\t\t\tcnr\t\t\tgain\t\t\tmargin\tp_e\t\t");

	/*for (int line=0;line<lines;line++) {
		fprintf(fp,"xt_n_%d\t",line);
	}*/
	fprintf(fp,"\n");

	for (i=0;i<DMTCHANNELS;i++) {
		fprintf(fp,"%d\t%d\t%4.2lf\t\t%4.2e  %4.2lf dB\t%4.2e  %4.2lf dB\t%4.2e  %5.3lf dB\t%3.1lf\t%4.2e\t",
									i,c->b[i],c->psd[i],c->snr[i],10*log10(c->snr[i]),
									c->cnr[i],10*log10(c->cnr[i]),c->gain[i],10*log10(c->gain[i]),
									c->gamma_m[i],c->symerr[i]);
		/*for (int line=0;line<lines;line++) {
			fprintf(fp,"%4.2lf\t",10*log10(c->xtalk_noise[line][i]/get_channel_gain(line_id,i)));
		}*/
		fprintf(fp,"\n");
	}
	
	fclose(fp);

}

int multiuser_weighted::run()
{

        int tone,user;
	int t_line=0;
	//int runs=0;
	//FILE* fp[lines];
	bool last_run=false;
	
        reset_data();


	//fp[0]=fopen("line0_order.txt","w");
	//fp[1]=fopen("line1_order.txt","w");

        

        //calc_normalised_cost_matrix();
        //print_cost_matrix(0,10);

        //calc_delta_p_debug=true;

	for (int user=0;user<lines;user++) {
		printf("rate target on line %d is %d\n",user,b_target[user]);
		if (b_target[user]!=NOT_SET) {
			t_line=user;
		}
	}


	while((abs(b_target[t_line] - b_total[t_line]) > 5) || last_run) {
		reset_data();
		
		for (tone=0;tone<DMTCHANNELS;tone++) {
                	for(user=0;user<lines;user++) {
                        	calc_delta_p_cost(tone,user);
	                }
	        }

		//print_cost_matrix(0,DMTCHANNELS-1);

		while(1) {
			min_cost(&tone,&user);
			if(all_tones_full) {
				break;
			}
			
	
			b[tone][user]++;        // add bit to min cos tone
			b_total[user]++;
			total_bits++;

			update_power(tone,user);        //update power on tone
			//fprintf(fp[user],"%d %d\n",tone,total_bits);

			if (last_run && graph_loading)
				write_current_stats(tone,user);

			if (b[tone][user] == MAXBITSPERTONE) {
				F[tone][user]=1;
			}

			if (b_total[user] == b_target[user]) {
				freeze_user(user);
				if (last_run)
					printf("just froze user %d cos he got his target\n",user);
			}			

			for (int user=0;user<lines;user++)
				calc_delta_p_cost(tone,user);   // update delta_p on all lines for this tone
		}
		
		/*for (int user=0;user<lines;user++) {
			printf("b_%d = %d\n",user,b_total[user]);
		}*/

		//printf("w = %lf\n",w[t_line]);
		//getchar();

		if (last_run) {
			printf("this is the last run so setting last_run to false\n");
			last_run=false;
			break;
		}
	
		if (abs(b_target[t_line] - b_total[t_line]) < 5) {
			printf("Met the rate target, setting the last_run to true\n");
			last_run=true;
			continue;
		}
	
		if (b_total[t_line] > b_target[t_line]) {	// if rate > target_rate then this is the min weight
		//	printf("Rate is greater than target, so increase the weight and we have found the min weight\n");
			w_min[t_line]=w[t_line];		// set the min weight
			if (w_max[t_line] < 0) { 			// dont yet know the max weight
		//		printf("dont yet know the max weight so increase the weight\n");
				w[t_line]*=2;
			}
		}
		else {						// this is the max weight
		//	printf("Rate is less than target so decrease the weight and we have found the max weight\n");
			w_max[t_line]=w[t_line];
			if (w_min[t_line] < 0) {
		//		printf("dont yet know the min weight so decrease the weight\n");
				w[t_line]/=2;
			}
		}

		if (w_max[t_line] > 0 && w_min[t_line] > 0) {
			w[t_line] = (w_max[t_line]+w_min[t_line])/2;	
		//	printf("Currently bisecting, new weight is %lf\n",w[t_line]);
		}


	}

	//for (int user=0;user<lines;user++)
	//	fclose(fp[user]);

        init_lines();

        calculate_snr();

        set_loading_algo();
	
	return 0;

}
