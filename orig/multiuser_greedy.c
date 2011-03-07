#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include "multiuser_load.h"
#include "multiuser_greedy.h"
#include "psd_vector.h"

//double cost_function(int channel,int line_id,double*p,double *old_p,int *b);

multiuser_greedy::multiuser_greedy()
{

	F = (int **)malloc(sizeof(int *)*DMTCHANNELS);
	cost = (double **)malloc(sizeof(double *)*DMTCHANNELS);
	p = (double **)malloc(sizeof(double *)*DMTCHANNELS);
	b = (int **)malloc(sizeof(int *)*DMTCHANNELS);
	p_used = new double[lines];
	_delta_p = new double**[DMTCHANNELS];
	b_total = new int[lines];
	b_target = new int[lines];
	normalised_cost = new double*[DMTCHANNELS];
	user_frozen = new bool [lines];
		
        for (int i=0;i<DMTCHANNELS;i++) {
                F[i] = (int *)malloc(sizeof(int)*lines);
		cost[i] = (double *)malloc(sizeof(double)*lines);
		p[i] = (double *)malloc(sizeof(double)*lines);
		b[i] = (int *)malloc(sizeof(int)*lines);
		_delta_p[i] = new double*[lines];
		normalised_cost[i] = new double[lines];
        }

        for (int i=0;i<DMTCHANNELS;i++) {
                for (int j=0;j<lines;j++) {
                        F[i][j]=0;
			b[i][j]=0;
			p[i][j]=0;
			_delta_p[i][j] = new double[lines];
                }
        }

	this->all_tones_full=false;
	this->p_budget=0.110;
	this->p_error=0.0001;
	this->spectral_mask_on=false;
	this->spectral_mask=dbmhz_to_watts(-36.5);
	this->calc_delta_p_debug=false;
	this->min_cost_debug=false;

	for (int user=0;user<lines;user++) {
		user_frozen[user]=false;
		b_total[user]=0;
		b_target[user]=NOT_SET;	// by default no rate target
		p_used[user]=0;
	}

	psd = new psd_vector;
	psd->_caching_on=false;
}


multiuser_greedy::~multiuser_greedy()
{
	for (int i=0;i<DMTCHANNELS;i++) {
		free(F[i]);
		free(cost[i]);
		free(p[i]);
		free(b[i]);
	}

	free(F);
	free(cost);
	free(p);
	free(b);
	delete p_used;
		
	delete psd;

}


int multiuser_greedy::run()
{

	int tone,user;

	reset_data();
		
	for (tone=0;tone<DMTCHANNELS;tone++) {
		for(user=0;user<lines;user++) {
			calc_delta_p_cost(tone,user);
		}
	}

	//calc_normalised_cost_matrix();
	//print_cost_matrix(0,10);

	//calc_delta_p_debug=true;

	while(1) {
		min_cost(&tone,&user);	
		if(all_tones_full)
			break;
		b[tone][user]++;	// add bit to min cos tone
		b_total[user]++;		
		
		if (b[tone][user] == MAXBITSPERTONE) {
			F[tone][user]=1;
		}

		if (b_total[user] == b_target[user]) {
			freeze_user(user);
		}
		
		update_power(tone,user);	//update power on tone
		
		for (int user=0;user<lines;user++)
			calc_delta_p_cost(tone,user);	// update delta_p on all lines for this tone
	}


	init_lines();
	
	calculate_snr();
	
	set_loading_algo();	

	return 0;
}

void multiuser_greedy::print_cost_matrix(int start, int end)
{

	for (int user=0;user<lines;user++) {
		for (int tone=start;tone<end;tone++) {
			printf("%2.1e ",cost[tone][user]);
		}
		printf("\n");
	}	
}

void multiuser_greedy::calc_normalised_cost_matrix()
{

	double max=0;

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			if (cost[tone][user] > max)
				max=cost[tone][user];
		}
	}

	for (int tone=0;tone<DMTCHANNELS;tone++) {
                for (int user=0;user<lines;user++) {
                        normalised_cost[tone][user]=cost[tone][user]/max;
                }
        }


}

int multiuser_greedy::update_power(int tone, int user)
{

	for (int i=0;i<lines;i++) {
		p[tone][i]+=_delta_p[tone][user][i];
		p_used[i]+=_delta_p[tone][user][i];
	}
	
	return 0;
}

void multiuser_greedy::init_cost_matrix()
{
	for (int tone=0;tone<DMTCHANNELS;tone++) {
                for(int user=0;user<lines;user++) {
                        calc_delta_p_cost(tone,user);
                }
        }

}	


// This function updates the delta_p and _delta_p arrays with the total amount of extra power needed to add an extra bit to a given channel and line
// delta_p stores the total amount for extra power required, _delta_p stores the vector of individual power increases for each line

double multiuser_greedy::calc_delta_p_cost(int channel,int line_id)
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

	//calculate_psd_vector(_b,g,channel,_p);
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

	free(_b);
	free(_p);
	free(_old_p);
	free(g);

	cost[channel][line_id] = new_p-old_p;	// this is the incremental delta_p i.e. total for all lines

	if (calc_delta_p_debug) {
		printf("cost = %6.4g\n",cost[channel][line_id]);
		getchar();
	}

	return 0;


}


void multiuser_greedy::freeze_user(int user)
{
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		F[tone][user]=1;
	}
	user_frozen[user]=true;
}

int multiuser_greedy::min_cost(int *channel,int *line)
{
	int tone,user;
	int min_tone,min_user;
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


	for (user=0;user<lines;user++) {
		if (user_full[user] == true)
			user_frozen[user]=true;
	}


	*channel=min_tone;
	*line=min_user;

	return 0;

}

bool multiuser_greedy::total_power_constraint_broken(int tone,int user)
{
	for (int i=0;i<lines;i++) {
		if (p_used[i] + _delta_p[tone][user][i] > p_budget) {
			//printf("adding a bit on line %d tone %d will break the power constraint on line %d\n",user,tone,i);
			return true;
		}
	}

	return false;
}

bool multiuser_greedy::spectral_mask_constraint_broken(int tone,int user)
{
	if(spectral_mask_on)
		if (p[tone][user] + _delta_p[tone][user][user] > spectral_mask)
			return true;
	return false;
}



bool multiuser_greedy::all_lines_under_budget()
{

	int tone,user;
	double power;

	for (user=0;user<lines;user++) {
		for (tone=0,power=0.0;tone<DMTCHANNELS;tone++) {
			power+=p[tone][user];	
		}
		if (power>p_budget) {
			printf("power is %lf p_used is %lf\n",power,p_used[user]);
			return false;
		}
	}

	return true;

}

double multiuser_greedy::power_used(int user)
{
	double _p=0.0;

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		_p+=p[tone][user];
	}

	return _p;

}

void multiuser_greedy::init_lines()
{

	int tone,user;
	struct line* current;

	for (user=0;user<lines;user++) {
		current=get_line(user);
                for (tone=0;tone<DMTCHANNELS;tone++) {
               		current->b[tone]=b[tone][user];
			current->psd[tone]=watts_to_dbmhz(p[tone][user]);         
                }
        }


}

void multiuser_greedy::set_loading_algo()
{
	int user;
	struct line* current;

	for (user=0;user<lines;user++) {
		current=get_line(user);
		strcpy(current->loading_algo,"multiuser_greedy");
	}

}


