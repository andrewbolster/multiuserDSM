#include "mipb_frac.h"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cfloat>
#include <cassert>
#include "multiuser_load.h"

mipb_frac::mipb_frac()
{

	_delta_p = new double**[DMTCHANNELS];
	_cost = new double*[DMTCHANNELS];
	_b = new double*[DMTCHANNELS];
	_p = new double*[DMTCHANNELS];
	_F = new int*[DMTCHANNELS];

	_wp = new double[lines];
	_w = new double[lines];
	_w_min = new double[lines];
	_w_max = new double[lines];
	_wxt = new double[lines];
	_p_used=new double[lines];
	_p_budget=new double[lines];
	_b_total=new double[lines];
	_rate_targets = new int[lines];

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		_delta_p[tone] = new double*[lines];
		_cost[tone] = new double[lines];
		_b[tone] = new double[lines];
		_p[tone] = new double[lines];
		_F[tone] = new int[lines];
		for (int user=0;user<lines;user++) {
			_b[tone][user]=0;
			_p[tone][user]=0;
			_F[tone][user]=0;
			_cost[tone][user]=0;
		}
	}

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			_delta_p[tone][user] = new double[lines];
			for (int user1=0;user1<lines;user1++) {
				_delta_p[tone][user][user1] = 0;
			}
		}
		
	}


	for (int user=0;user<lines;user++) {
		_wp[user]=1;
		_w[user]=0;
		_wxt[user]=1;
		_p_budget[user]=0.110;
		_p_used[user]=0;
		_rate_targets[user]=NOT_SET;
	}

	//_psd = new psd_vector;
	
	_p_ave=0.0;
	_last_bit_cost=0.0;

	_spectral_mask_on=true;
	_spectral_mask=dbmhz_to_watts(-30);

	_graph_loading=false;
	_greedy=false;

	_rate_tol=5;

	_bit_inc=1;
	_wp_graph = fopen("../data/mipb_wp_graph_most_recent.txt","w");
	if (_wp_graph == NULL) {
		exit(2);
	}
}


void mipb_frac::reset_data()
{

        for (int tone=0;tone<DMTCHANNELS;tone++) {
                for (int user=0;user<lines;user++) {
                        _delta_p[tone][user] = new double[lines];
                        for (int user1=0;user1<lines;user1++) {
                                _delta_p[tone][user][user1]=0;
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
        _p_ave=0;
}



mipb_frac::~mipb_frac()
{

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			delete[] _delta_p[tone][user];
		}
		delete[] _delta_p[tone];
	}
	delete[] _delta_p;


	//delete _psd;
	fclose(_wp_graph);

}


int mipb_frac::run()
{
	bool targets=false;
	for (int user=0;user<lines;user++) {
		if (_rate_targets[user] != NOT_SET)
			targets=true;
	}	
	
	
	if (targets) {
		for (int user=0;user<lines;user++) {
			if (_rate_targets[user] == NOT_SET) 
				continue;
			_w_min[user]=0;
			_w_max[user]=1e-5;
			while(1) {
				reset_data();
				load();
				
				print_vector(_w,"_w");
                                print_vector(_b_total,"_b_total");
				//getchar();

				if (_b_total[user] > _rate_targets[user]) {
					_w_max[user] = _w[user];
				}
				else if (_b_total[user] < _rate_targets[user]) {
					_w_min[user] = _w[user];
				}

				if (abs(_b_total[user] - _rate_targets[user]) < _rate_tol) {
					printf("Converged on line %d\n",user);
					break;
				}
				_w[user] = (_w_max[user]+_w_min[user])/2;	
			}
		}
	}
	else
		load();

	init_lines();
	
	calculate_snr();

	
	

}


void mipb_frac::load()
{

	int min_tone,min_user;
	FILE *lead_lag_fp;
	FILE *p_fp;

	lead_lag_fp=fopen("../data/mipb_lag_lead.txt","w");
	p_fp=fopen("../data/mipb_p_against_b.txt","w");

	if (p_fp == NULL || lead_lag_fp == NULL) {
		printf("Boom!\n");
		exit(2);
	}	
	
	init_cost_matrix();
	
	while(1) {
		min_cost(&min_tone,&min_user);

		//printf("Min tone/user was %d/%d\n",min_tone,min_user);
		//getchar();

		//if (min_tone == 115 && min_user == 1)
		//	getchar();


		if (_all_tones_full) {
			printf("All tones full\n");
			if (_graph_loading) {
				write_current_stats(min_tone,min_user);
			}	
			//for (int tone=0;tone<DMTCHANNELS;tone++) {
			//	print_vector(_F[tone],"F");
			//}
			//for (int tone=0;tone<DMTCHANNELS;tone++) {
			//	print_vector(_cost[tone],"cost");
			//}
			//print_vector(_wp,"wp");
			
			break;
		}
		_b[min_tone][min_user]+=_bit_inc;;        // add bit to min cos tone
                _b_total[min_user]+=_bit_inc;
                _total_bits+=_bit_inc;

		update_power(min_tone,min_user);

		if (_b[min_tone][min_user] >= MAXBITSPERTONE) {
			_F[min_tone][min_user] = 1;
		}

		//print_vector(_p[min_tone],"p on tone");
                //getchar();
		//print_vector(_b_total,"b_total");
		//printf("total bits = %d\n",_total_bits);
		//print_vector(_p_used,"p_used");

		if (!_greedy)
			update_wp();

		//print_vector(_wp,"wp");
                //getchar();

		fprintf(lead_lag_fp,"%lf\t",_total_bits);
		fprintf(p_fp,"%lf\t",_total_bits);
		for (int user=0;user<lines;user++) {
			fprintf(lead_lag_fp,"%6.4g\t",_p_used[user]-_p_ave);
			fprintf(p_fp,"%6.4g\t",_p_used[user]);
		}
		fprintf(lead_lag_fp,"\n");
		fprintf(p_fp,"\n");
		
		

		
		if ((int)_total_bits % 5 == 0) {
			if (_graph_loading) {
				write_current_stats(min_tone,min_user);
			}
		}
		

		for (int user=0;user<lines;user++) {
			calc_delta_p(min_tone,user);
		}
	
		
		for (int tone=0;tone<DMTCHANNELS;tone++) {
                	for (int user=0;user<lines;user++) {
		                _cost[tone][user] = cf(tone,user);
		                //_cost[tone][user] = var_cf(tone,user);
                	}
                }
	
		dump_cost_matrix();

		//print_vector(_cost[min_tone],"New costs");
		//getchar();

	}

	fclose(p_fp);
	fclose(lead_lag_fp);

	
	if (_graph_loading) {
		write_stats_totals();
	}

}



void mipb_frac::init_cost_matrix()
{
        for (int tone=0;tone<DMTCHANNELS;tone++) {
                for(int user=0;user<lines;user++) {
                        calc_delta_p(tone,user);
                        _cost[tone][user] = cf(tone,user);
                }
        }
	dump_cost_matrix();
}


int mipb_frac::min_cost(int *channel,int *line)
{

	int tone,user;
	int min_tone,min_user;
	double min=DBL_MAX;
	bool *user_full = new bool[lines];
	double cost=0.0;

	for (int user=0;user<lines;user++)
		user_full[user] = true;

	_all_tones_full=true;		// if i put this variable below this line, it segfaults.... compiler bug or weird bug in my code??

	for (tone=0;tone<DMTCHANNELS;tone++) {
		for (user=0;user<lines;user++) {
			if (_F[tone][user]==1) {
				continue;
			}
			user_full[user]=false;		// if we have reached here, the current current is not full!
			//assert(_cost[tone][user] > 0);
			//assert(_cost[tone][user] < DBL_MAX);
			//cost = _wp[user] * _delta_p[tone][user][user]; // new way 
			//cost = cf(tone,user);

			//if (cost < min) {
			if (_cost[tone][user] < min) {
				if (total_power_constraint_broken(tone,user) || spectral_mask_constraint_broken(tone,user)) {
					_F[tone][user] = 1;
					continue;
				}
				_all_tones_full=false;
				min=_cost[tone][user];
				//min=cost;
				min_tone=tone;
				min_user=user;
					
			}
		}
	}
	
	*channel=min_tone;
	*line=min_user;

	return 0;
}

int mipb_frac::update_wp()
{
	static int times;

	double a=1e4;
	double b=99;
	double c=91902.4;

	for (int user=0;user<lines;user++) {
		//double x=(_p_used[user]*100/_p_budget[user])-p_ave;
		double x = _p_used[user] - _p_ave;
		//_wp[user] = UNdB(2000/(100*(p_ave+0.00001))*(x));
		//_wp[user] = UNdB(10*(_p_used[user]*100/_p_budget[user] - _p_ave));
		//_wp[user] = pow(exp((_p_used[user] - _p_ave)),100);
		//_wp[user] = pow(exp(x),100);
		/*
		if (-c*x <= -709) {	// FIXME to use FP flags to test for exceptions
			_wp[user] = a*x+1.000001;
			return 0;
		}
		else if (-c*x >= 709) {
			_wp[user] = 0.000001;
		}
		else {
			_wp[user] = a*x+1.000001-(a*x+1)*exp(-c*x)/(exp(-c*x) + b);
		}
		*/
		/*
		if (_wp[user] < 0) {
			_wp[user] = 10*DBL_MIN;
			times++;
			printf("wp was less than 0 %d times\n",times);
		}*/
		//_wp[user] = UNdB(1e6*x);
		//fprintf(_wp_graph,"%g\t%g\n",_p_used[user]-_p_ave,_wp[user]);
		//
	
		// step
		/*
		if (x < 0) {
			_wp[user]=1e-10;
		}
		else if (x > 0) {
			_wp[user]=1e10;
		}
		else {
			_wp[user]=1;
		}*/

		// exponential 
		//_wp[user] = exp(100*x);

		//linear 
		//printf("x = %8.6g\n",x);
		double grad = 1/(2*fabs(x));
		//printf("grad = %8.6g\n",grad);
		_wp[user] = grad*x+1;
		//_wp[user] = 1e3*x+1;
		//_wxt[user] = grad/2*x+1;
		//printf("wp[%d] = %8.6g\n",user,_wp[user]);
		//getchar();

		/*
		_wp[user] = pow(10,1000000*x);

		if (!(_wp[user] <=DBL_MAX)) {
			printf("Whoops! wp = %g x = %g c*x = %g\n",_wp[user],x,c*x);
			_wp[user]=DBL_MAX;
		}
		if (_wp[user] < 0) {
			printf("Whoops! wp = %g x = %g c*x = %g\n",_wp[user],x,c*x);
		}
		*/
		assert(_wp[user] <= DBL_MAX);
		assert(_wp[user] >= 0);
	}

	return 0;

}


double mipb_frac::cf(int tone,int line_id)
{

	double cost=0.0;

	for (int user=0;user<lines;user++) {
		if (user == line_id) {
			cost+=_wp[user]*_delta_p[tone][line_id][user];		// power for the actual bit loaded
		}
		else {
			
			//if (_delta_p[tone][line_id][user] == 0) {	// victim line has no bits at the minute
				double b[lines];
				double p[lines];
				double bits = no_fext_bits(tone,user,minus_40_dbmhz_watts);
				p[line_id]=minus_36_5_dbmhz_watts;
				p[user]=minus_36_5_dbmhz_watts;
				calculate_b_vector_from_psd(p,9.95,tone,b);
				double lost_bits = bits - b[user];
				if (bits < _bit_inc) {
					//printf("Cant load anymore than the minimum bits on line %d tone %d\n",user,tone);
					//getchar();
					cost+=0;
				}
				else {
					cost+=_w[user]*lost_bits;
				}
				//printf("Loading line is %d\n",line_id);
				//printf("No FEXT bits = %6.4g\n",bits);
				//b[line_id] = _b[tone][line_id]+_bit_inc;		// loading users current value
				//p[user] = minus_40_dbmhz_watts;
				//calculate_b_vector_one_bit_set(p,9.95,tone,line_id,b);
				//printf("Bits after loading on line %d = %6.4g\n",line_id,b[user]);
				//cost+=_w[user]*(bits-b[user]);
				//printf("Adding to cost value of %6.4g for loading on line %d\n",_w[user]*(bits-b[user]),line_id);
				//getchar();
				//else {
				//	cost+=bits*bits;
				//}
			//}
			//else {
				cost+=_wp[user]*_delta_p[tone][line_id][user];		// power for the extra crosstalk incurred 0.37
			//}
		}	
	
	}
	//return _w[line_id]*cost;
	//printf("cost was %6.4g\n",cost);
	//getchar();
	return cost;
}

double mipb_frac::var_cf(int tone,int line_id)
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

bool mipb_frac::total_power_constraint_broken(int tone,int user)
{
        for (int i=0;i<lines;i++) {
                //if (p_used[i] + _delta_p[tone][user][i] > _p_budget[i]) {
                if (_p_used[i] + _delta_p[tone][user][i] > _p_budget[i]) {
                	//printf("adding a bit on line %d tone %d will break the power constraint on line %d\n",user,tone,i);
                        return true;
		}
	}

	return false;
}
              
int mipb_frac::update_power(int tone, int user)
{

        for (int i=0;i<lines;i++) {
                if (i == user)
			_last_bit_cost=_delta_p[tone][user][i];
		_p[tone][i]+=_delta_p[tone][user][i];
                _p_used[i]+=_delta_p[tone][user][i];
		_p_ave+=_delta_p[tone][user][i]/lines;
        }
	return 0;
}


int mipb_frac::calc_delta_p(int tone,int line_id)
{
	
	double old_p_tot=0.0,new_p_tot=0.0;

	double b[lines];
	double p[lines];
	double old_p[lines];
	double g[lines];


	for (int user=0;user<lines;user++) {
		g[user] = line_array[user]->gamma[line_array[user]->service[tone]];
		//old_p_tot+=_p[tone][user];
		old_p[user]=_p[tone][user];
		b[user] = _b[tone][user];
	}

	b[line_id]+=_bit_inc;;

	calculate_psd_vector(b,g,tone,p);

	for (int user=0;user<lines;user++) {
		if (p[user] < 0) {
			_F[tone][line_id]=1;
			return 1;
		}
		_delta_p[tone][line_id][user] = p[user] - old_p[user];
	}	
	
	return 0;
}




bool mipb_frac::spectral_mask_constraint_broken(int tone,int user)
{
        if(_spectral_mask_on)
                if (_p[tone][user] + _delta_p[tone][user][user] > _spectral_mask)
                        return true;
        return false;
}




void mipb_frac::init_lines()
{

        int tone,user;
        struct line* current;

        for (user=0;user<lines;user++) {
                current=get_line(user);
		current->is_frac=true;
		strcpy(current->loading_algo,"MIPB FRAC");
                for (tone=0;tone<DMTCHANNELS;tone++) {
                        current->b[tone]=_b[tone][user];
                        current->_b[tone]=_b[tone][user];
                        current->psd[tone]=watts_to_dbmhz(_p[tone][user]);
                }
        }


}

void mipb_frac::write_current_stats(int tone,int user)
{

	char fn[100];
	char fn1[100];
	char fn2[100];
	FILE *fp;
	FILE *fp1;
	FILE *fp2;

	sprintf(fn,"%s/data/mipb/stats/%d.txt",ROOT_DIR,(int)_total_bits);
	sprintf(fn1,"%s/data/mipb/stats/b_and_p_stats.txt",ROOT_DIR);
	sprintf(fn2,"%s/data/mipb/stats/cost_%d.txt",ROOT_DIR,(int)_total_bits);

	fp = fopen(fn,"w");

	if (fp==NULL) {
		printf("Fuck you!\n");	
		exit(2);
	}

	fp1 = fopen(fn1,"a");

	if (fp1==NULL) {
		printf("Fuck you!\n");	
		exit(2);
	}

	fp2 = fopen(fn2,"w");

	if (fp2==NULL) {
		printf("Fuck you!\n");	
		exit(2);
	}

	
	fprintf(fp,"#%dlines\n",lines);
	fprintf(fp,"#bit added was on tone %d line %d\n",tone,user);

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		fprintf(fp,"%d\t",tone);
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
			fprintf(fp,"%4.2lf\t",_b_total[tone]);
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

	fprintf(fp1,"%4.2lf\t",_total_bits);

	for (int user=0;user<lines;user++) {
		fprintf(fp1,"%4.2lf   ",_b_total[user]);
	}

	for (int user=0;user<lines;user++) {
		fprintf(fp1,"%6.4lg   ",_p_used[user]);
	}

	fprintf(fp1,"\n");


	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			fprintf(fp2,"%6.4g\t",_cost[tone][user]);
		}
		fprintf(fp2,"\n");
	}



	fclose(fp);
	fclose(fp1);
	fclose(fp2);

}

void mipb_frac::write_stats_totals()
{

        char fn[80];
        FILE *fp;

        sprintf(fn,"%s/data/mipb/stats/stats.txt",ROOT_DIR);

        fp = fopen(fn,"w");

        if (fp==NULL)
                exit(2);

        fprintf(fp,"#%dlines\n",lines);
        fprintf(fp,"#%dbits\n",(int)_total_bits);
        fprintf(fp,"#%dtones\n",DMTCHANNELS);

        fclose(fp);

}

void mipb_frac::dump_cost_matrix()
{
	char fn[100];
	static int bit=0;
	FILE *fp,*fp1;

	return;

	sprintf(fn,"/tmp/mipb/cost_%d.txt",bit);

	fp = fopen(fn,"w");

	if (fp == NULL)
		exit(2);


	fp1 = fopen("/tmp/mipb/cost.txt","w");
	if (fp1 == NULL)
		exit(2);

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%6.4g\t",_cost[tone][user]);
			fprintf(fp1,"%6.4g\t",_cost[tone][user]);
		}	
		fprintf(fp,"\n");
		fprintf(fp1,"\n");
	}

	fclose(fp);
	fclose(fp1);

	bit++;

	getchar();	// print graph here!

}

