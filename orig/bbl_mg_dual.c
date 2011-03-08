#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bbl_mg_dual.h"
#include "multiuser_load.h"
#include <algorithm>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#define inv 1
// this is where gamma is changed from dB to regular units
/*
double F(double g,int b) 
{
	return pow(10,(g+3)/10)*(pow(2,b)-1);
	//return pow(10,(g+3)/10)*pow(2,b-1);
}
double F(double g,double b) 
{
	return pow(10,(g+3)/10)*(pow(2,b)-1);
	//return pow(10,(g+3)/10)*pow(2,b-1);
}
*/
#define F(g,b) pow(10,(g+3)/10)*(pow(2,b)-1)
double *calculate_psd_vector(int *b,double *gamma,int channel);
// this is a binary predicate to make std::sort sort in descending order, why do i have to do this??? 
bool DESC(const bbl_entry &left,const bbl_entry &right)
{
	return left.metric > right.metric;
}
bbl_mg_dual::bbl_mg_dual()
{
	mg = new bbl_multiuser_greedy;
	rate_b_target = new int[lines];
	current_rate_b = new int[lines];
	bbl=new struct bbl_entry*[DMTCHANNELS];
	bbl_a=new struct bbl_entry*[DMTCHANNELS];
	bbl_b=new struct bbl_entry*[DMTCHANNELS];
	for (int i=0;i<DMTCHANNELS;i++) {
		bbl[i] = new struct bbl_entry[lines];	
		bbl_a[i] = new struct bbl_entry[lines];	
		bbl_b[i] = new struct bbl_entry[lines];	
	}
}
int bbl_mg_dual::run()
{
	init_bad_boys_league();
	print_bbl(bbl);
	getchar();
	print_bbl(bbl_a);
	getchar();
	print_bbl(bbl_b);
	for (int i=0;i<lines;i++) {
		rate_b_target[i] = 100;
	}
	while(1) {
		set_services();	
		mg->run();
		if (rate_b_targets_good())
			break;
		migrate_services();
		//print_bbl(bbl_a);
		//print_bbl(bbl_b);
		//getchar();
	}
 	mg->init_lines();		
	set_loading_algo();
	calculate_snr();
	return 0;	
}
void bbl_mg_dual::init_bad_boys_league(void)
{
	for (int i=0;i<DMTCHANNELS;i++) {
                for (int j=0;j<lines;j++) {
                        //this->bbl[j][i].metric = bbl_metric(j,i);
			//this->bbl[j][i].line_id = j;
			bbl[i][j] = bbl_metric(j,i);
                }
        }
	//print_bbl(bbl);
	for (int i=0;i<DMTCHANNELS;i++) {
		//std::sort(bbl[i],bbl[i]+lines,DESC);		
		std::sort(bbl[i],bbl[i]+lines);		// ascending order, this way lowest xtalkers get service A i.e low BER		
	}
	for (int i=0;i<DMTCHANNELS;i++) {
                for (int j=0;j<lines/2;j++) {
                        bbl_a[i][j] = bbl[i][j];
                }
        }
	for (int i=0;i<DMTCHANNELS;i++) {
                for (int j=lines/2,k=0;j<lines;j++,k++) {
                        bbl_b[i][j] = bbl[i][j];
                }
        }
}
struct bbl_entry bbl_mg_dual::bbl_metric(int line_id, int channel)
{
	double sum=0.0;
        int j;
	struct bbl_entry ret;
        for (j=0;j<lines;j++) {
                if (j != line_id) {     // just the xtalk gains
                        sum+=get_xtalk_gain(line_id,j,channel);
                }
        }
	ret.metric = (this->alpha*sum)/(this->beta*get_xtalk_gain(line_id,line_id,channel));
        ret.tone = channel;
	ret.line_id = line_id;
	return ret;
}
void bbl_mg_dual::print_bbl(struct bbl_entry ** l)
{
	//for(int i=0;i<DMTCHANNELS;i++) {
	for(int i=0;i<20;i++) {
		printf("Tone %d\t",i);
		for (int j=0;j<lines;j++) {
			printf("line %d: %4.2e\t",l[i][j].line_id,l[i][j].metric);
		}
		printf("\n");
	}
	printf("\n");
}
bool bbl_mg_dual::rate_b_targets_good()
{
	int rate_b[lines];
	int line_id;
	bool ret=true;
	for(int user=0;user<lines;user++) {
		rate_b[user]=0;
	}
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int i=lines-1;bbl_b[tone][i].line_id != -1 && i>=0;i--) {
			line_id=bbl_b[tone][i].line_id;
			rate_b[line_id]	+= mg->b[tone][line_id];
		}
	}
	for (int user=0;user<lines;user++) {
		if (abs(rate_b[user]-rate_b_target[user]) > e)
			ret = false;
		printf("line %d rate b = %d\n",user,rate_b[user]);
		current_rate_b[user]=rate_b[user];
	}
	return ret;
}
void bbl_mg_dual::migrate_services()
{
	struct bbl_entry *bbl_ent1;
	struct bbl_entry *bbl_ent2;
	for (int user=0;user<lines;user++) {
		if (current_rate_b[user] > rate_b_target[user]) {
			// find min metric value in bbl_b table for that user and move to bbl_a
			bbl_ent1 = find_bbl_b_min(user);
			printf("min tone for user %d was %d\n",user,bbl_ent1->tone);
		 	move_to_bbl_a(bbl_ent1);
		}
		else if (current_rate_b[user] < rate_b_target[user]){
			// find max metric value in bbl_a table for that user and move to bbl_b
			bbl_ent2 = find_bbl_a_max(user);
			printf("max tone for user %d was %d\n",user,bbl_ent2->tone);
			move_to_bbl_b(bbl_ent2);
		}
	}
}
struct bbl_entry *bbl_mg_dual::find_bbl_b_min(int user)
{
	double min=1e100;
	int min_tone;
	struct bbl_entry *ret = NULL;
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int i=lines-1;bbl_b[tone][i].line_id != -1 && i>=0;i--) {
			if (bbl_b[tone][i].line_id == user) {
				if (bbl_b[tone][i].metric < min) {
					min = bbl_b[tone][i].metric;
					min_tone = tone;
					ret = &bbl_b[tone][i];	
				}
			}
		}
	}
	return ret;
}
struct bbl_entry *bbl_mg_dual::find_bbl_a_max(int user)
{
	double max=0;
	int max_tone=-1;
	struct bbl_entry *ret = NULL;
	for (int tone=0;tone<DMTCHANNELS;tone++) {
                for (int i=0;bbl_a[tone][i].line_id != -1;i++) {
			if (bbl_a[tone][i].line_id == user) {
				if (bbl_a[tone][i].metric > max) {
					max = bbl_a[tone][i].metric;
					max_tone = tone;
					ret = &bbl_a[tone][i];
				}
			}
                }
        }
        return ret;
}
void bbl_mg_dual::move_to_bbl_a(struct bbl_entry *bbl_ent)
{
	int tone=bbl_ent->tone;
	for (int i=0;i<lines;i++) {
		if (bbl_a[tone][i].line_id == -1) {		// found a free spot starting from the top
			bbl_a[tone][i] = *bbl_ent;	// should assign free spot to bbl_ent
			break;
		}
	}
	null_entry(bbl_ent);		// nullify this entry in bbl_b
	std::sort(bbl_b[tone],bbl_b[tone]+lines);
}
void bbl_mg_dual::move_to_bbl_b(struct bbl_entry *bbl_ent)
{
	for (int i=lines-1;i>=0;i--) {
		if (bbl_b[bbl_ent->tone][i].line_id == -1) {		// found a free spot
			bbl_b[bbl_ent->tone][i] = *bbl_ent;	// assign to value of bbl_ent
			break;
		}
	}
	null_entry(bbl_ent);		// nullify this entry in bbl_a
}
void bbl_mg_dual::null_entry(struct bbl_entry* bbl_ent)
{
	bbl_ent->line_id=-1;
	bbl_ent->metric=0;
	bbl_ent->tone=-1;
}
void bbl_mg_dual::set_services()
{
	int tone,i,user;
	struct line *line_array[lines];
	for (user=0;user<lines;user++) {		// get pointers to all struct lines
		line_array[user]=get_line(user);
	}
	for (tone=0;tone<DMTCHANNELS;tone++) {	
		for(i=0;bbl_a[tone][i].line_id != -1 && i<lines;i++) {
			line_array[bbl_a[tone][i].line_id]->service[tone]=0;
		}	
	}
	for (tone=0;tone<DMTCHANNELS;tone++) {	
		for(i=lines-1;bbl_b[tone][i].line_id != -1 && i>=0;i--) {
			line_array[bbl_b[tone][i].line_id]->service[tone]=1;
		}	
	}
}
void bbl_mg_dual::set_loading_algo()
{
	int user;
	struct line* current;
	for (user=0;user<lines;user++) {
		current=get_line(user);
		strcpy(current->loading_algo,"bbl_mg_dual");
	}
}
