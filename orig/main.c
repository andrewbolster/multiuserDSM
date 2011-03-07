#include "multiuser_load.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "scenarios.h"
#include <signal.h>
//#include "loading.h"
#include <unistd.h>
#include <time.h>
#include <boost/regex.hpp>

#include "psd_cache.h"

bool isb_debug=true;

double minus_40_dbmhz_watts;
double minus_140_dbmhz_watts;
double minus_36_5_dbmhz_watts;
struct line **line_array;

bool write_results=false;
bool write_channel_matrix_and_exit=false;

int *rate_targets=NULL;
double *initial_weights=NULL;

int rate_tolerance;

psd_cache *cache=NULL;
double grad_step_size=NOT_SET;
double zeta=1.1;

void sig_handler(int signum)
{
	isb_debug=true;
	getchar();
}

void parse_rate_string(const char *rate_string,int *);
void parse_weight_string(const char * weight_string,double *init_weights);
int get_algo(const char * algostr);

int main(int argc,char **argv)
{

	struct sigaction new_action;
	
	new_action.sa_handler = sig_handler;
	new_action.sa_flags = 0;
	sigaction(3,&new_action,NULL);

	minus_36_5_dbmhz_watts=dbmhz_to_watts(-36.5);
	minus_40_dbmhz_watts=dbmhz_to_watts(-40);
	minus_140_dbmhz_watts=dbmhz_to_watts(-140);

	/*
	if (argc > 1) {
		if (strcmp(argv[1],"-r") == 0) {
			write_results=true;
		}
		else if (strcmp(argv[1],"-c") == 0) {
			write_channel_matrix_and_exit=true;
		}
	}*/

#ifdef VDSL_UPSTREAM
	DMTCHANNELS=1148;
#else
	DMTCHANNELS=480;
#endif
	beta_model=false;
	psd_caching=true;
	read_n=false;
	threads=1;
	rate_tolerance=80;
	zeta=0.1;
	int c,algo;
	char *scenario_file = NULL;
	char *algostr = NULL;

	rate_targets = new int[50];
	initial_weights = new double[50];

	for (int i=0;i<50;i++) {
		rate_targets[i]=NOT_SET;
		initial_weights[i] = NOT_SET;
	}

	while ((c = getopt(argc,argv,"wcs:a:bpn:k:t:r:W:S:R:z:")) !=-1) {
		switch(c) {
			case 'w':
				write_results=true;
				break;
			case 'c':
				write_channel_matrix_and_exit=true;
				break;
			case 's':
				scenario_file = optarg;
				break;
			case 'a':
				algostr = optarg;
				algo = get_algo(algostr);
				break;
			case 'b':
				beta_model=true;
				break;
			case 'p':
				psd_caching=false;
				break;
			case 'n':
				read_n=true;
				n = atoi(optarg);
				break;
			case 'k':
				DMTCHANNELS=atoi(optarg);
				break;
			case 't':
				threads=atoi(optarg);
				break;
			case 'r':
				parse_rate_string(optarg,rate_targets);
				break;
			case 'W':
				parse_weight_string(optarg,initial_weights);
				break;
			case 'S':
				grad_step_size = atof(optarg);
				break;
			case 'R':
				rate_tolerance=atoi(optarg);
				break;
			case 'z':
				zeta=atof(optarg);
				break;
			default:
				exit(1);
		}
	}

	if (scenario_file == NULL) {
		printf("Enter a scenario on the command line with the -s switch\n");
		exit(1);
	}
	if (algostr == NULL) {
		printf("Enter an algorithm on the command line with the -s switch\n");
		exit(1);
	}
	if (algo == -1) {
		printf("I don't know this algorithm\n");
		exit(1);
	}


	if (psd_caching && algo!=BB)  {
		cache = new psd_cache;
	}

	srand((unsigned)time(NULL));
	
	load_network(algo,scenario_file);

	//load_network(MIPB_RR,3lines_upstream1);
	
	if (psd_caching && algo!=BB)  {
		delete cache;
	}

        return 0;

}

void parse_rate_string(const char * rate_string,int *rate_targets)
{
	int str_index=0;
	while(1) {

		char line[3];
		char rate[6];		

		int line_index=0;
		while(1) {
			if (rate_string[str_index] == '=') {
				line[line_index] = '\0';
				break;
			}
			line[line_index] = rate_string[str_index];		
			line_index++;
			str_index++;
		}
		str_index++;
		
		//printf("line = %d\n",atoi(line));

		int rate_index=0;
		while(1) {
			if (rate_string[str_index] == ',' || rate_string[str_index] == '\0') {
				rate[rate_index] = '\0';
				break;
			}
			rate[rate_index] = rate_string[str_index];
			rate_index++;
			str_index++;
		}

		//printf("rate = %d\n",atoi(rate));

		rate_targets[atoi(line)]=atoi(rate);

		if (rate_string[str_index] == '\0')
			break;

		str_index++;
	}	

	//printf("%s\n",rate_string);
}


void parse_weight_string(const char * weight_string,double *init_weights)
{
	int str_index=0;
	while(1) {

		char line[3];
		char weight[15];

		int line_index=0;
		while(1) {
			if (weight_string[str_index] == '=') {
				line[line_index] = '\0';
				break;
			}
			line[line_index] = weight_string[str_index];		
			line_index++;
			str_index++;
		}
		str_index++;
		
		//printf("line = %d\n",atoi(line));

		int weight_index=0;
		while(1) {
			if (weight_string[str_index] == ',' || weight_string[str_index] == '\0') {
				weight[weight_index] = '\0';
				break;
			}
			weight[weight_index] = weight_string[str_index];
			weight_index++;
			str_index++;
		}

		//printf("rate = %d\n",atoi(rate));

		init_weights[atoi(line)]=atof(weight);

		if (weight_string[str_index] == '\0')
			break;

		str_index++;
	}	

}

int get_algo(const char * algostr)
{
	if (strcmp(algostr,"BB") == 0) {
		return BB;
	}
	else if (strcmp(algostr,"MIPB") == 0) {
		return MIPB;
	}
	else if (strcmp(algostr,"IWF") == 0) {
		return IWF;
	}
	else if (strcmp(algostr,"ISB") == 0) {
		return ISB_THREE_G;
	}
	else if (strcmp(algostr,"OSB") == 0) {
		return OSB;
	}
	else if (strcmp(algostr,"E_GREEDY") == 0) {
		return E_GREEDY;
	}
	else if (strcmp(algostr,"GREEDY") == 0) {
		return GREEDY;
	}
	else if (strcmp(algostr,"GREEDY_ADAPTIVE") == 0) {
		return GREEDY_ADAPTIVE;
	}
	else if (strcmp(algostr,"GREEDY_GRAD") == 0) {
		return GREEDY_GRAD;
	}
	else if (strcmp(algostr,"NSECTION") == 0) {
		return NSECTION;
	}
	else if (strcmp(algostr,"NGRAD") == 0) {
		return NGRAD;
	}
	else if (strcmp(algostr,"GREEDY_OLD") == 0) {
		return GREEDY_OLD;
	}
	return -1;
}
	


//int lines = 0;
        //double *channel_matrix;
	//int i,j,k,tot;
	//char tag[100];
	//double ave;

	/*for (k=0;k<2;k++) {
		for (j=0;j<lines;j++) {
			for (i=0;i<lines;i++) {
				printf("%e\t",*(channel_matrix + k + (DMTCHANNELS*j) + (DMTCHANNELS * lines * i)));
			}
		printf("\n");
		}
		printf("\n");
	}*/

	
	
/*
	am_load_ra(0,0.1,0);
	calculate_snr();
	write_line_stats(0,NULL);
	
	am_load_ra(1,0.1,0);
	calculate_snr();
	write_line_stats(1,NULL);
	write_line_stats(0,"after");
*/

	//find_rate_region();
/*
	for (i=0;i<10;i++) {
		//sprintf(tag,"iteration%d",i);
		printf("%d\t",i);
		printf("%d\t",am_load_ra(0,0.1,0));
		calculate_snr();
		//write_line_stats(0,tag);
		printf("%d\n",am_load_ra(1,0.1,0));	
		calculate_snr();
		//write_line_stats(1,tag);
		//am_load(2,0);
		//calculate_snr(channel_matrix,lines);
	}
*/

	//am_load_fm(0,1400,0);
	//dbw_cbr_1(0,9.95,9.95,700,700);
	//calculate_snr();	
/*	

	for (i=0;i<10;i++) {	
		dbw_cbr(0,9.95,9.95,400,400);
		calculate_snr();
		write_line_stats_dual(0,NULL);
		sprintf(tag,"/home/alastair/multiuser/scripts/graph-dual-bits-psd.sh -x ../data/line0_dualstats.txt -t %d",i);
		system(tag);
		dbw_cbr(1,9.95,9.95,400,400);
		calculate_snr();
	}

	write_line_stats_dual(0,NULL);

	free(channel_matrix);
*/

