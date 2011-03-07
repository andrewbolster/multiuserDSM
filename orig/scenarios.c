#include "multiuser_load.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include "bbl_mg_dual.h"
#include "scenarios.h"
//#include "osb_2.h"
#include "multiuser_greedy.h"
#include "multiuser_weighted.h"
#include "multiuser_new.h"
//#include "multiuser_fairer.h"
#include "isb_new.h"
#include "osb_new.h"
#include "osb_bb.h"
#include "osb_bb_frac.h"
#include "mipb.h"
#include "mipb_frac.h"
#include "isb3g.h"
#include "isb3g_frac.h"
#include "greedy_nsection.h"

double *channel_matrix;
int lines;
char scen[50];
extern bool write_results;
extern bool write_channel_matrix_and_exit;
int threads;

void load_network(int loading_algo,char *network)
{
	char file[50];
	sprintf(file,"%s.txt",network);
	struct timeval a,b;
	double running_time;

	
	lines = read_network(file);
	line_array = new struct line*[lines];

	strcpy(scen,network);

	channel_matrix = calc_channel_matrix(lines);
	check_normalised_xtalk_gains();
	//getchar();
	//print_channel_matrix(1);

	set_background_noise(-140);

	//set_alien_xtalk_model(MODEL_A);
	set_alien_xtalk_model(NONE);

	print_network(lines);

	if (write_channel_matrix_and_exit) {
		write_channel_matrix(network);
		exit(0);
	}

	//print_xtalkers(lines);

	for (int i=0;i<lines;i++) {
                line_array[i]=get_line(i);
                line_array[i]->set_line_single_qos(9.95);
#ifdef ADSL_DOWNSTREAM
		line_array[i]->no_fext_b_total = no_fext_ra(i,0.111,&line_array[i]->waterfill_level,&line_array[i]->active_channels);
                //line_array[i]->_no_fext_b_total = no_fext_ra_frac(i,0.111);
#endif
#ifdef VDSL_UPSTREAM
		line_array[i]->no_fext_b_total = no_fext_ra(i,0.02818,&line_array[i]->waterfill_level,&line_array[i]->active_channels);
                //line_array[i]->_no_fext_b_total = no_fext_ra_frac(i,0.02818);
#endif
		printf("Water fill level on line %d is %lf active on %d channels\n",i,watts_to_dbmhz(line_array[i]->waterfill_level),line_array[i]->active_channels);
        }


	switch(loading_algo) {
		case MULTIUSER_GREEDY: {
			char algo[50],tag[50],info[200];
                        strcpy(algo,"MULTIUSER_GREEDY");
                        sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			multiuser_greedy* mg;
			mg = new multiuser_greedy;
			gettimeofday(&a,NULL);
			mg->run();
			gettimeofday(&b,NULL);
			calculate_snr();
			print_summary();
			break;
		}

		case BB: {
			char algo[50],tag[255],info[200];
			strcpy(algo,"BB");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			osb_bb* o;
			o = new osb_bb(threads);
			o->_p_tol = 0.001;
			o->_rate_tol = rate_tolerance;
			o->_min_sl = 50;
			o->_sw = 2e-5;
			//o->_rate_targets[0]=530;
			//o->_rate_targets[0]=14500;
			//o->_w[0]=0.10;
			//o->_w[1]=0.24;	
			//o->_rate_targets[2]=3000;
			/*o->_rate_targets[2]=3000;
			o->_w[0]=0.22;
			o->_w[1]=0.34-0.22;
			*/

			for (int user=0;user<lines;user++) {
				o->_rate_targets[user]=rate_targets[user];
				if (initial_weights[user]!=NOT_SET)
					o->_w[user] = initial_weights[user];
			}

			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				o->_p_budget[user]=0.02818;
			}
			#endif
			o->_dynamic_lambda=true;
			sprintf(info,"p tol=%lf\n\nmin sl=%lf\n\ndynamic lambda=%d\n\n",o->_p_tol,o->_min_sl,o->_dynamic_lambda);
			gettimeofday(&a,NULL);
			o->run();	
			//o->bisect_l();	
			gettimeofday(&b,NULL);
			calculate_snr();
			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,info);
			}
			print_summary();
			print_network(lines);
			delete o;
			delete_network();
                        break;
		}
		case BB_RR: {
			char algo[50],tag[50],info[200];
			strcpy(algo,"BB");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			
			FILE *rr_fp;
			char rr_fn[200];
			sprintf(rr_fn,"../data/%s_rate_region_3000.txt",tag);
			if ((rr_fp=fopen(rr_fn,"w")) == NULL) {
				printf("Feck off\n");
				exit(2);
			}

			fprintf(rr_fp,"0\t%d\n",line_array[1]->no_fext_b_total);
			fclose(rr_fp);

			gettimeofday(&a,NULL);
			for (double w=0.00;w<=0.34;w+=0.02) {		
				osb_bb* o;
				o = new osb_bb(threads);
				o->_p_tol = 0.001;
				o->_min_sl = 50;
				o->_rate_tol = 10;
				o->_w[0]=w;
				o->_w[1]=0.34-w;
				o->_rate_targets[2]=3000;		// rate target on short line
				for (int user=0;user<lines;user++) {
					o->_p_budget[user]=0.02818;	// vdsl power budget
				}
				o->_dynamic_lambda=true;
				sprintf(info,"p tol=%lf\n\nmin sl=%lf\n\ndynamic lambda=%d\n\n",o->_p_tol,o->_min_sl,o->_dynamic_lambda);
				o->run();	
				//o->bisect_l();	
				calculate_snr();
				/*if (write_results) {
					write_all_line_stats(tag);
					write_all_b_and_p_stats(tag);
					write_channel_matrix(tag);
					generate_graphs(tag);
					create_results_summary_file(tag,info);
				}*/
				print_summary();
				print_network(lines);
				delete o;
				printf("Finished w = %lf\n",w);
				//getchar();
				
				if ((rr_fp=fopen(rr_fn,"a")) == NULL) {
					printf("Feck off\n");
					exit(2);
				}
				fprintf(rr_fp,"%d\t%d\n",line_array[0]->rate[0],line_array[1]->rate[0]);
				fclose(rr_fp);
			}
			gettimeofday(&b,NULL);

			if ((rr_fp=fopen(rr_fn,"a")) == NULL) {
				printf("Feck off\n");
				exit(2);
			}
			fprintf(rr_fp,"%d\t0\n",line_array[0]->no_fext_b_total);
			
			fclose(rr_fp);
			delete_network();
                        break;
		}
		case BB_FRAC: {
			char algo[50],tag[50],info[200];
			strcpy(algo,"BB_FRAC");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			osb_bb_frac* o;
                        o = new osb_bb_frac;
                        o->_p_tol = 0.001;
                        o->_inc=0.1;	// cant really make this smaller than 0.001
			//o->_min_sl = 50;
			o->_w[0]=0.765;
			o->_w[1]=1-o->_w[0];
			o->_dynamic_lambda=true;
			sprintf(info,"p tol=%lf\n\nmin sl=%lf\n\ndynamic lambda=%d\n\n",o->_p_tol,o->_min_sl,o->_dynamic_lambda);
			gettimeofday(&a,NULL);
			o->run();	
			//o->bisect_l();	
                        gettimeofday(&b,NULL);
			calculate_snr();
			if (write_results) {
				write_all_line_stats(tag);
                        	write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
                        	create_results_summary_file(tag,info);
			}
                        print_summary();
                        print_network(lines);
			delete_network();
                        break;
		}
		case BB_FRAC_RR: {
			char algo[50],tag[50],info[200];
			strcpy(algo,"BB_FRAC");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			FILE *rr_fp;
			char rr_fn[200];
			sprintf(rr_fn,"../data/%s_rate_region.txt",tag);
			if ((rr_fp=fopen(rr_fn,"w")) == NULL) {
				printf("Feck off\n");
				exit(2);
			}

			fprintf(rr_fp,"0\t%lf\n",line_array[1]->_no_fext_b_total);
			
			for (double w=0.3;w<=0.95;w+=0.05) {
				osb_bb_frac* o;
				o = new osb_bb_frac;
				o->_p_tol = 0.001;
				o->_inc=0.01;	// cant really make this smaller than 0.001
				//o->_min_sl = 50;
				o->_w[0]=w;
				o->_w[1]=1-o->_w[0];
				o->_dynamic_lambda=true;
				//sprintf(info,"p tol=%lf\n\nmin sl=%lf\n\ndynamic lambda=%d\n\n",o->_p_tol,o->_min_sl,o->_dynamic_lambda);
				gettimeofday(&a,NULL);
				o->run();	
				//o->bisect_l();	
				gettimeofday(&b,NULL);
				calculate_snr();
				/*if (write_results) {
					write_all_line_stats(tag);
					write_all_b_and_p_stats(tag);
					write_channel_matrix(tag);
					generate_graphs(tag);
					create_results_summary_file(tag,info);
				}*/
				print_summary();
				fprintf(rr_fp,"%lf\t%lf\n",line_array[0]->_rate[0],line_array[1]->_rate[0]);
				delete o;
			}
			
			fprintf(rr_fp,"%lf\t0\n",line_array[0]->_no_fext_b_total);
			
			fclose(rr_fp);
                        //print_network(lines);
			//delete_network();
                        break;
		}

		case MULTIUSER_NEW: {
			char algo[50],tag[50];
			strcpy(algo,"mn");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			multiuser_new *mn;
			mn = new multiuser_new;
			mn->rate_targets_off=true;
			//mn->graph_loading=true;
			
			//mn->_p_budget[0]=0.110283;
			//mn->_p_budget[1]=0.111285;
			//mn->_p_budget[2]=0.110801;
			//mn->_p_budget[3]=0.108905;
			//mn->_p_budget[4]=0.111994;
			//mn->_p_budget[5]=0.111207;
			//mn->_p_budget[6]=0.111810;
			//mn->_p_budget[7]=0.110757;
			gettimeofday(&a,NULL);
			mn->run_a();
			gettimeofday(&b,NULL);
			calculate_snr();
			if (write_results) {
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			print_network(lines);
			print_summary();
			delete_network();
			break;
		}
		case MIPB: {
			char algo[50],tag[300];
			strcpy(algo,"mipb3g");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			mipb *m;
			m = new mipb(threads);
			//m->_w[0]=UNdB(-40);
			//m->_w[1]=1;			
			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				m->_p_budget[user]=0.02818;	// vdsl power budget
			}
			#endif
			//m->_is_frac=true;
			//m->_rate_targets[0]=14500;
			//m->_w[0]=5e-6;
			//m->_rate_targets[2]=3000;
			//m->_rate_targets[0]=760;
			//m->_w[2]=1.5258789062e-04;
			gettimeofday(&a,NULL);
			m->_graph_loading=false;
			m->_greedy=false;
			m->_old=false;
			//m->_num_threads=threads;
			//m->_show_solution=true;
			m->run();
			gettimeofday(&b,NULL);
			calculate_snr();
			//check_all_margins(0.1);
			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			delete m;
			print_network(lines);
			print_summary();
			delete_network();
			break;
		}
		case E_GREEDY: {
			char algo[50],tag[300];
			strcpy(algo,"e_greedy");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			mipb *m;
			m = new mipb(threads);
			//m->_w[0]=UNdB(-40);
			//m->_w[1]=1;			
			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				m->_p_budget[user]=0.02818;	// vdsl power budget
			}
			#endif
			//m->_is_frac=true;
			//m->_rate_targets[0]=14500;
			//m->_w[0]=5e-6;
			//m->_rate_targets[2]=3000;
			//m->_rate_targets[0]=6000;
			//m->_w[2]=1.5258789062e-04;
			m->_rate_tol=20;
			for (int user=0;user<lines;user++) {
				m->_rate_targets[user]=rate_targets[user];
			}
			gettimeofday(&a,NULL);
			m->_graph_loading=false;
			m->_greedy=true;
			m->_old=false;
			m->_simple_search=false;
			m->_search1=true;
			m->_search2=false;
			//m->_num_threads=threads;
			//m->_show_solution=true;
			m->run();
			gettimeofday(&b,NULL);
			calculate_snr();
			//check_all_margins(0.1);
			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			delete m;
			print_network(lines);
			print_summary();
			delete_network();
			break;
		}
		case GREEDY: {
			char algo[50],tag[300];
			strcpy(algo,"greedy_bisection");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			mipb *m;
			m = new mipb(threads);
			//m->_w[0]=UNdB(-40);
			//m->_w[1]=1;			
			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				m->_p_budget[user]=0.02818;	// vdsl power budget
			}
			#endif
			//m->_is_frac=true;
			//m->_rate_targets[0]=14500;
			//m->_w[0]=5e-6;
			//m->_rate_targets[2]=3000;
			//m->_rate_targets[0]=6000;
			//m->_w[2]=1.5258789062e-04;
			for (int user=0;user<lines;user++) {
				m->_rate_targets[user]=rate_targets[user];
			}
			gettimeofday(&a,NULL);
			m->_rate_tol=rate_tolerance;
			m->_graph_loading=false;
			m->_greedy=true;
			m->_old=true;
			m->_simple_search=false;
			m->_search1=true;
			m->_search2=false;

			m->_rate_search_type=BISECTION;

			if (grad_step_size != NOT_SET)
				m->_sw=grad_step_size;
			//m->_num_threads=threads;
			//m->_show_solution=true;
			m->run();
			gettimeofday(&b,NULL);
			calculate_snr();
			//check_all_margins(0.1);
			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			delete m;
			print_network(lines);
			print_summary();
			delete_network();
			break;
		}
		case GREEDY_GRAD: {
			char algo[50],tag[300];
			strcpy(algo,"greedy_grad");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			mipb *m;
			m = new mipb(threads);
			//m->_w[0]=UNdB(-40);
			//m->_w[1]=1;			
			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				m->_p_budget[user]=0.02818;	// vdsl power budget
			}
			#endif
			//m->_is_frac=true;
			//m->_rate_targets[0]=14500;
			//m->_w[0]=5e-6;
			//m->_rate_targets[2]=3000;
			//m->_rate_targets[0]=6000;
			//m->_w[2]=1.5258789062e-04;
			for (int user=0;user<lines;user++) {
				m->_rate_targets[user]=rate_targets[user];
			}
			gettimeofday(&a,NULL);
			m->_rate_tol=rate_tolerance;
			m->_graph_loading=false;
			m->_greedy=true;
			m->_old=true;
			m->_simple_search=false;
			m->_search1=true;
			m->_search2=false;

			m->_rate_search_type=GRAD;

			if (grad_step_size != NOT_SET)
				m->_sw=grad_step_size;
			//m->_num_threads=threads;
			//m->_show_solution=true;
			m->run();
			gettimeofday(&b,NULL);
			calculate_snr();
			//check_all_margins(0.1);
			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			delete m;
			print_network(lines);
			print_summary();
			delete_network();
			break;
		}
		case GREEDY_ADAPTIVE: {
			char algo[50],tag[300];
			strcpy(algo,"greedy_adaptive_grad");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			mipb *m;
			m = new mipb(threads);
			//m->_w[0]=UNdB(-40);
			//m->_w[1]=1;			
			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				m->_p_budget[user]=0.02818;	// vdsl power budget
			}
			#endif
			//m->_is_frac=true;
			//m->_rate_targets[0]=14500;
			//m->_w[0]=5e-6;
			//m->_rate_targets[2]=3000;
			//m->_rate_targets[0]=6000;
			//m->_w[2]=1.5258789062e-04;
			for (int user=0;user<lines;user++) {
				m->_rate_targets[user]=rate_targets[user];
			}
			gettimeofday(&a,NULL);
			m->_rate_tol=rate_tolerance;
			m->_graph_loading=false;
			m->_greedy=true;
			m->_old=true;
			m->_simple_search=false;
			m->_search1=true;
			m->_search2=false;

			m->_rate_search_type=ADAPTIVE_GRAD;
			m->_sw=1e-8;
			if (grad_step_size != NOT_SET)
				m->_sw=grad_step_size;
			//m->_num_threads=threads;
			//m->_show_solution=true;
			m->run();
			gettimeofday(&b,NULL);
			calculate_snr();
			//check_all_margins(0.1);
			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			delete m;
			print_network(lines);
			print_summary();
			delete_network();
			break;
		}
		case GREEDY_OLD: {
			char algo[50],tag[300];
			strcpy(algo,"greedy_old");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			mipb *m;
			m = new mipb(threads);
			//m->_w[0]=UNdB(-40);
			//m->_w[1]=1;			
			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				m->_p_budget[user]=0.02818;	// vdsl power budget
			}
			#endif
			//m->_is_frac=true;
			//m->_rate_targets[0]=14500;
			//m->_w[0]=5e-6;
			//m->_rate_targets[2]=3000;
			//m->_rate_targets[0]=6000;
			//m->_w[2]=1.5258789062e-04;
			for (int user=0;user<lines;user++) {
				m->_rate_targets[user]=rate_targets[user];
			}
			gettimeofday(&a,NULL);
			m->_rate_tol=rate_tolerance;
			m->_graph_loading=false;
			m->_greedy=true;
			m->_old=true;
			m->_simple_search=false;
			m->_search1=true;
			m->_search2=false;

			m->_rate_search_type=OLD;
			m->_zeta=zeta;
			m->_sw=1e-8;
			if (grad_step_size != NOT_SET)
				m->_sw=grad_step_size;
			//m->_num_threads=threads;
			//m->_show_solution=true;
			m->run();
			gettimeofday(&b,NULL);
			calculate_snr();
			//check_all_margins(0.1);
			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			delete m;
			print_network(lines);
			print_summary();
			delete_network();
			break;
		}
			
		case NSECTION: {
			char algo[50],tag[300];
			strcpy(algo,"greedy_nsection");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			greedy_nsection *ns;
			ns = new greedy_nsection(threads);
			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				ns->_p_budget[user]=0.02818;	// vdsl power budget
			}
			#endif
			for (int user=0;user<lines;user++) {
				ns->_rate_targets[user]=rate_targets[user];
			}
			gettimeofday(&a,NULL);
			ns->nsection();
			gettimeofday(&b,NULL);
			calculate_snr();
			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			delete ns;
			print_network(lines);
			print_summary();
			delete_network();
			break;
		}
		case NGRAD: {
			char algo[50],tag[300];
			strcpy(algo,"greedy_ngrad_search");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			greedy_nsection *ns;
			ns = new greedy_nsection(threads);
			ns->_grad_search=true;
			ns->_nsection=false;
			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				ns->_p_budget[user]=0.02818;	// vdsl power budget
			}
			#endif
			for (int user=0;user<lines;user++) {
				ns->_rate_targets[user]=rate_targets[user];
			}
			gettimeofday(&a,NULL);
			ns->nsection();
			gettimeofday(&b,NULL);
			calculate_snr();
			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			delete ns;
			print_network(lines);
			print_summary();
			delete_network();
			break;
		}
		case MIPB_RR: {
			char algo[50],tag[50];
			strcpy(algo,"mipb_testing_3");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			FILE *rr_fp;
			char rr_fn[200];
			sprintf(rr_fn,"../data/%s_rate_region.txt",tag);
			if ((rr_fp=fopen(rr_fn,"w")) == NULL) {
				printf("Feck off\n");
				exit(2);
			}

			//fprintf(rr_fp,"%d\t0\n",line_array[0]->no_fext_b_total);
			//for (double w=1e-7;w<=1e-5;w*=2) {
			for (int r=6000;r<=14500;r+=500) {
				mipb *m;
				m = new mipb(1);
				//m->_w[0]=w;
				//m->_w[1]=1;
				for (int user=0;user<lines;user++) {
					m->_p_budget[user]=0.02818;	// vdsl power budget
				}
				gettimeofday(&a,NULL);
				//m->_graph_loading=true;
				m->_rate_targets[0]=r;
				m->_rate_targets[2]=3000;
				m->run();
				gettimeofday(&b,NULL);
				calculate_snr();
				//check_all_margins(0.1);
				/*if (write_results) {
					write_all_line_stats(tag);
					write_all_b_and_p_stats(tag);
					write_channel_matrix(tag);
					generate_graphs(tag);
					create_results_summary_file(tag,NULL);
				}*/
				print_network(lines);
				print_summary();
				//getchar();
				fprintf(rr_fp,"%d\t%d\t%lf\t%lf\n",line_array[0]->rate[0],line_array[1]->rate[0],m->_w[0],m->_w[1]);
				delete m;
			}			

			//fprintf(rr_fp,"0\t%d\n",line_array[1]->no_fext_b_total);
			
			fclose(rr_fp);
			
			delete_network();
			break;
		}
		case MIPB_FRAC: {
			char algo[50],tag[50];
			strcpy(algo,"mipb3g_frac");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			mipb *m;
			m = new mipb(1);
			m->_is_frac=true;
			//m->_w[0]=0.1;
			//m->_w[1]=10;
			//m->_rate_targets[0]=700;
			//m->_bit_inc=0.1;
			gettimeofday(&a,NULL);
			m->_graph_loading=false;
			m->_greedy=false;
			m->_show_solution=true;
			m->_old=false;
			m->run();
			gettimeofday(&b,NULL);
			calculate_snr();
			//check_all_margins(0.1);
			if (write_results) {
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			print_network(lines);
			print_summary();
			delete_network();
			break;
		}
		case MIPB_FRAC_RR: {
			char algo[50],tag[200],rr_fn[200];
			strcpy(algo,"mipb_frac3g");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			
			FILE *rr_fp;
			sprintf(rr_fn,"../data/%s_rate_region.txt",tag);
			if ((rr_fp=fopen(rr_fn,"w")) == NULL) {
				printf("Feck off\n");
				exit(2);
			}

			fprintf(rr_fp,"%lf\t0\n",line_array[0]->_no_fext_b_total);
			
			for(int r=350;r<=850;r+=25) {	
				mipb_frac *m;
				m = new mipb_frac;
				//m->_w[0]=w;
				//m->_w[1]=1;
				m->_bit_inc=0.1;
				gettimeofday(&a,NULL);
				//m->_graph_loading=true;
				//m->_greedy=true;
				m->_rate_targets[0]=r;
				m->run();
				gettimeofday(&b,NULL);
				calculate_snr();
				//check_all_margins(0.1);
				/*if (write_results) {
					write_all_line_stats(tag);
					write_all_b_and_p_stats(tag);
					write_channel_matrix(tag);
					generate_graphs(tag);
					create_results_summary_file(tag,NULL);
				}*/
				//print_network(lines);
				print_summary();
				fprintf(rr_fp,"%lf\t%lf\t%lf\t%lf\n",line_array[0]->_rate[0],line_array[1]->_rate[0],m->_w[0],m->_w[1]);
				delete m;
			}
			
			fprintf(rr_fp,"0\t%lf\n",line_array[1]->_no_fext_b_total);
			
			fclose(rr_fp);
		
			delete_network();
			break;
		}
		case ISB_THREE_G: {
			char algo[50],tag[300],info[200];
			strcpy(algo,"ISB3g");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			isb3g* i;
			i = new isb3g(threads);
			//i->_threadpool = false;
			i->_p_tol = 0.006;
			i->_min_sl = 5;
			i->_rate_tol=rate_tolerance;;
			//i->_w[0]=0.5;
			//i->_w[1]=0.5;
			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				i->_p_budget[user]=0.02818;	// vdsl power budget
			}
			#endif
			//i->_dynamic_lambda=true;
			sprintf(info,"p tol=%lf\n\nmin sl=%lf\n\ndynamic lambda=%d\n\n",i->_p_tol,i->_min_sl,i->_dynamic_lambda);
			gettimeofday(&a,NULL);
			//i->run();	
			for (int user=0;user<lines;user++) {
				i->_rate_targets[user]=rate_targets[user];
				if (initial_weights[user]!=NOT_SET)
					i->_w[user] = initial_weights[user];

			}
			i->run();	
			gettimeofday(&b,NULL);
			calculate_snr();
			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,info);
			}
			print_summary();
			print_network(lines);
			delete i;
			delete_network();
                        break;
		}
		case ISB_THREE_G_RR: {
			char algo[50],tag[50],info[200];
			strcpy(algo,"ISB3g");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			
			FILE *rr_fp;
			char rr_fn[200];
			sprintf(rr_fn,"../data/%s_rate_region.txt",tag);
			if ((rr_fp=fopen(rr_fn,"w")) == NULL) {
				printf("Feck off\n");
				exit(2);
			}

			fprintf(rr_fp,"0\t%d\n",line_array[1]->no_fext_b_total);

			for (double w=0.3;w<0.99;w+=0.01) {
				isb3g* i;
				i = new isb3g(threads);
				i->_p_tol = 0.002;
				i->_min_sl = 5;
				i->_w[0]=w;
				i->_w[1]=1-w;
				//i->_dynamic_lambda=true;
				sprintf(info,"p tol=%lf\n\nmin sl=%lf\n\ndynamic lambda=%d\n\n",i->_p_tol,i->_min_sl,i->_dynamic_lambda);
				gettimeofday(&a,NULL);
				i->run();	
				//i->bisect_l();	
				gettimeofday(&b,NULL);
				calculate_snr();
				/*if (write_results) {
					write_all_line_stats(tag);
					write_all_b_and_p_stats(tag);
					write_channel_matrix(tag);
					generate_graphs(tag);
					create_results_summary_file(tag,info);
				}*/
				print_summary();
				fprintf(rr_fp,"%d\t%d\n",line_array[0]->rate[0],line_array[1]->rate[0]);
				//print_network(lines);
				delete i;
			}

			fprintf(rr_fp,"%d\t0\n",line_array[0]->no_fext_b_total);
			
			fclose(rr_fp);			
			delete_network();
                        break;
		}
		case ISB_THREE_G_FRAC: {
			char algo[50],tag[50],info[200];
			strcpy(algo,"ISB3g_frac");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			isb3g_frac* i;
			i = new isb3g_frac;
			i->_p_tol = 0.001;
			i->_min_sl = 50;
			i->_w[0]=0.5;
			i->_w[1]=0.5;
			i->_bit_inc=0.01;
			i->_dynamic_lambda=true;
			sprintf(info,"p tol=%lf\n\nmin sl=%lf\n\ndynamic lambda=%d\n\n",i->_p_tol,i->_min_sl,i->_dynamic_lambda);
			gettimeofday(&a,NULL);
			i->run();	
			//i->bisect_l();	
			gettimeofday(&b,NULL);
			calculate_snr();
			if (write_results) {
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,info);
			}
			print_summary();
			print_network(lines);
			delete i;
			delete_network();
                        break;
		}
		case ISB_THREE_G_FRAC_RR: {
			char algo[50],tag[50],info[200],rr_fn[200];
			strcpy(algo,"ISB3g_frac");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			
			FILE *rr_fp;
			sprintf(rr_fn,"../data/%s_rate_region.txt",tag);
			if ((rr_fp=fopen(rr_fn,"w")) == NULL) {
				printf("Feck off\n");
				exit(2);
			}

			fprintf(rr_fp,"0\t%lf\n",line_array[1]->_no_fext_b_total);

			for (double w=0.1;w<=0.9;w+=0.1) {
				isb3g_frac* i;
				i = new isb3g_frac;
				i->_p_tol = 0.001;
				i->_min_sl = 50;
				i->_w[0]=w;
				i->_w[1]=1-w;
				i->_bit_inc=0.01;
				i->_dynamic_lambda=true;
				sprintf(info,"p tol=%lf\n\nmin sl=%lf\n\ndynamic lambda=%d\n\n",i->_p_tol,i->_min_sl,i->_dynamic_lambda);
				gettimeofday(&a,NULL);
				i->run();	
				//i->bisect_l();	
				gettimeofday(&b,NULL);
				calculate_snr();
				/*if (write_results) {
					write_all_line_stats(tag);
					write_all_b_and_p_stats(tag);
					write_channel_matrix(tag);
					generate_graphs(tag);
					create_results_summary_file(tag,info);
				}*/
				print_summary();
				fprintf(rr_fp,"%lf\t%lf\n",line_array[0]->_rate[0],line_array[1]->_rate[0]);
				delete i;
			}
			fprintf(rr_fp,"%lf\t0\n",line_array[0]->_no_fext_b_total);
			
			fclose(rr_fp);

			delete_network();
                        break;
		}
		case IWF: {
			double* iwf_p;
			iwf_p = new double[lines];
			int rate_tol=rate_tolerance;
			#ifdef ADSL_DOWNSTREAM
			for (int user=0;user<lines;user++) {
				iwf_p[user]=0.111;
			}
			#endif
			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				iwf_p[user]=0.02818;
			}
			#endif

			char algo[50],tag[300];
                        strcpy(algo,"iwf");
                        sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);		
			calculate_snr();
			gettimeofday(&a,NULL);
			int rate_iters=0;	
			int iters=0;	
			//for (int iters=0;iters<50;iters++) {	
			
			while(1) {
				do {
					for (int user=0;user<lines;user++) {
						if (rate_targets[user]==NOT_SET) {
							am_load_ra(user,iwf_p[user],0);
						} 
						else {
							am_load_fm(user,rate_targets[user],iwf_p[user],0);
						}
						calculate_snr();
					}
					/*
					am_load_fm(2,rate_targets[2],iwf_p[2],0);
						calculate_snr();
					am_load_fm(0,rate_targets[0],iwf_p[0],0);
						calculate_snr();
					//am_load_fm(1,rate_targets[1],iwf_p[1],0);
					//	calculate_snr();
					printf("Iter %d:\n",iters);
					for (int user=0;user<lines;user++) {
						printf("rate[%d] = %d\n",user,line_array[user]->rate[0]);
					}
					*/
					iters++;

					if (iters > 30) {
						printf("Margins will not converge\n");
						exit(0);
					}
					
				}while(check_all_margins(2));
			
				bool targets_met=true;

				rate_iters++;
				printf("Rate iters = %d\n",rate_iters);
				for (int user=0;user<lines;user++) {
					if (rate_targets[user] == NOT_SET)
						continue;
					if (abs(line_array[user]->rate[0]-rate_targets[user]) > rate_tol) {
						targets_met=false;
					}
					printf("Rate[%d] = %d\n",user,line_array[user]->rate[0]);
					printf("\n");
				}

				if (targets_met)
					break;

				if (rate_iters > 10) {
					printf("Can't reach target rates\n");
					exit(0);
				}
				/*
				for (int user=0;user<lines;user++) {
					//if ((rate_targets[user] == NOT_SET) && targets_met=false)
					if ((rate_targets[user] == NOT_SET) && (targets_met == false))
						iwf_p[user]/=1.1;
				}
				*/		

			}

			printf("iters = %d\n",iters);
		
			//}
			/*
			double *p_budget = new double[lines];
			int iwf_tol=20;
			double iwf_eps=0.001;
			for (int user=0;user<lines;user++) {
				p_budget[user] = 0.110;
			}
	
			while (1) {	
				int iters=0;
				do {
					for (int user=0;user<lines;user++) {
						am_load_ra(user,p_budget[user],0);
						calculate_snr();
					}
					iters++;
				}while (iters < 20);


				bool done=true;

				for (int user=0;user<lines;user++) {
					if (rate_targets[user] == NOT_SET)
						continue;
					if (abs(rate_targets[user]-line_array[user]->rate[0]) > iwf_tol) {
						done=false;
						break;
					}
				}

				if (done)
					break;

				for (int user=0;user<lines;user++) {
					printf("Rate[%d] = %d\n",user,line_array[user]->rate[0]);
				}	

				//getchar();

				for (int user=0;user<lines;user++) {
					if (rate_targets[user] == NOT_SET)
						continue;
					if ((line_array[user]->rate[0] - rate_targets[user]) > iwf_tol)
						//p_budget[user] -= iwf_eps;
						p_budget[user]/=1.1;
					if ((rate_targets[user] - line_array[user]->rate[0]) > iwf_tol)
						//p_budget[user] += iwf_eps;
						p_budget[user]*=1.1;
					if (p_budget[user] > 0.110)
						p_budget[user]=0.110;
				}

			}
			*/
			gettimeofday(&b,NULL);
			/*do {
				for (int user=0;user<lines;user++) {
					margin_tighten(user);
					calculate_snr();
				}
			}while(check_all_margins(0.1));
			*/

			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			print_network(lines);
			print_summary();
			delete_network();
			break;
		}
		case ISB: {
			//char tag[100];	
			char algo[50],tag[50];
                        strcpy(algo,"ISB");
                        sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			isb* i;
			i = new isb;
			calculate_snr();
			i->run();
			calculate_snr();
			if (write_results) {
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,NULL);
			}
			print_network(lines);
			print_summary();
			delete i;
			delete_network();
			break;
		}
		case OSB: {
			char algo[50],tag[255],info[200];
			strcpy(algo,"OSB");
			sprintf(tag,"%s_%s_%d",network,algo,DMTCHANNELS);
			osb_bb* o;
			o = new osb_bb(threads);
			o->_p_tol = 0.001;
			o->_rate_tol = 10;
			o->_min_sl = 50;
			o->_sw = 2e-5;
			o->_branch_and_bound=false;
			#ifdef VDSL_UPSTREAM
			for (int user=0;user<lines;user++) {
				o->_p_budget[user]=0.02818;
			}
			#endif

			for (int user=0;user<lines;user++) {
				o->_rate_targets[user]=rate_targets[user];
				if (initial_weights[user]!=NOT_SET)
					o->_w[user] = initial_weights[user];
			}

			o->_dynamic_lambda=true;
			sprintf(info,"p tol=%lf\n\nmin sl=%lf\n\ndynamic lambda=%d\n\n",o->_p_tol,o->_min_sl,o->_dynamic_lambda);
			gettimeofday(&a,NULL);
			o->run();	
			//o->bisect_l();	
			gettimeofday(&b,NULL);
			calculate_snr();
			if (write_results) {
				for (int user=0;user<lines;user++) {
					char rate[10];
					sprintf(rate,"_%d",line_array[user]->rate[0]);
					strcat(tag,rate);
				}
				write_all_line_stats(tag);
				write_all_b_and_p_stats(tag);
				write_channel_matrix(tag);
				generate_graphs(tag);
				create_results_summary_file(tag,info);
			}
			print_summary();
			print_network(lines);
			delete o;
			delete_network();
                        break;
		}

	}

	running_time = ((b.tv_sec - a.tv_sec)*1e6 + (b.tv_usec - a.tv_usec))/1e6;

	printf("Approx running time = %lf seconds\n",running_time);
	printf("%lf\n",running_time);

	printf("Data rates\n");

	for (int user=0;user<lines;user++) {
		printf("%d\t",line_array[user]->rate[0]);
	}
	printf("\n");

}

#ifdef UNDEF

void single_qos_8_lines_1(int scenario)
{

	lines = read_network("8lines_1.txt");
	//int max_rate1,max_rate0;
	line_array = new struct line*[lines];
	//double rate1_target=1;
	
	strcpy(scen,"8_lines_1");


	channel_matrix = calc_channel_matrix(lines);	

	write_channel_matrix("8lines_1");

	set_background_noise(-140);
	set_alien_xtalk_model(MODEL_A);
	//set_alien_xtalk_model(NONE);

	print_network(lines);

	print_channel_matrix(1);
	
	for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
		line_array[i]->set_line_single_qos(9.95);
		//line_array[i]->no_fext_b_total = no_fext_ra(i,0.110,&line_array[i]->waterfill_level);
	}

	switch(scenario) {
		case ISB: {
			//char tag[100];	
			isb* i;
			i = new isb;
			calculate_snr();
			i->rate_target[1]=200;
			i->rate_target[2]=200;
			i->rate_target[3]=200;
			i->rate_target[4]=200;
			i->rate_target[5]=700;
			i->rate_target[6]=700;
			i->rate_target[7]=700;
			i->e=3;
			/*
			i->w[1]=0.7;
			i->w[2]=0.7;
			i->w[3]=0.7;
			i->w[4]=1.3;
			i->w[5]=1.3;
			i->w[6]=1.3;
			i->w[7]=1.3;
			*/
			i->run();
			delete i;
			calculate_snr();
			//write_all_line_stats("8lines_1_isb");
			//write_all_p_stats("8lines_1_isb");
			write_all_b_and_p_stats("8lines_1_isb");
			print_summary();
			delete_network();
			break;
		}
		case OSB: {
			//char tag[100];	
			osb* i;
			i = new osb;
			calculate_snr();
			i->rate_target[1]=125;
			i->rate_target[2]=125;
			i->rate_target[3]=125;
			i->rate_target[4]=250;
			i->rate_target[5]=250;
			i->rate_target[6]=250;
			i->rate_target[7]=250;
			i->e=3;
			i->run();
			delete i;
			calculate_snr();
			write_all_line_stats("8lines_1_osb");
			print_summary();
			delete_network();
			break;
		}
		case BB: {
			osb_bb* o;
			o = new osb_bb;
			o->run();
			calculate_snr();
			print_summary();
			write_all_line_stats("8lines_1_bb");
			write_all_b_and_p_stats("8lines_1_bb");
			delete_network();
			break;
		}
		case MULTIUSER_WEIGHTED: {
			multiuser_weighted *mw;
			mw = new multiuser_weighted;
			
			mw->b_target[0]=1855;
			mw->b_target[1]=1855;
			mw->b_target[2]=1855;
			mw->b_target[3]=1855;
			mw->b_target[4]=2505;
			mw->b_target[5]=2505;
			mw->b_target[6]=2505;
			mw->b_target[7]=2505;
		/*	
			mw->w[0]=1;
			mw->w[1]=1;
			mw->w[2]=1;
			mw->w[3]=1;
			mw->w[4]=12;
			mw->w[5]=12;
			mw->w[6]=12;
			mw->w[7]=12;
		*/	
			
			mw->w[0]=6.5;
			mw->w[1]=6.5;
			mw->w[2]=6.5;
			mw->w[3]=6.5;
			mw->w[4]=34;
			mw->w[5]=34;
			mw->w[6]=34;
			mw->w[7]=34;
			
			mw->w_updates_off=true;
			//mw->rate_targets_off=true;
			mw->run_a();
			delete mw;
			calculate_snr();
			write_all_line_stats("8lines_1_mw");
			write_all_b_and_p_stats("8lines_1_mw");
			//delete mw;
			print_summary();
			delete_network();
			break;
		}
		case MULTIUSER_GREEDY: {
			multiuser_greedy *mg;
			mg = new multiuser_greedy;	
			mg->b_target[4]=660;
			mg->b_target[5]=660;
			mg->b_target[6]=660;
			mg->b_target[7]=660;
			mg->run();
			calculate_snr();
			print_summary();
			delete_network();
			break;
		}
		
		case MULTIUSER_NEW: {
			multiuser_new *mn;
			mn = new multiuser_new;
			mn->rate_targets_off=true;
			mn->graph_loading=true;
			mn->run_a();
			calculate_snr();
			write_all_b_and_p_stats("8lines_1_mn");
			print_summary();
			delete_network();
			break;
		}

	}


}

void single_qos_8_lines_2(int scenario)
{

	//struct line** line_array;
	lines = read_network("8lines_2.txt");
	//int max_rate1,max_rate0;
	line_array = new struct line*[lines];
	//double rate1_target=1;

	channel_matrix = calc_channel_matrix(lines);	

	write_channel_matrix("8lines_2");

	strcpy(scen,"8_lines_2");
	
	set_background_noise(-140);
	set_alien_xtalk_model(MODEL_A);
	//set_alien_xtalk_model(NONE);

	print_network(lines);
	print_case(2,3);

	print_channel_matrix(1);

	for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
		line_array[i]->set_line_single_qos(9.95);
	}

	switch(scenario) {
		case ISB: {
			//char tag[100];	
			isb* i;
			i = new isb;
			calculate_snr();
			i->rate_target[1]=200;
			i->rate_target[2]=200;
			i->rate_target[3]=200;
			i->rate_target[4]=200;
			i->rate_target[5]=700;
			i->rate_target[6]=700;
			i->rate_target[7]=700;
			i->e=3;
			//i->w[0]=2;
			//i->w[2]=0.7;
			//i->w[3]=0.7;
			//i->w[4]=1.3;
			//i->w[5]=0.001;
			//i->w[6]=1.3;
			//i->w[7]=1.3;
			i->p_tol=0.002;
			i->run();
			delete i;
			calculate_snr();
			write_all_line_stats("8lines_2_isb");
			write_all_b_and_p_stats("8lines_2_isb");
			print_summary();
			delete_network();
			break;
		}
		case OSB: {
			//char tag[100];	
			osb* i;
			i = new osb;
			calculate_snr();
			i->rate_target[1]=125;
			i->rate_target[2]=125;
			i->rate_target[3]=125;
			i->rate_target[4]=250;
			i->rate_target[5]=250;
			i->rate_target[6]=250;
			i->rate_target[7]=250;
			i->e=3;
			i->run();
			delete i;
			calculate_snr();
			write_all_line_stats("8lines_2_osb");
			print_summary();
			delete_network();
			break;
		}
		case MULTIUSER_WEIGHTED: {
			multiuser_weighted *mw;
			mw = new multiuser_weighted;
			/*	
			mw->b_target[0]=52;
			mw->b_target[1]=120;
			mw->b_target[2]=410;
			mw->b_target[3]=770;
			mw->b_target[4]=468;
			mw->b_target[5]=1055;
			mw->b_target[6]=464;
			mw->b_target[7]=507;
			*/

			// This set of targets and weights work with w updates off and rate targets on

			mw->b_target[0]=88;
			mw->b_target[1]=235;
			mw->b_target[2]=383;
			mw->b_target[3]=532;
			mw->b_target[4]=789;
			mw->b_target[5]=1289;
			mw->b_target[6]=785;
			mw->b_target[7]=460;
			
			/*
			mw->w[0]=4;
			mw->w[1]=16;
			mw->w[2]=4;
			mw->w[3]=0;
			mw->w[4]=8192;
			mw->w[5]=8192;
			mw->w[6]=8192;
			mw->w[7]=32;
			*/
			
			//mw->w_updates_off=true;
			//mw->rate_targets_off=true;
			mw->e=1;
			mw->run_a();
			delete mw;
			calculate_snr();
			write_all_line_stats("8lines_2_mw");
			
			//delete mw;
		/*	
			mw->w_updates_off=true;
			//mw->rate_targets_off=true;
			mw->run_a();
			delete mw;
			calculate_snr();
			write_all_line_stats("8lines_2_mw_intial_weights");
		*/
			write_all_b_and_p_stats("8lines_2_mw");
			print_summary();
			delete_network();

			break;
		}

		case MULTIUSER_NEW: {
			multiuser_new *mn;
			mn = new multiuser_new;
			mn->rate_targets_off=true;
			//mn->b_target[0]=872;
			//mn->b_target[1]=1025;
			//mw->w[0]=0;
			//mw->w[1]=74;	
			mn->run_a();
			calculate_snr();
			print_summary();
			delete_network();
			break;
		}

		case MULTIUSER_GREEDY: {
			multiuser_greedy *mg;
			mg = new multiuser_greedy;	
			mg->b_target[4]=660;
			mg->b_target[5]=660;
			mg->b_target[6]=660;
			mg->b_target[7]=660;
			mg->run();
			calculate_snr();
			print_summary();
			delete_network();
			break;
		}

	}


}

void single_qos_8_lines_3(int scenario)
{

	//struct line** line_array;
	lines = read_network("8lines_3.txt");
	//int max_rate1,max_rate0;
	line_array = new struct line*[lines];
	//double rate1_target=1;

	channel_matrix = calc_channel_matrix(lines);	

	write_channel_matrix("8lines_3");

	strcpy(scen,"8_lines_3");
	
	set_background_noise(-140);
	set_alien_xtalk_model(MODEL_A);
	//set_alien_xtalk_model(NONE);

	print_network(lines);
	print_case(2,3);

	print_channel_matrix(1);

	for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
		line_array[i]->set_line_single_qos(9.95);
		//line_array[i]->no_fext_b_total = no_fext_ra(i,0.110,&line_array[i]->waterfill_level);
	}

	switch(scenario) {
		case ISB: {
			//char tag[100];	
			isb* i;
			i = new isb;
			calculate_snr();
			i->rate_target[1]=200;
			i->rate_target[2]=200;
			i->rate_target[3]=200;
			i->rate_target[4]=200;
			i->rate_target[5]=700;
			i->rate_target[6]=700;
			i->rate_target[7]=700;
			i->p_tol=0.002;
			i->e=3;
			//i->w[0]=1;
			//i->w[2]=0.7;
			//i->w[3]=0.7;
			//i->w[4]=1.3;
			//i->w[5]=0.001;
			//i->w[6]=1.3;
			//i->w[7]=1.3;
			i->run();
			delete i;
			calculate_snr();
			check_all_margins(0.1);
			write_all_b_and_p_stats("8lines_3_isb");
			print_summary();
			delete_network();
			break;
		}
		case OSB: {
			//char tag[100];	
			osb* i;
			i = new osb;
			calculate_snr();
			i->rate_target[1]=125;
			i->rate_target[2]=125;
			i->rate_target[3]=125;
			i->rate_target[4]=250;
			i->rate_target[5]=250;
			i->rate_target[6]=250;
			i->rate_target[7]=250;
			i->e=3;
			i->run();
			delete i;
			calculate_snr();
			write_all_line_stats("8lines_2_osb");
			print_summary();
			delete_network();
			break;
		}
		case BB: {
			osb_bb* o;
			o = new osb_bb;
			o->run();
			calculate_snr();
			print_summary();
			write_all_line_stats("8lines_3_bb");
			write_all_b_and_p_stats("8lines_3_bb");
			delete_network();
			break;
		}
		case MULTIUSER_WEIGHTED: {
			multiuser_weighted *mw;
			mw = new multiuser_weighted;
			/*	
			mw->b_target[0]=52;
			mw->b_target[1]=120;
			mw->b_target[2]=410;
			mw->b_target[3]=770;
			mw->b_target[4]=468;
			mw->b_target[5]=1055;
			mw->b_target[6]=464;
			mw->b_target[7]=507;
			*/


			mw->b_target[0]=628;
			mw->b_target[1]=860;
			mw->b_target[2]=732;
			mw->b_target[3]=1025;
			mw->b_target[4]=1252;
			mw->b_target[5]=1810;
			mw->b_target[6]=1252;
			mw->b_target[7]=753;
			
			
			mw->w[0]=6000;
			mw->w[1]=2000;
			mw->w[2]=2000;
			mw->w[3]=2000;
			mw->w[4]=10;
			mw->w[5]=10;
			mw->w[6]=10;
			mw->w[7]=0;
			
			
			mw->w_updates_off=true;
			//mw->rate_targets_off=true;
			mw->e=1;
			mw->run_a();
			delete mw;
			calculate_snr();
			write_all_line_stats("8lines_3_mw");
			write_all_b_and_p_stats("8lines_3_mw");
			
			//delete mw;
		/*	
			mw->w_updates_off=true;
			//mw->rate_targets_off=true;
			mw->run_a();
			delete mw;
			calculate_snr();
			write_all_line_stats("8lines_2_mw_intial_weights");
		*/
			write_all_p_stats("8lines_3_mw");
			write_all_b_stats("8lines_3_mw");
			print_summary();
			delete_network();

			break;
		}
		case MULTIUSER_GREEDY: {
			multiuser_greedy *mg;
			mg = new multiuser_greedy;	
			mg->b_target[4]=660;
			mg->b_target[5]=660;
			mg->b_target[6]=660;
			mg->b_target[7]=660;
			mg->run();
			calculate_snr();
			print_summary();
			delete_network();
			break;
		}

		case MULTIUSER_NEW: {
			multiuser_new *mn;
			mn = new multiuser_new;
			//mn->wb[1]=UNdB(-10);	
			//mn->wxt[7]=UNdB(200);	
			//mn->wxt[6]=2;	
			//mn->wxt[7]=2;	
			//mn->wb[0]=1e-150;
			//mn->greedy=true;
			//mn->p_lead[1]=10;
			
			//mn->_p_budget[5]=0.0005;
			//mn->wp_offset[5]=5000;
			//mn->wxt[0]=UNdB(15);
			//mn->wxt[1]=UNdB(2000);
			/*
			mn->wb[0]=UNdB(30);
			mn->wb[1]=UNdB(0);
			mn->wb[2]=UNdB(30);
			mn->wb[3]=UNdB(30);
			mn->wb[4]=UNdB(30);
			mn->wb[5]=UNdB(30);
			mn->wb[6]=UNdB(30);
			mn->wb[7]=UNdB(30);
			*/
			//mn->wxt[2]=UNdB(15);
			//mn->wxt[3]=UNdB(15);
			mn->w_updates_off=true;
			mn->rate_targets_off=true;
			//mn->graph_loading=true;
			mn->e=1;
			mn->run_a();
			delete mn;
			calculate_snr();
			check_all_margins(0.1);
			write_all_line_stats("8lines_3_mn");
			
			//delete mw;
		/*	
			mw->w_updates_off=true;
			mw->rate_targets_off=true;
			mw->run_a();
			delete mw;
			calculate_snr();
			write_all_line_stats("8lines_2_mw_intial_weights");
		*/
			write_all_b_and_p_stats("8lines_3_mn");
			print_summary();
			delete_network();

			break;
		}
		

	}


}


void single_qos_16_lines_3(int scenario)
{

	//struct line** line_array;
	lines = read_network("16lines_3.txt");
	//int max_rate1,max_rate0;
	line_array = new struct line*[lines];
	//double rate1_target=1;

	channel_matrix = calc_channel_matrix(lines);	

	write_channel_matrix("16lines_3");

	strcpy(scen,"16_lines_3");
	
	set_background_noise(-140);
	set_alien_xtalk_model(MODEL_A);
	//set_alien_xtalk_model(NONE);

	print_network(lines);
	print_case(2,3);

	print_channel_matrix(1);

	for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
		line_array[i]->set_line_single_qos(9.95);
		//line_array[i]->no_fext_b_total = no_fext_ra(i,0.110,&line_array[i]->waterfill_level);
	}

	switch(scenario) {
		case ISB: {
			//char tag[100];	
			isb* i;
			i = new isb;
			calculate_snr();
			i->run();
			delete i;
			calculate_snr();
			write_all_b_and_p_stats("16lines_3_isb");
			print_summary();
			delete_network();
			break;
		}

		case MULTIUSER_GREEDY: {
			multiuser_greedy *mg;
			mg = new multiuser_greedy;	
			mg->run();
			calculate_snr();
			print_summary();
			delete_network();
			break;
		}

		case MULTIUSER_NEW: {
			multiuser_new *mn;
			mn = new multiuser_new;
			
			mn->w_updates_off=true;
			mn->rate_targets_off=true;
			mn->e=1;
			mn->run_a();
			delete mn;
			calculate_snr();
			write_all_line_stats("16lines_3_mn");
			
			write_all_b_and_p_stats("16lines_3_mn");
			print_summary();
			delete_network();

			break;
		}
		

	}


}




void single_qos_3_lines_1(int scenario)
{

	struct line** line_array;
	lines = read_network("3lines_1.txt");
	int max_rate1,max_rate0;
	line_array = new struct line*[lines];
	double rate1_target=4;

	channel_matrix = calc_channel_matrix(lines);	

	write_channel_matrix("3lines_1");

	set_background_noise(-140);
	set_alien_xtalk_model(MODEL_A);
	//set_alien_xtalk_model(NONE);

	print_network(lines);


	for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
		line_array[i]->set_line_single_qos(9.95);
	}

	switch(scenario) {
		case SINGLE_IWF: {
			calculate_snr();
			
			for (int i=0;i<5;i++) {
				am_load_fm(1,rate1_target*250,0.110,BOTH_HALVES);
				calculate_snr();
				am_load_ra(0,0.110,BOTH_HALVES);
				calculate_snr();
			}

			write_line_stats(0,"2_lines_near_far_IWF");	
			write_line_stats(1,"2_lines_near_far_IWF");	
			print_summary();
			delete_network();		
			break;
		}
		case OSB: {
			//char tag[100];
			osb* o;
			o = new osb;
			calculate_snr();
			//o->rate_target[1]=rate1_target*250;
			//o->rate_target[2]=rate1_target*250;
			o->e=1;
			o->run();
			delete o;
			calculate_snr();		
			print_summary();
			delete_network();
			break;
		}
		case ISB: {
			//char tag[100];	
			isb* i;
			i = new isb;
			calculate_snr();
			i->rate_target[1]=rate1_target*250;
			i->rate_target[2]=rate1_target*250;
			i->e=1;
			i->run();
			delete i;
			calculate_snr();
			write_all_line_stats("3lines_1_isb");
			print_summary();
			delete_network();
			break;
		}
		case RATE_REGION_IWF: {
			int j=0;
			calculate_snr();

			am_load_ra(1,0.110,BOTH_HALVES);

			calculate_snr();
			max_rate1=line_array[1]->rate[0];
	
			printf("max rate on RT line is %d = approx %d Mbps\n",max_rate1,max_rate1/250);

			printf("0\t%lf\n",(double)(max_rate1*FRAMERATE)/1e6);
		
			for(int rate1=max_rate1/250;rate1>=0;rate1--) {
				//hjunnnnnnnnnnnnnnnnnnnnnnnnnnnnnnp0lkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkko
				//comment by rambo
				do {
					am_load_fm(1,rate1*250,0.110,BOTH_HALVES);	
					calculate_snr();
					am_load_ra(0,0.110,BOTH_HALVES);
					calculate_snr();
					j++;
				//} while(check_all_margins(0.01));
				} while(j<5);
				printf("%lf\t%lf\n",(double)(line_array[0]->rate[0]*FRAMERATE/1e6),
							(double)(line_array[1]->rate[0]*FRAMERATE/1e6));

			}
	
			break;
		}
		case RATE_REGION_OSB: {
			struct osb_2_params* p;
			p = new osb_2_params;
			calculate_snr();	
			am_load_ra(0,0.110,BOTH_HALVES);			
			calculate_snr();
			max_rate0=line_array[0]->rate[0];

			printf("%lf\t0\n",(double)(max_rate0*FRAMERATE/1e6));

			for (int rate0=max_rate0/250;rate0>=0;rate0--) {
				p->osb_2_params_init();
				p->rate0_target=rate0*250;
				osb_2(p);
				calculate_snr();
				printf("%lf\t%lf\n",(double)(line_array[0]->rate[0]*FRAMERATE/1e6),
                                                        (double)(line_array[1]->rate[0]*FRAMERATE/1e6));
			}
				
			break;
		}
		case MULTIUSER_GREEDY: {
			multiuser_greedy *mg;
			mg = new multiuser_greedy;	
			mg->b_target[1]=rate1_target*250;
			mg->run();
			calculate_snr();
			write_line_stats(0,"2_lines_near_far_mg");
			write_line_stats(1,"2_lines_near_far_mg");
			print_summary();
			delete_network();
			break;
		}
		case MULTIUSER_WEIGHTED: {
			multiuser_weighted *mw;
			mw = new multiuser_weighted;
			mw->b_target[1]=rate1_target*250;
			mw->b_target[2]=rate1_target*250;
			mw->run_a();
			calculate_snr();
			write_all_line_stats("3lines_1_mw");
			delete mw;
			print_summary();
			delete_network();
			break;
		}
		/*
		case MULTIUSER_FAIRER: {
			multiuser_fairer *mf;
			mf = new multiuser_fairer;
			mf->b_target[1]=rate1_target*250;
			mf->run();
			calculate_snr();
			write_line_stats(0,"2_lines_near_far_mf");
			write_line_stats(1,"2_lines_near_far_mf");
			print_summary();
			delete_network();
			break;
		}*/
		
	}
}


void xt_power_test()
{

	int b[2];
	double g[2];
	double p[2];
	struct line** line_array;
	lines = read_network("2lines_near_far.txt");
	line_array = new struct line* [lines];
	FILE *fp,*fp1;
	char fn[60];
	char fn1[60];


	channel_matrix = calc_channel_matrix(lines);	

	set_background_noise(-140);
	
	set_alien_xtalk_model(NONE);
	
	for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
		line_array[i]->set_line_single_qos(9.95);
		g[i]=9.95;
	}

	vector_next(b,2,15,true);

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		sprintf(fn,"../data/xt_test/line0_%d.txt",tone);
		sprintf(fn1,"../data/xt_test/line1_%d.txt",tone);
		fp=fopen(fn,"w");
		fp1=fopen(fn1,"w");
		while(vector_next(b,2,15,false)) {
			calculate_psd_vector(b,g,tone,p);
			fprintf(fp,"%d\t%d\t%4.2e\t%4.2lf\n",b[0],b[1],watts_to_dbmhz(p[0]),nxt(1,0,p[1],tone));	
			fprintf(fp1,"%d\t%d\t%4.2e\t%4.2lf\n",b[0],b[1],watts_to_dbmhz(p[1]),nxt(0,1,p[0],tone));	
		}
		
		fclose(fp);
		fclose(fp1);
	}


}



void single_qos_2_lines_near_far(int scenario)
{

	lines = read_network("2lines_near_far.txt");
	int max_rate1,max_rate0;
	line_array = new struct line*[lines];
	double rate1_target=4;

	channel_matrix = calc_channel_matrix(lines);	

	write_channel_matrix("2lines_near_far");

	set_background_noise(-140);
	//set_alien_xtalk_model(MODEL_A);
	set_alien_xtalk_model(NONE);

	print_network(lines);


	for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
		line_array[i]->set_line_single_qos(9.95);
	}

	switch(scenario) {
		case SINGLE_IWF: {
			calculate_snr();
			
			for (int i=0;i<5;i++) {
				am_load_fm(1,rate1_target*250,0.110,BOTH_HALVES);
				calculate_snr();
				am_load_ra(0,0.110,BOTH_HALVES);
				calculate_snr();
			}

			write_line_stats(0,"2_lines_near_far_IWF");	
			write_line_stats(1,"2_lines_near_far_IWF");	
			print_summary();
			delete_network();		
			break;
		}
		case OSB: {
			osb *o;
			o = new osb;
			o->p_tol=0.001;
			o->rate_target[1]=1025;
			//o->w[0]=20;
			//o->w[1]=1;
			calculate_snr();	
			o->run();
			calculate_snr();
			print_summary();
			delete_network();
			break;
		}
		case BB: {
			osb_bb* o;
			o = new osb_bb;
			o->run();
			print_summary();
			delete_network();
			break;
		}
		case OSB_OLD: {	
			char tag[100];	
			struct osb_2_params* p;
			p = new osb_2_params;
			calculate_snr();
			p->osb_2_params_init();
			//p->rate1_target=rate1_target*250;
			p->rate1_target=1025;
			p->e=1;
			//p->w=0.999999999501;
			osb_2(p);
			calculate_snr();
			sprintf(tag,"2_lines_near_far_OSB_%d_%d",p->rate0_target,p->rate1_target);
			write_line_stats(0,tag);
			write_line_stats(1,tag);
			print_summary();
			delete_network();
			break;
		}
		case ISB: {
			//char tag[100];	
			isb* i;
			i = new isb;
			calculate_snr();
			//i->rate_target[1]=810;
			i->p_tol=0.001;
			i->e=1;
			i->w[1]=50;
			i->run();
			delete i;
			calculate_snr();
			print_summary();
			delete_network();
			break;
		}
		case RATE_REGION_IWF: {
			int j=0;
			calculate_snr();

			am_load_ra(1,0.110,BOTH_HALVES);

			calculate_snr();
			max_rate1=line_array[1]->rate[0];
	
			printf("max rate on RT line is %d = approx %d Mbps\n",max_rate1,max_rate1/250);

			printf("0\t%lf\n",(double)(max_rate1*FRAMERATE)/1e6);
		
			for(int rate1=max_rate1/250;rate1>=0;rate1--) {
				//hjunnnnnnnnnnnnnnnnnnnnnnnnnnnnnnp0lkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkko
				//comment by rambo
				do {
					am_load_fm(1,rate1*250,0.110,BOTH_HALVES);	
					calculate_snr();
					am_load_ra(0,0.110,BOTH_HALVES);
					calculate_snr();
					j++;
				//} while(check_all_margins(0.01));
				} while(j<5);
				printf("%lf\t%lf\n",(double)(line_array[0]->rate[0]*FRAMERATE/1e6),
							(double)(line_array[1]->rate[0]*FRAMERATE/1e6));

			}
	
			break;
		}
		case RATE_REGION_OSB: {
			struct osb_2_params* p;
			p = new osb_2_params;
			calculate_snr();	
			am_load_ra(0,0.110,BOTH_HALVES);			
			calculate_snr();
			max_rate0=line_array[0]->rate[0];

			printf("%lf\t0\n",(double)(max_rate0*FRAMERATE/1e6));

			for (int rate0=max_rate0/250;rate0>=0;rate0--) {
				p->osb_2_params_init();
				p->rate0_target=rate0*250;
				osb_2(p);
				calculate_snr();
				printf("%lf\t%lf\n",(double)(line_array[0]->rate[0]*FRAMERATE/1e6),
                                                        (double)(line_array[1]->rate[0]*FRAMERATE/1e6));
			}
				
			break;
		}
		case MULTIUSER_GREEDY: {
			multiuser_greedy *mg;
			mg = new multiuser_greedy;	
			mg->b_target[1]=1025;
			mg->run();
			calculate_snr();
			write_line_stats(0,"2_lines_near_far_mg");
			write_line_stats(1,"2_lines_near_far_mg");
			print_summary();
			delete_network();
			break;
		}
		case MULTIUSER_WEIGHTED: {
			multiuser_weighted *mw;
			mw = new multiuser_weighted;
			//mw->rate_targets_off=true;
			mw->b_target[0]=872;
			mw->b_target[1]=1025;
			//mw->w[0]=0;
			//mw->w[1]=74;	
			mw->run_a();
			calculate_snr();
			write_line_stats(0,"2_lines_near_far_mw");
			write_line_stats(1,"2_lines_near_far_mw");
			print_summary();
			delete_network();
			break;
		}

		case MULTIUSER_NEW: {
			multiuser_new *mn;
			mn = new multiuser_new;
			mn->rate_targets_off=true;
			//mn->b_target[0]=872;
			//mn->b_target[1]=1025;
			//mw->w[0]=0;
			//mw->w[1]=74;	
			mn->graph_loading=true;
			mn->run_a();
			calculate_snr();
			write_line_stats(0,"2_lines_near_far_mn");
			write_line_stats(1,"2_lines_near_far_mn");
			print_summary();
			delete_network();
			break;
		}
		/*
		case MULTIUSER_FAIRER: {
			multiuser_fairer *mf;
			mf = new multiuser_fairer;
			mf->b_target[1]=rate1_target*250;
			mf->run();
			calculate_snr();
			write_line_stats(0,"2_lines_near_far_mf");
			write_line_stats(1,"2_lines_near_far_mf");
			print_summary();
			delete_network();
			break;
		}*/
		
	}
}

void single_qos_2_lines_near_far_1(int scenario)
{

	struct line** line_array;
	lines = read_network("2lines_near_far1.txt");
	int max_rate1,max_rate0;
	line_array = new struct line*[lines];
	double rate1_target=2;

	channel_matrix = calc_channel_matrix(lines);	
	
	write_channel_matrix("2lines_near_far_1");

	set_background_noise(-140);
	//set_alien_xtalk_model(MODEL_A);
	set_alien_xtalk_model(NONE);

	print_network(lines);


	for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
		line_array[i]->set_line_single_qos(9.95);
	}

	switch(scenario) {
		case SINGLE_IWF: {
			calculate_snr();
			
			for (int i=0;i<5;i++) {
				am_load_fm(1,rate1_target*250,0.110,BOTH_HALVES);
				calculate_snr();
				am_load_ra(0,0.110,BOTH_HALVES);
				calculate_snr();
			}

			write_line_stats(0,"2_lines_near_far_1_IWF");	
			write_line_stats(1,"2_lines_near_far_1_IWF");	
			print_summary();
			delete_network();		
			break;
		}
		case OSB: {
			char tag[100];
			struct osb_2_params* p;
			p = new osb_2_params;
			calculate_snr();
			p->osb_2_params_init();
			p->rate1_target=rate1_target*250;
			p->e=1;
			osb_2(p);
			calculate_snr();
			sprintf(tag,"2_lines_near_far_1_OSB_%d_%d",p->rate0_target,p->rate1_target);
			write_line_stats(0,tag);
			write_line_stats(1,tag);
			print_summary();
			delete_network();
			break;
		}
		case RATE_REGION_IWF: {
			int j=0;
			calculate_snr();

			am_load_ra(1,0.110,BOTH_HALVES);

			calculate_snr();
			max_rate1=line_array[1]->rate[0];
	
			printf("max rate on RT line is %d = approx %d Mbps\n",max_rate1,max_rate1/250);

			printf("0\t%lf\n",(double)(max_rate1*FRAMERATE)/1e6);
		
			for(int rate1=max_rate1/250;rate1>=0;rate1--) {
				//hjunnnnnnnnnnnnnnnnnnnnnnnnnnnnnnp0lkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkko
				//comment by rambo
				do {
					am_load_fm(1,rate1*250,0.110,BOTH_HALVES);	
					calculate_snr();
					am_load_ra(0,0.110,BOTH_HALVES);
					calculate_snr();
					j++;
				//} while(check_all_margins(0.01));
				} while(j<5);
				printf("%lf\t%lf\n",(double)(line_array[0]->rate[0]*FRAMERATE/1e6),
							(double)(line_array[1]->rate[0]*FRAMERATE/1e6));

			}
	
			break;
		}
		case RATE_REGION_OSB: {
			struct osb_2_params* p;
			p = new osb_2_params;
			calculate_snr();	
			am_load_ra(0,0.110,BOTH_HALVES);			
			calculate_snr();
			max_rate0=line_array[0]->rate[0];

			printf("%lf\t0\n",(double)(max_rate0*FRAMERATE/1e6));

			for (int rate0=max_rate0/250;rate0>=0;rate0--) {
				p->osb_2_params_init();
				p->rate0_target=rate0*250;
				osb_2(p);
				calculate_snr();
				printf("%lf\t%lf\n",(double)(line_array[0]->rate[0]*FRAMERATE/1e6),
                                                        (double)(line_array[1]->rate[0]*FRAMERATE/1e6));
			}
				
			break;
		}
		case MULTIUSER_GREEDY: {
			multiuser_greedy *mg;
			mg = new multiuser_greedy;	
			mg->b_target[1]=rate1_target*250;
			mg->run();
			calculate_snr();
			write_line_stats(0,"2_lines_near_far_1_mg");
			write_line_stats(1,"2_lines_near_far_1_mg");
			print_summary();
			delete_network();
			break;
		}
		case MULTIUSER_WEIGHTED: {
			multiuser_weighted *mw;
			mw = new multiuser_weighted;
			mw->b_target[1]=rate1_target*250;
			mw->run_a();
			calculate_snr();
			write_line_stats(0,"2_lines_near_far_1_mw");
			write_line_stats(1,"2_lines_near_far_1_mw");
			print_summary();
			delete_network();
			break;
		}
		
	}
}

void single_qos_2_lines_near_far_2(int scenario)
{

	struct line** line_array;
	lines = read_network("2lines_near_far2.txt");
	int max_rate1,max_rate0;
	line_array = new struct line*[lines];
	double rate1_target=10.5;

	channel_matrix = calc_channel_matrix(lines);	
	
	write_channel_matrix("2lines_near_far_2");

	set_background_noise(-140);
	//set_alien_xtalk_model(MODEL_A);
	set_alien_xtalk_model(NONE);

	print_network(lines);


	for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
		line_array[i]->set_line_single_qos(9.95);
	}

	switch(scenario) {
		case SINGLE_IWF: {
			calculate_snr();
			
			for (int i=0;i<5;i++) {
				am_load_fm(1,rate1_target*250,0.110,BOTH_HALVES);
				calculate_snr();
				am_load_ra(0,0.110,BOTH_HALVES);
				calculate_snr();
			}

			write_line_stats(0,"2_lines_near_far_2_IWF");	
			write_line_stats(1,"2_lines_near_far_2_IWF");	
			print_summary();
			delete_network();		
			break;
		}
		case OSB: {
			char tag[100];
			struct osb_2_params* p;
			p = new osb_2_params;
			calculate_snr();
			p->osb_2_params_init();
			p->rate1_target=rate1_target*250;
			p->e=5;
			osb_2(p);
			calculate_snr();
			sprintf(tag,"2_lines_near_far_2_OSB_%d_%d",p->rate0_target,p->rate1_target);
			write_line_stats(0,tag);
			write_line_stats(1,tag);
			print_summary();
			delete_network();
			break;
		}
		case RATE_REGION_IWF: {
			int j=0;
			calculate_snr();

			am_load_ra(1,0.110,BOTH_HALVES);

			calculate_snr();
			max_rate1=line_array[1]->rate[0];
	
			printf("max rate on RT line is %d = approx %d Mbps\n",max_rate1,max_rate1/250);

			printf("0\t%lf\n",(double)(max_rate1*FRAMERATE)/1e6);
		
			for(int rate1=max_rate1/250;rate1>=0;rate1--) {
				//hjunnnnnnnnnnnnnnnnnnnnnnnnnnnnnnp0lkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkko
				//comment by rambo
				do {
					am_load_fm(1,rate1*250,0.110,BOTH_HALVES);	
					calculate_snr();
					am_load_ra(0,0.110,BOTH_HALVES);
					calculate_snr();
					j++;
				//} while(check_all_margins(0.01));
				} while(j<5);
				printf("%lf\t%lf\n",(double)(line_array[0]->rate[0]*FRAMERATE/1e6),
							(double)(line_array[1]->rate[0]*FRAMERATE/1e6));

			}
	
			break;
		}
		case RATE_REGION_OSB: {
			struct osb_2_params* p;
			p = new osb_2_params;
			calculate_snr();	
			am_load_ra(0,0.110,BOTH_HALVES);			
			calculate_snr();
			max_rate0=line_array[0]->rate[0];

			printf("%lf\t0\n",(double)(max_rate0*FRAMERATE/1e6));

			for (int rate0=max_rate0/250;rate0>=0;rate0--) {
				p->osb_2_params_init();
				p->rate0_target=rate0*250;
				osb_2(p);
				calculate_snr();
				printf("%lf\t%lf\n",(double)(line_array[0]->rate[0]*FRAMERATE/1e6),
                                                        (double)(line_array[1]->rate[0]*FRAMERATE/1e6));
			}
				
			break;
		}
		case MULTIUSER_GREEDY: {
			multiuser_greedy *mg;
			mg = new multiuser_greedy;	
			mg->b_target[1]=rate1_target*250;
			mg->run();
			calculate_snr();
			write_line_stats(0,"2_lines_near_far_2_mg");
			write_line_stats(1,"2_lines_near_far_2_mg");
			print_summary();
			delete_network();
			break;
		}
		case MULTIUSER_WEIGHTED: {
			multiuser_weighted *mw;
			mw = new multiuser_weighted;
			mw->b_target[1]=rate1_target*250;
			mw->run_a();
			calculate_snr();
			write_line_stats(0,"2_lines_near_far_2_mw");
			write_line_stats(1,"2_lines_near_far_2_mw");
			print_summary();
			delete_network();
			break;
		}
		
	}
}
void single_qos_2_lines_CO_short_long(int scenario)
{

	struct line** line_array;
	lines = read_network("2lines_CO_short_long.txt");
	int max_rate1,max_rate0;
	line_array = new struct line*[lines];
	double rate1_target=8.8;

	channel_matrix = calc_channel_matrix(lines);	

	write_channel_matrix("2lines_short_long");

	set_background_noise(-140);
	//set_alien_xtalk_model(MODEL_A);
	set_alien_xtalk_model(NONE);

	print_network(lines);


	for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
		line_array[i]->set_line_single_qos(9.95);
	}

	switch(scenario) {
		case SINGLE_IWF: {
			calculate_snr();
			
			for (int i=0;i<5;i++) {
				am_load_fm(1,(int)rate1_target*250,0.110,BOTH_HALVES);
				calculate_snr();
				am_load_ra(0,0.110,BOTH_HALVES);
				calculate_snr();
			}

			write_line_stats(0,"2_lines_CO_short_long_IWF");	
			write_line_stats(1,"2_lines_CO_short_long_IWF");	
			print_summary();
			delete_network();		
			break;
		}
		case OSB: {
			char tag[100];
			struct osb_2_params* p;
			p = new osb_2_params;
			calculate_snr();
			p->osb_2_params_init();
			p->rate1_target=(int)(rate1_target*250);
			p->e=5;
			osb_2(p);
			calculate_snr();
			sprintf(tag,"2_lines_CO_short_long_OSB_%d_%d",p->rate0_target,p->rate1_target);
			write_line_stats(0,tag);
			write_line_stats(1,tag);

			print_summary();
			delete_network();
			break;
		}

		case ISB: {
			//char tag[100];	
			isb* i;
			i = new isb;
			calculate_snr();
			//i->rate_target[1]=810;
			i->e=1;
			i->run();
			delete i;
			calculate_snr();
			print_summary();
			delete_network();
			break;
		}

		case RATE_REGION_IWF: {
			int j=0;
			calculate_snr();

			am_load_ra(1,0.110,BOTH_HALVES);

			calculate_snr();
			max_rate1=line_array[1]->rate[0];
	
			printf("max rate on RT line is %d = approx %d Mbps\n",max_rate1,max_rate1/250);

			printf("0\t%lf\n",(double)(max_rate1*FRAMERATE)/1e6);
		
			for(int rate1=max_rate1/250;rate1>=0;rate1--) {
				//hjunnnnnnnnnnnnnnnnnnnnnnnnnnnnnnp0lkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkko
				//comment by rambo
				do {
					am_load_fm(1,rate1*250,0.110,BOTH_HALVES);	
					calculate_snr();
					am_load_ra(0,0.110,BOTH_HALVES);
					calculate_snr();
					j++;
				//} while(check_all_margins(0.01));
				} while(j<5);
				printf("%lf\t%lf\n",(double)(line_array[0]->rate[0]*FRAMERATE/1e6),
							(double)(line_array[1]->rate[0]*FRAMERATE/1e6));

			}
	
			break;
		}
		case RATE_REGION_OSB: {
			struct osb_2_params* p;
			p = new osb_2_params;
			calculate_snr();	
			am_load_ra(0,0.110,BOTH_HALVES);			
			calculate_snr();
			max_rate0=line_array[0]->rate[0];

			printf("%lf\t0\n",(double)(max_rate0*FRAMERATE/1e6));

			for (int rate0=max_rate0/250;rate0>=0;rate0--) {
				p->osb_2_params_init();
				p->rate0_target=rate0*250;
				osb_2(p);
				calculate_snr();
				printf("%lf\t%lf\n",(double)(line_array[0]->rate[0]*FRAMERATE/1e6),
                                                        (double)(line_array[1]->rate[0]*FRAMERATE/1e6));
			}
				
			break;
		}
		case MULTIUSER_GREEDY: {
			multiuser_greedy *mg;
			mg = new multiuser_greedy;	
			mg->b_target[1]=(int)rate1_target*250;
			mg->run();
			calculate_snr();
			write_line_stats(0,"2_lines_CO_short_long_mg");
			write_line_stats(1,"2_lines_CO_short_long_mg");
			print_summary();
			delete_network();
			break;
		}
		case MULTIUSER_WEIGHTED: {
			multiuser_weighted *mw;
			mw = new multiuser_weighted;
			mw->b_target[1]=(int)(rate1_target*250);
			mw->run_a();
			delete mw;
			calculate_snr();
			write_line_stats(0,"2_lines_CO_short_long_mw");
			write_line_stats(1,"2_lines_CO_short_long_mw");
			print_summary();
			delete_network();
			break;
		}

		case MULTIUSER_NEW: {
			multiuser_new *mn;
			mn = new multiuser_new;
			mn->rate_targets_off=true;
			//mn->b_target[0]=872;
			//mn->b_target[1]=1025;
			//mw->w[0]=0;
			//mw->w[1]=74;	
			mn->run_a();
			calculate_snr();
			write_line_stats(0,"2_lines_CO_short_long_mn");
			write_line_stats(1,"2_lines_CO_short_long_mn");
			print_summary();
			delete_network();
			break;
		}
		
	}
}


void single_qos_50_equal_lines(int scenario)
{
	//struct line** line_array;
        //lines = read_network("50_lines.txt");
        lines = read_network("10_lines.txt");
	line_array = new struct line*[lines];

        channel_matrix = calc_channel_matrix(lines);

        set_background_noise(-140);
        set_alien_xtalk_model(MODEL_A);
        //set_alien_xtalk_model(NONE);

        print_network(lines);

	//printf("%d lines read\n",lines);

        for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
                line_array[i]->set_line_single_qos(9.95);
        }

	switch(scenario) {
		case MULTIUSER_GREEDY: {
			multiuser_greedy *mg;
                        mg = new multiuser_greedy;
                        mg->run();
			print_summary();
			break;
		}
		case MULTIUSER_WEIGHTED: {
			multiuser_weighted *mw;
			mw = new multiuser_weighted;
			mw->run_a();
			delete mw;
			calculate_snr();
			print_summary();
			break;
		}
		case ISB: {
			isb *i;
			i = new isb;
			calculate_snr();
			i->run();
			delete i;
			print_summary();
			delete_network();
		}
			
	}


}

void multiuser_greedy_vs_iwf(int scenario)
{

	struct line** line_array;
	lines = read_network("multiuser_greedy_vs_iwf.txt");
	line_array = new struct line*[lines];
	int rate0_target=4;

	channel_matrix = calc_channel_matrix(lines);	
	write_channel_matrix(NULL);

	set_background_noise(-140);
	//set_alien_xtalk_model(MODEL_A);
	set_alien_xtalk_model(NONE);

	print_network(lines);

	for (int i=0;i<lines;i++) {
		line_array[i]=get_line(i);
		line_array[i]->set_line_single_qos(9.95);
	}

	switch(scenario) {
		case SINGLE_IWF: {
			calculate_snr();
			int iters=0;
			do{
				am_load_fm(0,rate0_target*250+20,0.110,BOTH_HALVES);			
				calculate_snr();
				am_load_ra(1,0.110,BOTH_HALVES);
				calculate_snr();
				iters++;
			
			}while(iters<10);
			

			do {
                                for (int i=0;i<lines;i++) {
					margin_tighten(i);
                                	calculate_snr();
				}
                        } while(check_all_margins(2));
			

			write_line_stats(0,"multiuser_greedy_vs_IWF_IWF");	
			write_line_stats(1,"multiuser_greedy_vs_IWF_IWF");	
			print_summary();
			delete_network();		
			break;
		}
		case MULTIUSER_GREEDY: {
			multiuser_greedy *mg;
                        mg = new multiuser_greedy;
			calculate_snr();
			mg->b_target[0]=rate0_target*250;
			mg->run();
			calculate_snr();
			write_line_stats(0,"multiuser_greedy_vs_IWF_mg");	
			write_line_stats(1,"multiuser_greedy_vs_IWF_mg");	
			print_summary();
			delete_network();
			break;
		}
	}


}
#ifdef NOTDEF

/* 11 lines total
   10 ADSL operating in fixed-margin mode loaded iteratively (3km)
   1 dual qos line loaded with dbw_cbr (3km)
   
*/

void scenario1()
{

	int i,j,iters=50;

	lines = read_network("scenario1.txt");

        //list_print(list_head);
        print_network(lines);
        //print_xtalkers(lines);

        channel_matrix = calc_channel_matrix(lines);    // returns pointer to malloced memory with chan_matrix
	
	for (i=0;i<lines;i++) {
                set_psd_on_line(i,0);           // turn off line i
        }

	calculate_snr();

	
	for (i=0;i<lines;i++) {
                set_psd_on_line(i,0);           // turn off line i
        }
	
	i=0;
	do {
		//if (i % 2 == 0 || i == 0) {
		if (1) {
			printf("Loading from line 0 up\n");
			dbw_cbr(0,9.95,6.05,400,400);
			calculate_snr();
			for (j=1;j<lines;j++) {
				am_load_fm(j,800,0);
				calculate_snr();
			}
		}
		else {
			printf("loading from line %d down\n",lines-1);
			for (j=lines-1;j>=1;j--) {
				am_load_fm(j,800,0);
				calculate_snr();
			}
			dbw_cbr(0,9.95,6.05,400,400);
			calculate_snr();
		}
		i++;
		//printf("%d\n",i);
	}while(check_all_margins(1));
	
	printf("Converged after %d iterations\n",i);
	
	write_line_stats(1,"s1");
	write_line_stats_dual(0,"s1");
	write_line_stats(10,"s1");


}



void scenario2()
{


	int i,j;
        lines = read_network("scenario2.txt");

        //list_print(list_head);
        print_network(lines);
        //print_xtalkers(lines);

        channel_matrix = calc_channel_matrix(lines);    // returns pointer to malloced memory with chan_matrix

        for (i=0;i<lines;i++) {
                set_psd_on_line(i,0);           // turn off line i
        }

        calculate_snr();

	dbw_cbr(0,9.95,6.05,400,400);
	calculate_snr();

	check_all_margins(0.1);	

	write_line_stats_dual(0,"s2");

}	


/* two lines only used to test osb versus iwf */

void scenario3()
{

	int i=0,j;
        lines = read_network("scenario3.txt");

        //list_print(list_head);
        print_network(lines);
        //print_xtalkers(lines);

        channel_matrix = calc_channel_matrix(lines);    // returns pointer to malloced memory with chan_matrix

	print_channel_matrix(1);

	calculate_snr();
/*
	do {
		printf("iteration %d\n",i++);
		am_load_ra(0,0.1,0);
		calculate_snr();
		am_load_fm(1,461,0);
		calculate_snr();	
	}while(i<10);
	//}while(check_all_margins(2));
	
	write_line_stats(0,"s3iwf");
	write_line_stats(1,"s3iwf");
*/
	//find_rate_region();
	//osb_find_rate_region();

	//osb_2(0.5);
	//osb_2_new();
	isb_2();
/*
	do {
		//am_load_fm(1,989,0);
		am_load_fm(0,1495,0);
		calculate_snr();
		am_load_fm(1,989,0);
		calculate_snr();
	}while(check_all_margins(0.1));
		
*/
	//calculate_snr();

	//print2_line_stats(0,1);

	//check_all_margins(0.2);

//calculate_snr();

	//write_line_stats(0,"s3osb");
	//write_line_stats(1,"s3osb");

	
}

/* a 3km line and a 2km line trying to test osb_dual_qos*/

void scenario4()
{

	int i=0,j;
	double p1,p2;
        FILE *fp;
	//bbl_mg_dual *a;
	//multiuser_greedy *a;
	weighted_multiuser_greedy *a;
	struct line* current;
	//lines = read_network("scenario4.txt");
	lines = read_network("scenario5.txt");

        //list_print(list_head);
        print_network(lines);
        
	channel_matrix = calc_channel_matrix(lines);    // returns pointer to malloced memory with chan_matrix
        
	//print_xtalkers(lines);


	print_channel_matrix(1);

	//write_channel_matrix();
	//exit(0);


	set_background_noise(-140);

	set_alien_xtalk_model(MODEL_A);

	//calculate_snr();

	for (i=0;i<lines;i++) {
		current=get_line(i);
		current->set_line_single_qos(9.95);
	}


	//a = new bbl_mg_dual;
	//a = new multiuser_greedy;
	a = new weighted_multiuser_greedy;
	
	a->run();

	return;
/*
	do {
		printf("iteration %d\n",i++);
		am_load_ra(0,0.1,0);
		calculate_snr();
		am_load_fm(1,461,0);
		calculate_snr();	
	}while(i<10);
	//}while(check_all_margins(2));
	
	write_line_stats(0,"s3iwf");
	write_line_stats(1,"s3iwf");
*/
	//find_rate_region();
	//osb_find_rate_region();

//	fp=fopen("dbw_rate_region.txt","w");
	
//	if (fp==NULL)
//		exit(2);

//	for (p1=0.0;p1<=0.1001;p1+=0.005) {
//        	for (p2=0.0;p2<=0.005;p2+=0.00005) {

			i=0;
			do {
				i++;
				//dbw_vbr_cbr(0,9.95,6.02,200);
				//am_load_ra(0,0.1,0);
				//am_load_fm(0,1740,0);
				//dbw_vbr_cbr(1,9.95,6.02,400);
				//dbw_cbr(0,9.95,6.02,1583,200);
				dbw_vbr_cbr_exhaustive(0,9.95,6.02,100,0.1);
				//dbw_vbr_cbr_exhaustive(0,9.95,6.02,600);
				calculate_snr();
				//dbw_vbr_cbr_exhaustive(1,9.95,6.02,100,p2);
				//calculate_snr();
				//print_line_stats_dual(0);
				//dbw_vbr_cbr(1,9.95,6.02,200);
				dbw_cbr_exhaustive(1,9.95,6.02,1250,100);	// service0 at 5mbps
				//am_load_ra(1,0.1,0);
				calculate_snr();
			//} while(check_all_margins(2));
			} while(i<10);
			

//			fprintf(fp,"p1=%lf\tp2=%lf\t0,0=%d\t0,1=%d\t1,0=%d\t1,1=%d\n",p1,p2,linerate(0,0),linerate(0,1),linerate(1,0),linerate(1,1));

//			if ((linerate(0,1) == 100) && (linerate(1,1) == 100)) {
//				printf("%d\t%d\n",linerate(0,0),linerate(1,0));
//			}
	
/*			
			do {
				margin_tighten(0);
				calculate_snr();
				margin_tighten(1);
				calculate_snr();
			} while(check_all_margins(0.1));
*/
//		}
//	}

//	fclose(fp);


	//osb_2_dual();

	//write_line_stats_dual(0,"s4osb_RT_5mbps");
	//write_line_stats_dual(1,"s4osb_RT_5mbps");

	//osb_2();
	//asb_2();

	
	//calculate_snr();
	
	//write_line_stats_dual(0,"s4osbdual_no_alien_xtalk");		
	//write_line_stats_dual(1,"s4osbdual_no_alien_xtalk");		

	//write_line_stats(0,"s4isbdual_as_single_qos");		
	//write_line_stats(1,"s4isbdual_as_single_qos");		

/*
	do {
		//am_load_fm(1,989,0);
		am_load_fm(0,1495,0);
		calculate_snr();
		am_load_fm(1,989,0);
		calculate_snr();
	}while(check_all_margins(0.1));
		
*/
	//calculate_snr();

	//print2_line_stats(0,1);

	//check_all_margins(0.2);

	//calculate_snr();

	//write_line_stats(0,"s3osb");
	//write_line_stats(1,"s3osb");

	
}



void scenario5()
{

        int i=0,j;
        lines = read_network("scenario5.txt");

        //list_print(list_head);
        print_network(lines);
        print_xtalkers(lines);

        channel_matrix = calc_channel_matrix(lines);    // returns pointer to malloced memory with chan_matrix

        print_channel_matrix(1);

	write_channel_matrix();

	set_background_noise(-140);

	set_alien_xtalk_model(MODEL_A);
       
	//osb_2();
 
	calculate_snr();

	find_rate_region();
	osb_find_rate_region();
	//am_load_ra(0,0.110,0);
	exit(0);

	do {
		am_load_ra(0,0.110,0);
		calculate_snr();
		am_load_ra(1,0.110,0);
		calculate_snr();
	} while(check_all_margins(0.5));
		//i++;
	//} while(i<10);

	do {
		margin_tighten(0);
		calculate_snr();
		margin_tighten(1);
		calculate_snr();
	} while(check_all_margins(0.1));

	//print_line_stats(0);


}	

void scenario7()
{

	lines = read_network("scenario5.txt");

        //list_print(list_head);
        print_network(lines);
        //print_xtalkers(lines);

        channel_matrix = calc_channel_matrix(lines);    // returns pointer to malloced memory with chan_matrix

        print_channel_matrix(1);

        calculate_snr();

	//single_user_genetic();
	//gsl_test();
}

#endif

#endif
