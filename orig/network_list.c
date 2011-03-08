#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "multiuser_load.h"
#include <sys/time.h>
#include <time.h>
struct line *list_head = NULL;
bool read_n;
int n;
int read_network(const char * network_file) {
	FILE *fp;
        char buf[255];
	char fn[255];
	char num[10];
	int i,j,nt,lt,lines=0;
	sprintf(fn,"%s/%s",ROOT_DIR,network_file);
        if ((fp = fopen(fn,"r")) == NULL) {
                printf("Cant open network file %s\n",fn);
                exit(1);
        }
	while (1) {
		fgets(buf,(sizeof(buf) - 1),fp);
		if (buf[0] == '#') {	// comment in scenario file
			if (feof(fp)) 
				break;
			continue;
		}
		if(feof(fp))
			break;
		for (i=0,j=0;i<10;i++,j++) {
			if (buf[i] == '.') {
				num[j] = '\0';
				i++;
				break;
			}
			else {
				num[j] = buf[i];
			}
		}
		lt = (double)atoi(num);
		for (j=0;(buf[i] != '\n');i++,j++) {
			num[j] = buf[i];
		}
		num[j] = '\0';
		nt = (double)atoi(num);	
		add_line(lt,nt);
		lines++;
		if (read_n) {
			if (lines == n)
				break;
		}
	}
	return lines;
}
void print_network(int lines)
{
	int i;
	struct line *current = list_head;
	for (i=0;i<lines;i++) {
		while (current->line_id != i)
			current = current->next;
		print_line(current);
	}
}
void print_network_to_file(int lines,FILE *fp)
{
	int i;
	struct line *current = list_head;
	for (i=0;i<lines;i++) {
		while (current->line_id != i)
			current = current->next;
		print_line_to_file(current,fp);
		fprintf(fp,"\n\n");
	}
}
void print_xtalkers(int lines)
{
	int i,j;
	for (i=0;i<lines;i++) {
        	for (j=0;j<lines;j++) {
                	if (i != j)	// victim line_id is j, crosstalker line_id is i
                        	print_case(i,j); 
                        }
        }
}
void print_case(int x, int v)
{
 	struct line* xtalker = list_head;
        struct line* victim = list_head;
	xtalker=get_line(x);
	victim=get_line(v);
        printf("xtalker id = %d\tlt=%.0lf nt=%.0lf\n",x,xtalker->lt,xtalker->nt);
        printf("victim id = %d\tlt=%.0lf nt=%.0lf\n",v,victim->lt,victim->nt);
        print_line(victim);
        print_line(xtalker);
	printf("upstream: case %d\n",get_case(xtalker,victim,UPSTREAM));
	printf("downstream: case %d\n",get_case(xtalker,victim,DOWNSTREAM));
	printf("channel 0 xtalk gain = %4.2lf\n",10*log10(get_xtalk_gain(x,v,0)));
        printf("\n");
}
void print_line(struct line* line) {
	int metres=0;
	printf("%d: ",line->line_id);
        for (metres=0;line->lt-metres>0;metres+=200) {
                printf("-");
        }
        printf("|");
        for (metres=0;line->length-metres-200>0;metres+=200) {
                printf("-");
        }
        printf("|");
	printf("\tLength = %.0lf metres\n",line->length);
}
void print_line_to_file(struct line* line, FILE *fp) {
	int metres=0;
	fprintf(fp,"%d: ",line->line_id);
        for (metres=0;line->lt-metres>0;metres+=200) {
                fprintf(fp,"-");
        }
        fprintf(fp,"\\textbar");
        for (metres=0;line->length-metres-200>0;metres+=200) {
                fprintf(fp,"-");
        }
        fprintf(fp,"\\textbar");
	fprintf(fp,"\tLength = %.0lf metres\n",line->length);
}
void add_line(double lt, double nt)
{
	static int i=0;
	struct line* current;
	if (list_head == NULL) {
		i=0;
		/*
		if ((list_head = (struct line*)malloc(sizeof(struct line))) == NULL) {
			printf("cannot malloc line!\n");
			exit(1);		
		}
		*/
		list_head = new line;
		list_head->lt = lt;
		list_head->nt = nt;
		list_head->length = nt - lt;
		list_head->line_id = i;
		list_head->alien_xtalk_model = NOT_SET;
		list_head->background_noise = NOT_SET;
		list_head->next = NULL;
		i++;
		return;
	}
	current = list_head;
	while (current->next != NULL)
		current = current->next;	// find end of list	
	/*
	if ((current->next = (struct line*)malloc(sizeof(struct line))) == NULL) {
		printf("cannot malloc line!\n");
		exit(1);
	}
	*/
	current->next = new line;
	current = current->next;
	current->lt = lt;
        current->nt = nt;
        current->length = nt - lt;
	current->line_id = i;
	current->alien_xtalk_model = NOT_SET;
	current->background_noise = NOT_SET;
        current->next = NULL;
	i++;
	return;
}
void delete_network()
{
	struct line* current = list_head;
	struct line* next;
	while (current != NULL) {
		next = current->next;
		//free(current);
		delete current;
		current=next;
	}	
	list_head=NULL;	
	free(channel_matrix);
}
struct line *get_line(int id)
{
	struct line *current = list_head;
	while (current->line_id != id) {
		current=current->next;	
		if (current == NULL)
			return NULL;
	}
	return current;
}
void list_print(struct line* current)
{
	while (current != NULL) {
		printf("head = %p line[%d] = %lf %lf next = %p\n",current,current->line_id,current->lt,current->nt,current->next);
		current = current->next;
	}
}
void write_all_line_stats(const char * tag)
{
	for (int user=0;user<lines;user++) {
		write_line_stats(user,tag);
	}
}
void write_all_p_stats(const char *tag)
{
	char fn[255];
	FILE *fp;
	if (tag != NULL)
		sprintf(fn,"%s/data/p_stats%s.txt",ROOT_DIR,tag);
	else
		sprintf(fn,"%s/data/p_stats.txt",ROOT_DIR);		
	if ((fp = fopen(fn,"w")) == NULL) {
		printf("cannot open file to write\n");
		exit(1);
	}
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%g,",dbmhz_to_watts(line_array[user]->psd[tone]));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}
void write_all_b_stats(const char *tag)
{
	char fn[255];
	FILE *fp;
	if (tag != NULL)
		sprintf(fn,"%s/data/b_stats%s.txt",ROOT_DIR,tag);
	else
		sprintf(fn,"%s/data/b_stats.txt",ROOT_DIR);		
	if ((fp = fopen(fn,"w")) == NULL) {
		printf("cannot open file to write\n");
		exit(1);
	}
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int user=0;user<lines;user++) {
			fprintf(fp,"%d,",line_array[user]->b[tone]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}
void write_all_b_and_p_stats(const char *tag)
{
	char fn[300];
	FILE *fp;
	if (tag != NULL)
		sprintf(fn,"%s/data/b_and_p_stats_%s.txt",ROOT_DIR,tag);
	else
		sprintf(fn,"%s/data/b_and_p_stats.txt",ROOT_DIR);		
	if ((fp = fopen(fn,"w")) == NULL) {
		printf("cannot open file to write\n");
		exit(1);
	}
	fprintf(fp,"#%d\n",lines);
	double b_tot=0.0,p_tot=0.0;
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		double b_ave=0.0,b_var=0.0;
		double p_ave=0.0,p_var=0.0;
		//fprintf(fp,"%d ",tone);
#ifdef VDSL_UPSTREAM
		if (tone == LOWERCHANNELS) {
			for (int times=0;times<2;times++) {
				if (times == 0)
					fprintf(fp,"%6.4lf ",freq_matrix[tone-1]+CHANNEL_BANDWIDTH);
				else 
					fprintf(fp,"%6.4lf ",freq_matrix[tone]-CHANNEL_BANDWIDTH);
				for (int user=0;user<lines;user++) {
					fprintf(fp,"%d %s ",0,"-inf");
				}
				fprintf(fp,"\n");
			}
		}
#endif
		fprintf(fp,"%6.4lf ",freq_matrix[tone]);
		for (int user=0;user<lines;user++) {
			b_ave+=line_array[user]->b[tone];
			b_tot+=line_array[user]->b[tone];
			p_ave+=dbmhz_to_watts(line_array[user]->psd[tone]);
			p_tot+=dbmhz_to_watts(line_array[user]->psd[tone]);
			if (line_array[user]->is_frac != true) {
				fprintf(fp,"%d %4.2lf ",line_array[user]->b[tone]
							,line_array[user]->psd[tone]);
			}
			else {
				fprintf(fp,"%6.4lf %4.2lf ",line_array[user]->_b[tone]
							,line_array[user]->psd[tone]);
			}
		}
		b_ave/=lines;
		p_ave/=lines;
		fprintf(fp,"%4.2lf %6.4g ",b_ave,p_ave*1e3);
		for (int user=0;user<lines;user++) {
			b_var+=pow(line_array[user]->b[tone]-b_ave,2);
			p_var+=pow((dbmhz_to_watts(line_array[user]->psd[tone])-p_ave)*1e3,2);
		}
		b_var/=lines;
		p_var/=lines;
		fprintf(fp,"%4.2lf %6.4g %6.4g %4.0lf %6.4g ",b_var,p_var,b_ave/(p_ave*1e3),b_tot,p_tot*1e3);
		if (tone < lines) {
			fprintf(fp,"%d %d",line_array[tone]->b_total,line_array[tone]->no_fext_b_total);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}
void write_line_stats(int line_id,const char *tag)
{
	int i;
	struct line* c;
	char fn[300];
	FILE *fp;
	c = get_line(line_id);
	if (tag != NULL)
		sprintf(fn,"%s/data/line%d_stats_%s.txt",ROOT_DIR,line_id,tag);
	else
		sprintf(fn,"%s/data/line%d_stats.txt",ROOT_DIR,line_id);		
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
	fprintf(fp,"#channels = %d\n",DMTCHANNELS);
	fprintf(fp,"#p_total = %lf mW\n",c->p_total*1e3);
	fprintf(fp,"#b_total = %d  data_rate = %d approx\n",c->b_total,c->b_total*FRAMERATE);
	fprintf(fp,"#background noise = %4.2lf dBm/Hz\n",c->background_noise);
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
void write_line_stats_dual(int line_id,const char *tag)
{
        int i;
        struct line* c;
        char fn[255];
        FILE *fp;
        c = get_line(line_id);
        if (tag != NULL)
                sprintf(fn,"%s/data/line%d_dualstats%s.txt",ROOT_DIR,line_id,tag);
        else
                sprintf(fn,"%s/data/line%d_dualstats.txt",ROOT_DIR,line_id);
        if ((fp = fopen(fn,"w")) == NULL) {
                printf("cannot open file to write\n");
                exit(1);
        }
	if (c->is_dual == 0) {
		printf("this is a single qos line and you are trying to print dual qos stats\n");
		exit(1);
	}
        fprintf(fp,"#Line %d: Service0 rate = %d Service1 rate = %d\n",c->line_id,c->rate[0],c->rate[1]);
	fprintf(fp,"#Service0 gamma = %4.2lf Service1 gamma = %4.2lf\n",c->gamma[0],c->gamma[1]); 
	fprintf(fp,"#length = %lf lt = %lf nt = %lf\n",c->length,c->lt,c->nt);
        fprintf(fp,"#p_total = %lf mW\n",c->p_total*1e3);
        fprintf(fp,"#b_total = %d  data_rate = %d approx\n",c->b_total,c->b_total*FRAMERATE);
	fprintf(fp,"#background noise = %4.2lf dBm/Hz\n",c->background_noise);
        fprintf(fp,"#tone\tbits\tpsd\t\tsnr\t\t\tcnr\t\t\tgain\t\t\tmargin\tp_e\t\ts0\ts1\n");
        for (i=0;i<DMTCHANNELS;i++) {
                fprintf(fp,"%d\t%d\t",i,c->b[i]);
		fprintf(fp,"%4.2lf\t\t%4.2e  %4.2lf dB\t%4.2e  %4.2lf dB\t%4.2e  %5.3lf dB\t%3.1lf\t%4.2e\t",
                                                                        c->psd[i],c->snr[i],10*log10(c->snr[i]),
                                                                        c->cnr[i],10*log10(c->cnr[i]),c->gain[i],10*log10(c->gain[i]),
                                                                        c->gamma_m[i],c->symerr[i]);
        	if (c->service[i] == 0)
                        fprintf(fp,"%d\t%d\n",c->b[i],0);
                else if (c->service[i] == 1)
                        fprintf(fp,"%d\t%d\n",0,c->b[i]);
                else {
                        printf("ouch!");
                        exit(1);
                }
	}
        fclose(fp);
}
void print_line_stats(int line_id)
{
        int i;
        struct line* c;
        c = get_line(line_id);
	/*
        if (tag != NULL)
                sprintf(fn,"%s/data/line%d_stats%s.txt",ROOT_DIR,line_id,tag);
        else
                sprintf(fn,"%s/data/line%d_stats.txt",ROOT_DIR,line_id);
        if ((fp = fopen(fn,"w")) == NULL) {
                printf("cannot open file to write\n");
                exit(1);
        }
	*/
	printf("#length = %lf lt = %lf nt = %lf\n",c->length,c->lt,c->nt);
        printf("#p_total = %lf mW\n",c->p_total*1e3);
        printf("#b_total = %d  data_rate = %d approx\n",c->b_total,c->b_total*FRAMERATE);
        printf("#tone\tbits\tpsd\t\tsnr\t\t\tcnr\t\t\tgain\t\t\tmargin\tp_e\n");
        for (i=0;i<DMTCHANNELS;i++) {
                if (c->b[i])
			printf("%d\t%d\t%4.2lf\t\t%4.2e  %4.2lf dB\t%4.2e  %4.2lf dB\t%4.2e  %5.3lf dB\t%3.1lf\t%4.2e\n",
                                                                        i,c->b[i],c->psd[i],c->snr[i],10*log10(c->snr[i]),
                                                                        c->cnr[i],10*log10(c->cnr[i]),c->gain[i],10*log10(c->gain[i]),
                                                                        c->gamma_m[i],c->symerr[i]);
        }
        //fclose(fp);
}
void print2_line_stats(int line1, int line2)
{
	int k;
	struct line* l1;
	struct line* l2;
	l1=get_line(line1);
	l2=get_line(line2);
	printf("\tline%d\tline%d\n",line1,line2);
	printf("p_tot\t%4.2lf\t%4.2lf\n",l1->p_total*1e3,l2->p_total*1e3);
	printf("b_tot\t%d\t%d\n",l1->b_total,l2->b_total);
	for (k=0;k<DMTCHANNELS;k++) {
		printf("%d\t%d\t%d\t%4.2lf\t%4.2lf\t%4.2lf\t%4.2lf\n",k,l1->b[k],l2->b[k],l1->psd[k],l2->psd[k],l1->gamma_m[k],l2->gamma_m[k]);
	}
}
void print_line_stats_dual(int line_id)
{
	int i;
        struct line* c;
        c = get_line(line_id);
        printf("#length = %lf lt = %lf nt = %lf\n",c->length,c->lt,c->nt);
        printf("#p_total = %lf mW\n",c->p_total*1e3);
        printf("#b_total = %d  data_rate = %d approx\n",c->b_total,c->b_total*FRAMERATE);
        printf("#tone\tbits\tpsd\t\tsnr\t\t\tcnr\t\t\tgain\t\t\tmargin\tp_e\t\ts0\ts1\n");
        for (i=0;i<DMTCHANNELS;i++) {
                printf("%d\t%d\t",i,c->b[i]);
                printf("%4.2lf\t\t%4.2e  %4.2lf dB\t%4.2e  %4.2lf dB\t%4.2e  %5.3lf dB\t%3.1lf\t%4.2e\t",
                                                                        c->psd[i],c->snr[i],10*log10(c->snr[i]),
                                                                        c->cnr[i],10*log10(c->cnr[i]),c->gain[i],10*log10(c->gain[i]),
                                                                        c->gamma_m[i],c->symerr[i]);
                if (c->service[i] == 0)
                        printf("%d\t%d\n",c->b[i],0);
                else if (c->service[i] == 1)
                        printf("%d\t%d\n",0,c->b[i]);
                else {
                        printf("ouch!");
                        exit(1);
                }
        }
}
void print_channel_matrix(int tones)
{
	int i,j,k;
	for (k=0;k<tones;k++) {
                for (j=0;j<lines;j++) {
                        for (i=0;i<lines;i++) {
                                printf("%6.4lf  ",10*log10(*(channel_matrix + k + (DMTCHANNELS*j) + (DMTCHANNELS * lines * i))));
                        }
                printf("\n");
                }
                printf("\n");
        }
}
void write_channel_matrix(const char * tag)
{
	return;
	int i,j,k;
	FILE *fp;
	char fn[300];
	sprintf(fn,"%s/data/channel_matrix_%s.txt",ROOT_DIR,tag);
	fp=fopen(fn,"w");
	if (fp==NULL)
		exit(2);
	fprintf(fp,"#\t");
	for (j=0;j<lines;j++) {
		for (i=0;i<lines;i++) {
			fprintf(fp,"%d%d\t\t\t",j,i);
        	}
        }
	fprintf(fp,"\n");
	for (k=0;k<DMTCHANNELS;k++) {
		fprintf(fp,"%d\t",k);
                for (j=0;j<lines;j++) {
                        for (i=0;i<lines;i++) {
                                //fprintf(fp,"%4.2e %4.2lf\t\t",get_xtalk_gain(j,i,k),10*log10(*(channel_matrix + k + (DMTCHANNELS*j) + (DMTCHANNELS * lines * i))));
                                fprintf(fp,"%4.2e \t",get_xtalk_gain(j,i,k));
                                //fprintf(fp,"%e\t",*(channel_matrix + k + (DMTCHANNELS*j) + (DMTCHANNELS * lines * i)));
                        }
                }
                fprintf(fp,"\n");
        }
	fclose(fp);
}
void print_summary()
{
	int i,j;
	struct line *current;
	int total_cap=0;
	double _total_cap=0;
	printf("\n\nSummary:\n");
	for (i=0;i<lines;i++) {
		current=get_line(i);
		printf("Line %d:\n",current->line_id);
		printf("Loading algorithm = %s\n",current->loading_algo);
		//printf("Background noise = %4.2lf\n",current->background_noise);
		printf("Alien xtalk model = %d\n",current->alien_xtalk_model);
		printf("p_total = %lf mW\n",current->p_total*1e3);
		for (j=0;j<8;j++) {
			if (current->rate[j] != 0) {
				if (!current->is_frac) {
					printf("service%d rate = %d = %lf Mbps gamma = %4.2lf\n",j,current->rate[j],
											   (double)(current->rate[j]*FRAMERATE/1e6),current->gamma[j]);
					total_cap+=current->rate[j];
				}
				else {
					printf("service%d rate = %lf = %lf Mbps gamma = %4.2lf\n",j,current->_rate[j],
											   (double)(current->_rate[j]*FRAMERATE/1e6),current->gamma[j]);
					_total_cap+=current->_rate[j];
				}
			}
			else 
				continue;
		}
		if (!current->is_frac)
			printf("no fext rate = %d = %lf Mbps\n",current->no_fext_b_total,(double)(current->no_fext_b_total*FRAMERATE/1e6));
		else 
			printf("no fext rate = %lf = %lf Mbps\n",current->_no_fext_b_total,(double)(current->_no_fext_b_total*FRAMERATE/1e6));
		//printf("\n");
	}
	if (!current->is_frac) {
		printf("Total bundle capacity = %d = %lf Mbps\n",total_cap,total_cap*FRAMERATE/1e6);
		printf("Average = %lf\n",total_cap*FRAMERATE/(1e6*lines));
	}
	else {
		printf("Total bundle capacity = %lf = %lf Mbps\n",_total_cap,_total_cap*FRAMERATE/1e6);
		printf("Average = %lf\n",_total_cap*FRAMERATE/(1e6*lines));
	}
}
void create_results_summary_file(char *tag, char *other_info)
{
#ifndef XC
	int i,j;
	struct line *current;
	int total_cap=0;
	double _total_cap=0;
	FILE * fp;
	char temp_fn[400];	
	char pdflatex_buf[400];
	char cp_buf[400];
	srand((unsigned int)time(NULL));
	int key = rand() % 10000;
	//sprintf(temp_fn,"/tmp/%s_%d.tex",tag,key);
	sprintf(temp_fn,"/tmp/%s.tex",tag);
	fp = fopen(temp_fn,"w");
	if (fp == NULL) {
		printf("Cannot open tex results file\n");
		exit(2);
	}
	fprintf(fp,"%s","\\documentclass{article}\n");
	fprintf(fp,"%s","\\usepackage{graphicx}\n");
	fprintf(fp,"%s","\\usepackage{fullpage}\n");
	fprintf(fp,"%s","\\begin{document}\n");
	print_network_to_file(lines,fp);
	fprintf(fp,"\n\nSummary:\n\n");
	for (i=0;i<lines;i++) {
		current=get_line(i);
		fprintf(fp,"Line %d:\n\n",current->line_id);
		fprintf(fp,"Loading algorithm = %s\n\n",current->loading_algo);
		//printf("Background noise = %4.2lf\n",current->background_noise);
		fprintf(fp,"Alien xtalk model = %d\n\n",current->alien_xtalk_model);
		fprintf(fp,"p total = %lf mW\n\n",current->p_total*1e3);
		for (j=0;j<8;j++) {
			if (current->rate[j] != 0) {
				if (!current->is_frac) {
					fprintf(fp,"service%d rate = %d = %lf Mbps gamma = %4.2lf\n\n",j,current->rate[j],
											   (double)(current->rate[j]*FRAMERATE/1e6),current->gamma[j]);
					total_cap+=current->rate[j];
				}
				else {
					fprintf(fp,"service%d rate = %lf = %lf Mbps gamma = %4.2lf\n\n",j,current->_rate[j],
											   (double)(current->rate[j]*FRAMERATE/1e6),current->gamma[j]);
					_total_cap+=current->_rate[j];
				}
			}
			else 
				continue;
		}
		fprintf(fp,"no fext rate = %d = %lf Mbps\n\n",current->no_fext_b_total,(double)(current->no_fext_b_total*FRAMERATE/1e6));
		//printf("\n");
	}
	if (!current->is_frac) 
		fprintf(fp,"Total bundle capacity = %d\n\n\n",total_cap);
	else
		fprintf(fp,"Total bundle capacity = %lf\n\n\n",_total_cap);
	fprintf(fp,"DMTCHANNELS = %d\n\n\n",DMTCHANNELS);
	if (other_info != NULL)
		fprintf(fp,"%s\n\n",other_info);
	fprintf(fp,"\\includegraphics[width=\\textwidth]{%s/graphs/b_and_p_stats_%s.png}\n\n",ROOT_DIR,tag);
	//fprintf(fp,"\\includegraphics[width=\\textwidth]{%s/graphs/channel_matrix_%s.png}\n\n",ROOT_DIR,tag);
	fprintf(fp,"%s","\\end{document}\n\n");
	fclose(fp);
	printf("Running pdflatex...\n");
	sprintf(pdflatex_buf,"pdflatex -output-directory=/tmp %s",temp_fn);
	if (system(pdflatex_buf) == 0)
		printf("Ok\n");
	printf("Copying pdf file\n");
	//sprintf(cp_buf,"cp /tmp/%s_%d.pdf %s/results/",tag,key,ROOT_DIR);	
	sprintf(cp_buf,"cp /tmp/%s.pdf %s/results/",tag,ROOT_DIR);	
	if (system(cp_buf) == 0)
		printf("Ok\n");
	printf("Copying to most recent file\n");
	//sprintf(cp_buf,"cp /tmp/%s_%d.pdf %s/results/most_recent.pdf",tag,key,ROOT_DIR);	
	sprintf(cp_buf,"cp /tmp/%s.pdf %s/results/most_recent.pdf",tag,ROOT_DIR);	
	if (system(cp_buf) == 0)
		printf("Ok\n");
#endif
}
int linerate(int line_id,int service)
{
	struct line* current;
	current=get_line(line_id);
	return current->rate[service];
}
void generate_graphs(char *tag) 
{
#ifndef XC
	char prog_string1[400];
	char prog_string2[400];
	sprintf(prog_string1,"../scripts/graph-channel-matrix.sh -p ../data/channel_matrix_%s.txt",tag);
	sprintf(prog_string2,"../scripts/graph-all-bits.sh -p ../data/b_and_p_stats_%s.txt",tag);
	/*printf("Generating channel matrix graph....\n");
	if (system(prog_string1) == 0)
		printf("Ok\n");
	*/	
	printf("Generating bits/psd graph....\n");
	if (system(prog_string2) == 0)
		printf("Ok\n");
#endif
}
