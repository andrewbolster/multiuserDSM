#include "multiuser_load.h"
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>


double nearest_integer_bit_p(int tone, int user, double p)
{

	double p_up,p_down;
	double init_b;

	init_b = no_fext_bits(tone,user,p);

	//printf("init b = %lf\n",init_b);

	p_down=single_psd(tone,user,(int)init_b);
	p_up=single_psd(tone,user,(int)init_b+1);
	
	if ((p-p_down) < (p_up-p)) {	// difference in p(floor(b)) is less than difference in p(floor(b+1)) 
		return p_down;
	}
	else {
		return p_up;
	}
	
}

double single_psd(int tone,int user,int b)
{
	double p;
	p = (UNdB(12.95)*(pow(2,(double)b)-1)*minus_140_dbmhz_watts)/get_channel_gain(user,tone);
	return p;
}

double no_fext_bits(int tone,int user,double p)
{
	return log2(1+(p*get_channel_gain(user,tone))/(minus_140_dbmhz_watts*UNdB(12.95)));
}

double nxt(int xtalker,int victim,double p,int tone)
{

	return 10*log10(p*get_xtalk_gain(xtalker,victim,tone)/get_channel_gain(victim,tone));

}

bool vector_next(int *vec,int len,int max,bool reset)
{
        int i;

        if (reset) {
                for (i=0;i<len;i++) {
                        vec[i]=0;
                }
                return true;
        }


        for (i=len-1;i>=0;i--) {
                if (vec[i]+1 > max) {
                        vec[i]=0;
                }
                else {
                        vec[i]++;
                        return true;
                }
        }

        return false;

}

void print_vector(int *vec,const char * tag)
{
	printf("%s = [ ",tag);
                for (int user=0;user<lines;user++) {
                        printf("%d ",vec[user]);
        }
     	printf("]\n");

}

void print_vector(double *vec,const char * tag)
{
	printf("%s = [ ",tag);
                for (int user=0;user<lines;user++) {
                        printf("%12.10e ",vec[user]);
        }
     	printf("]\n");

}

void print_vector_hires(double *vec,const char * tag)
{
	printf("%s = [ ",tag);
                for (int user=0;user<lines;user++) {
                        printf("%50.48e ",vec[user]);
        }
     	printf("]\n");

}

void print_vector_to_file(double *vec,const char * tag,FILE *fp) 
{
	fprintf(fp,"%s = [",tag);
                for (int user=0;user<lines;user++) {
                        fprintf(fp,"%4.2e ",vec[user]);
        }
     	fprintf(fp,"]\n");

}

void print_vector_psd(double *vec,const char * tag)
{
	printf("%s = [ ",tag);
                for (int user=0;user<lines;user++) {
                        printf("%.20lf ",watts_to_dbmhz(vec[user]));
        }
     	printf("]\n");

}



int calculate_b_vector_one_bit_set(double *p,double *g,int channel,int line_id,double *_b)             // line_id is the line where the bit is set, the other lines have p set at -40
{

        double noise=0;
	struct line *current = get_line(line_id);
        // first calculate noise into line line_id

        for (int user=0;user<lines;user++) {
                if (user!=line_id)
                        noise+=p[user]*get_xtalk_gain(user,line_id,channel);
        }
        noise+=minus_140_dbmhz_watts+current->alien_xtalk[channel];

        // now calculate p on line line_id required to support _b[line_id]

        p[line_id]=pow(10,(g[line_id]+3)/10)*(pow(2,_b[line_id])-1)*noise/get_channel_gain(line_id,channel);

	//printf("p required to support b[%d]=%lf is %4.2lf\n",line_id,_b[line_id],watts_to_dbmhz(p[line_id])); 

        for (int user=0;user<lines;user++) {
                noise=minus_140_dbmhz_watts+current->alien_xtalk[channel];	
                if (user!=line_id) {
                        for (int user1=0;user1<lines;user1++) {	// calculate noise into line user
                                if (user1==user)
                                        ;
                                else
                                        noise+=p[user1]*get_xtalk_gain(user1,user,channel);
                        }
			// calculate b on line user
          //              printf("b[%d] = %4.2lf\n",user, log2( 1+p[user]*get_channel_gain(user,channel)/(pow(10,(g[user]+3)/10)*noise) ));
			_b[user]=log2( 1+p[user]*get_channel_gain(user,channel)/(pow(10,(g[user]+3)/10)*noise) );
                }
		else {
	//		printf("b[%d] = %lf\n",user,_b[user]);
		}
        }

	return 0;

}


int calculate_b_vector_one_bit_set(double *p,double g,int channel,int line_id,double *_b)             // line_id is the line where the bit is set, the other lines have p set at -40
{

        double noise=0;
	struct line *current = get_line(line_id);
        // first calculate noise into line line_id

        for (int user=0;user<lines;user++) {
                if (user!=line_id)
                        noise+=p[user]*get_xtalk_gain(user,line_id,channel);
        }
        noise+=minus_140_dbmhz_watts+current->alien_xtalk[channel];

        // now calculate p on line line_id required to support _b[line_id]

        p[line_id]=pow(10,(g+3)/10)*(pow(2,_b[line_id])-1)*noise/get_channel_gain(line_id,channel);

	//printf("p required to support b[%d]=%lf is %4.2lf\n",line_id,_b[line_id],watts_to_dbmhz(p[line_id])); 

        for (int user=0;user<lines;user++) {
                noise=minus_140_dbmhz_watts+current->alien_xtalk[channel];	
                if (user!=line_id) {
                        for (int user1=0;user1<lines;user1++) {	// calculate noise into line user
                                if (user1==user)
                                        ;
                                else
                                        noise+=p[user1]*get_xtalk_gain(user1,user,channel);
                        }
			// calculate b on line user
          //              printf("b[%d] = %4.2lf\n",user, log2( 1+p[user]*get_channel_gain(user,channel)/(pow(10,(g[user]+3)/10)*noise) ));
			_b[user]=log2( 1+p[user]*get_channel_gain(user,channel)/(pow(10,(g+3)/10)*noise) );
                }
		else {
	//		printf("b[%d] = %lf\n",user,_b[user]);
		}
        }

	return 0;

}




int calculate_b_vector_from_psd(double *p,double *g,int channel,double *_b)
{

	//struct line **line_array;

	//line_array = new struct line* [lines];

	//for (int user=0;user<lines;user++) {
	//	line_array[user] = get_line(user);
	//}

	for (int user=0;user<lines;user++) {
		double noise=minus_140_dbmhz_watts+line_array[user]->alien_xtalk[channel];
		
		for (int user1=0;user1<lines;user1++) {
			if (user1!=user)
				noise+=p[user1]*get_xtalk_gain(user1,user,channel);
		}
		_b[user]=log2(1+p[user]*get_channel_gain(user,channel)/(pow(10,g[user]/10)*noise));
	}
	
	//delete[] line_array;

	return 0;
}

int calculate_b_vector_from_psd(double *p,double g,int channel,double *_b)
{

	//struct line **line_array;

	//line_array = new struct line* [lines];

	//for (int user=0;user<lines;user++) {
	//	line_array[user] = get_line(user);
	//}

	for (int user=0;user<lines;user++) {
		double noise=minus_140_dbmhz_watts+line_array[user]->alien_xtalk[channel];
		
		for (int user1=0;user1<lines;user1++) {
			if (user1!=user)
				noise+=p[user1]*get_xtalk_gain(user1,user,channel);
		}
		_b[user]=log2(1+p[user]*get_channel_gain(user,channel)/(pow(10,g/10)*noise));
	}
	
	//delete[] line_array;

	return 0;
}

int calculate_b_vector_flat_psd(double p,double *g,int channel,int *_b)
{

	for (int user=0;user<lines;user++) {
		double noise=dbmhz_to_watts(-140)+alien_xtalk(user,channel);
		
		for (int user1=0;user1<lines;user1++) {
			if (user1!=user)
				noise+=p*get_xtalk_gain(user1,user,channel);
		}
		
		_b[user]=log2(1+p*get_channel_gain(user,channel)/(pow(10,(g[user]+3)/10)*noise));	
		if (_b[user] > MAXBITSPERTONE) 
			_b[user] = MAXBITSPERTONE;
	}


	return 0;
}

int calculate_b_vector_flat_psd(double p,double *g,int channel,double *_b)
{

	for (int user=0;user<lines;user++) {
		double noise=dbmhz_to_watts(-140)+alien_xtalk(user,channel);
		
		for (int user1=0;user1<lines;user1++) {
			if (user1!=user)
				noise+=p*get_xtalk_gain(user1,user,channel);
		}
		
		_b[user]=log2(1+p*get_channel_gain(user,channel)/(pow(10,g[user]/10)*noise));	
		if (_b[user] > MAXBITSPERTONE) 
			_b[user] = MAXBITSPERTONE;
	}


	return 0;
}

void print_matrix(gsl_matrix * m,const char * str)
{
        int i,j;
        int rows = m->size1;
        int cols = m->size2;

        printf("%s\n",str);
        for (i=0;i<cols;i++) {
                for (j=0;j<rows;j++) {
                        printf("%g\t",gsl_matrix_get(m,i,j));
                }
                printf("\n");
        }


}

void print_vec(gsl_vector *v,const char *str)
{
        int i;
        int length=v->size;

        printf("%s\n",str);
        for (i=0;i<length;i++) {
                printf("%g\n",gsl_vector_get(v,i));
        }

}

bool any_neg(double *v)
{
	for (int user=0;user<lines;user++) {
		if (v[user] < 0) {
			return true;
		}
	}

	return false;
}

void dump_vector(double * vector,const char *tag)
{
	char fn[300];
        FILE *fp;

        sprintf(fn,"/tmp/vector-%s.txt",tag);

        fp = fopen(fn,"w");

        if (fp==NULL)
                exit(2);

        for (int tone=0;tone<DMTCHANNELS;tone++) {
#ifdef VDSL_UPSTREAM
                if (tone == LOWERCHANNELS) {
                        for (int times=0;times<2;times++) {
                                if (times == 0) {
                                        fprintf(fp,"%6.4lf\t",freq_matrix[tone-1]+CHANNEL_BANDWIDTH);
                                }
                                else {
                                        fprintf(fp,"%6.4lf\t",freq_matrix[tone]-CHANNEL_BANDWIDTH);
                                }

                                fprintf(fp,"%d ",0);
                                fprintf(fp,"\n");
                        }
                }
#endif
                fprintf(fp,"%6.4lf\t",freq_matrix[tone]);
                fprintf(fp,"%6.4g ",vector[tone]);
                fprintf(fp,"\n");
        }

        fclose(fp);

}

void dump_int_vector(int *vector,const char *tag)
{
	char fn[300];
        FILE *fp;

        sprintf(fn,"/tmp/vector-%s.txt",tag);

        fp = fopen(fn,"w");

        if (fp==NULL)
                exit(2);

        for (int tone=0;tone<DMTCHANNELS;tone++) {
#ifdef VDSL_UPSTREAM
                if (tone == LOWERCHANNELS) {
                        for (int times=0;times<2;times++) {
                                if (times == 0) {
                                        fprintf(fp,"%6.4lf\t",freq_matrix[tone-1]+CHANNEL_BANDWIDTH);
                                }
                                else {
                                        fprintf(fp,"%6.4lf\t",freq_matrix[tone]-CHANNEL_BANDWIDTH);
                                }

                                fprintf(fp,"%d ",0);
                                fprintf(fp,"\n");
                        }
                }
#endif
                fprintf(fp,"%6.4lf\t",freq_matrix[tone]);
                fprintf(fp,"%d ",vector[tone]);
                fprintf(fp,"\n");
        }

        fclose(fp);

}


