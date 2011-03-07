#ifndef PSD_VECTOR_H
#define PSD_VECTOR_H

#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <Eigen/Core>
#include <Eigen/LU>
USING_PART_OF_NAMESPACE_EIGEN

struct cache_entry{
	cache_entry()
	{
		hash=0;
		p_vec=NULL;
		b_vec=NULL;
		next=NULL;
	}
	unsigned int hash;
	double *p_vec;
	int *b_vec;
	struct cache_entry *next;
};

class psd_vector {
	
	public:
	psd_vector();
	~psd_vector();
	int calc(int *,double *,int,double*);
	void print_cache_size()
	{
		printf("current cache size is %lf bytes\n",_total_cache_size*(double)sizeof(cache_entry)/1024);
	}
	bool _caching_on;
	
	private:
	int _cache_hits;
	int _cache_misses;
	double _compares;
	double _total_cache_size;
	int _calls;
	double *_cache_size;
	double _collisions;	
	
	struct cache_entry **_cache_head;
	struct cache_entry **_cache_tail;
	
	struct cache_entry * _cache;

	pthread_mutex_t _lock;


	struct line **_line_array;	

#ifdef GSL
	gsl_matrix *A;
	gsl_vector *B;
	gsl_vector *X;
	gsl_matrix *invA;
	gsl_permutation *p;
	gsl_vector *tau;
#endif

#ifdef EIGEN
	MatrixXd A;
	VectorXd B;
	VectorXd x;
#endif

#ifdef MKL


#endif

	int calculate_psd_vector(int *,double *,int,double *);
	/*double F(double g,int b)
	{
        	return pow(10,(g+3)/10)*(pow(2,b)-1);
        }*/
	
	unsigned int murmurhash(const unsigned char *,int,unsigned int);
	unsigned int murmurhash2(const void *, int len, unsigned int);
	
	int insert_cache_entry(int *,int,double *,unsigned int);

	

};

#endif
