#ifndef PSD_CACHE_H
#define PSD_CACHE_H

#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include "psd_vector.h"

USING_PART_OF_NAMESPACE_EIGEN


struct cache_ent{
	cache_ent()
	{
		hash=0;
		p_vec=NULL;
		b_vec=NULL;
		next=NULL;
	}
	unsigned int hash;
	float *p_vec;
	char *b_vec;
	struct cache_entry *next;
};


class psd_cache {
	
	public:
	psd_cache();
	~psd_cache();
	int check(int *,double *,int,double*);
	void print_cache_size()
	{
		printf("current cache size is %lf bytes\n",_total_cache_size*(double)sizeof(cache_ent)/1024);
	}
	bool _caching_on;
	
	int insert_cache_entry(int *,int,double *);
	int replace_cache_entry(int *,int,double *);
	
	pthread_mutex_t *_tone_locks;
	
	private:
	double _cache_hits;
	double _cache_misses;
	double _compares;
	double _total_cache_size;
	int _calls;
	double *_cache_size;
	double _collisions;	
	
	struct cache_ent **_cache_head;
	struct cache_ent **_cache_tail;
	
	struct cache_ent * _cache;

	pthread_mutex_t _lock;


	struct line **_line_array;	

	unsigned int murmurhash(const unsigned char *,int,unsigned int);
	unsigned int murmurhash2(const void *, int len, unsigned int);
	

};

int calculate_psd_vector(int *b,double *gamma,int channel,double *_p,psd_cache *cache);

#endif
