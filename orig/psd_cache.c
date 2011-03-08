#include "multiuser_load.h"
#include "psd_cache.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <float.h>
//#include <mkl_lapack.h>
#include <Eigen/Core>
#include <Eigen/LU>
USING_PART_OF_NAMESPACE_EIGEN
#define inv 1
#define F(g,b) pow(10,(g+3)/10)*(pow(2,b)-1)
#define TONE_CACHE_SIZE 80000
enum {
	CACHE_HIT,
	CACHE_MISS,
	CACHE_COLLISION
};
psd_cache::psd_cache()
{
	_cache_hits=0;
	_cache_misses=0;
	_total_cache_size=0;
	_calls=0;
	_compares=0;
	_collisions=0;
	_cache_head = new struct cache_ent *[DMTCHANNELS];
	_cache_tail = new struct cache_ent *[DMTCHANNELS];
	_cache_size = new double [DMTCHANNELS];
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		_cache_size[tone]=0;
		_cache_head[tone]=new struct cache_ent [TONE_CACHE_SIZE];	
	}
	_line_array = new struct line *[lines];
	for (int user=0;user<lines;user++) {
		_line_array[user]=get_line(user);
	} 	
	_tone_locks = new pthread_mutex_t[DMTCHANNELS];
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		pthread_mutex_init(&(_tone_locks[tone]),NULL);
	}
}
psd_cache::~psd_cache()
{
	/*
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		struct cache_ent *current=_cache_head[tone];
		struct cache_ent *next;
		while(current !=NULL) {
			next=current->next;
			delete[] current->p_vec;
			delete current;
			current=next;	
		}
	}
	delete[] _cache_head;
	delete[] _cache_tail;
	*/
	printf("This psd_cache class managed %lf cache hits against %lf cache misses and %lf collisions\n",_cache_hits,_cache_misses,_collisions);	
	printf("calls = %d\n",_calls);
	printf("total cache size = %.0lf entries\n",_total_cache_size);
	printf("total cache in bytes = %lf kbytes\n",_total_cache_size*(double)sizeof(cache_ent)/1024);
	/*
	for (int tone=0;tone<DMTCHANNELS;tone++) {
		printf("Cache size on channel %d is %.0lf\n",tone,_cache_size[tone]);
	}
	*/
	//printf("total compares = %lf\n",_compares);
	//printf("total collisions = %lf\n",_collisions);
}
int psd_cache::check(int *b,double *gamma,int channel,double *p)
{
	//struct cache_entry* cache=_cache_head[channel];
	//bool debug=false;
	//double scale_fac = pow(2,-15);
	//pthread_mutex_lock(&_lock);
	//pthread_mutex_lock(&(_tone_locks[channel]));
	unsigned int hash=0;
	/*
	if (b[0] == 15 && b[1] == 12 && channel == 0) {
		debug=true; 
		printf("cache head on channel 0 is %p\n",cache);
		getchar();
	}
	*/
	//hash = murmurhash((unsigned char *)b,lines*sizeof(int),hash);
	//
	int b_hash[lines];
        for (int user=0;user<lines;user++) {
                b_hash[user] = (int)b[user];
        }
	hash = murmurhash2((void *)b_hash,lines*sizeof(int),hash);
//		int index = (int)((double)hash*scale_fac);
	int index = hash % TONE_CACHE_SIZE;
	if (_cache_head[channel][index].p_vec == NULL) {
		/*
		psd_cache::calculate_psd_cache(b,gamma,channel,p);
		_calls++;
		_cache_misses++;
		_cache_head[channel][index].p_vec = new double [lines];
		_cache_head[channel][index].b_vec = new int [lines];
		_total_cache_size++;
		_cache_size[channel]++;
		_cache_head[channel][index].hash=hash;
		memcpy(_cache_head[channel][index].p_vec,p,sizeof(double)*lines);
		memcpy(_cache_head[channel][index].b_vec,b,sizeof(int)*lines);
		pthread_mutex_unlock(&_lock);
		return 0;
		*/
		//pthread_mutex_unlock(&(_tone_locks[channel]));
		_cache_misses++;
		return CACHE_MISS;
	}
	else {
		if (hash != _cache_head[channel][index].hash) {	// different hashes, index collision after mod operation
			//printf("Collision!\n");
			//psd_cache::calculate_psd_cache(b,gamma,channel,p);
			_collisions++;
			//_cache_misses++;
			//pthread_mutex_unlock(&_lock);
			//pthread_mutex_unlock(&(_tone_locks[channel]));
			return CACHE_COLLISION;
		}
		for (int user=0;user<lines;user++) {
			if (b[user] != _cache_head[channel][index].b_vec[user]) { // proper hash collision 
				// eg (8 8 8 8 12 0 0 0) == (6 7 9 8 13 7 13 13) 
				//psd_cache::calculate_psd_cache(b,gamma,channel,p);
				//_collisions++;
				//pthread_mutex_unlock(&_lock);
				//pthread_mutex_unlock(&(_tone_locks[channel]));			
				_collisions++;
				//_cache_misses++;
				return CACHE_COLLISION;
			}
		}
		//memcpy(p,_cache_head[channel][index].p_vec,sizeof(double)*lines);
		_cache_hits++;
		for (int user=0;user<lines;user++) {
                        p[user] = (double)_cache_head[channel][index].p_vec[user];
                }
		//pthread_mutex_unlock(&_lock);
		//pthread_mutex_unlock(&(_tone_locks[channel]));
		return CACHE_HIT;
	}
/*
	if (debug) {
		printf("hash = %u\n",hash);
	}
*/
	/*	
	if (_caching_on) {
		while (cache != NULL) {			// right now, assume gammas are all the same!	
			if (cache->hash == hash) {
				memcpy(p,cache->p_vec,sizeof(double)*lines);
				_cache_hits++;
				return 1;		// cache hit
			}
			_compares++;
			cache=cache->next;
		}
	}
	psd_cache::calculate_psd_cache(b,gamma,channel,p);
	_calls++;
	if (_caching_on) {
		insert_cache_entry(b,channel,p,hash);		// thats right! assuming gamma is the same for all entries!
		_cache_misses++;
	}	
	return 0;
	*/
}
int psd_cache::insert_cache_entry(int *b,int channel,double *p)
{
	int b_hash[lines];
        for (int user=0;user<lines;user++) {
                b_hash[user] = (int)b[user];
        }
	unsigned int hash = 0;
	hash = murmurhash2((void *)b,lines*sizeof(int),hash);
	int index = hash % TONE_CACHE_SIZE;
	//_cache_misses++;
	_cache_head[channel][index].p_vec = new float [lines];
	_cache_head[channel][index].b_vec = new char [lines];
	_total_cache_size++;
	_cache_size[channel]++;
	_cache_head[channel][index].hash=hash;
	//memcpy(_cache_head[channel][index].p_vec,p,sizeof(float)*lines);
	//memcpy(_cache_head[channel][index].b_vec,b,sizeof(char)*lines);
	for (int user=0;user<lines;user++) {
                _cache_head[channel][index].p_vec[user] = (float)p[user];
                _cache_head[channel][index].b_vec[user] = (char)b[user];
        }
	return 0;
}
int psd_cache::replace_cache_entry(int *b,int channel,double *p)
{
	int b_hash[lines];
        for (int user=0;user<lines;user++) {
                b_hash[user] = (int)b[user];
        }
	unsigned int hash = 0;
	hash = murmurhash2((void *)b,lines*sizeof(int),hash);
	int index = hash % TONE_CACHE_SIZE;
	//_cache_misses++;
	//_cache_head[channel][index].p_vec = new float [lines];
	//_cache_head[channel][index].b_vec = new char [lines];
	//_total_cache_size++;
	//_cache_size[channel]++;
	_cache_head[channel][index].hash=hash;
	//memcpy(_cache_head[channel][index].p_vec,p,sizeof(float)*lines);
	//memcpy(_cache_head[channel][index].b_vec,b,sizeof(char)*lines);
	for (int user=0;user<lines;user++) {
                _cache_head[channel][index].p_vec[user] = (float)p[user];
                _cache_head[channel][index].b_vec[user] = (char)b[user];
        }
	return 0;
}
unsigned int psd_cache::murmurhash(const unsigned char *data,int len,unsigned int h)
{
	const unsigned int m = 0x7fd652ad;
	const int r = 16;
	h += 0xdeadbeef;
	while(len >= 4)
	{
		h += *(unsigned int *)data;
		h *= m;
		h ^= h >> r;
		data += 4;
		len -= 4;
	}
	switch(len)
	{
	case 3:
		h += data[2] << 16;
	case 2:
		h += data[1] << 8;
	case 1:
		h += data[0];
		h *= m;
		h ^= h >> r;
	};
	h *= m;
	h ^= h >> 10;
	h *= m;
	h ^= h >> 17;
	return h;
}
unsigned int psd_cache::murmurhash2( const void * key, int len, unsigned int seed )
{
	// 'm' and 'r' are mixing constants generated offline.
	// They're not really 'magic', they just happen to work well.
	const unsigned int m = 0x5bd1e995;
	const int r = 24;
	// Initialize the hash to a 'random' value
	unsigned int h = seed ^ len;
	// Mix 4 bytes at a time into the hash
	const unsigned char * data = (const unsigned char *)key;
	while(len >= 4)
	{
		unsigned int k = *(unsigned int *)data;
		k *= m; 
		k ^= k >> r; 
		k *= m; 
		h *= m; 
		h ^= k;
		data += 4;
		len -= 4;
	}
	// Handle the last few bytes of the input array
	switch(len)
	{
	case 3: h ^= data[2] << 16;
	case 2: h ^= data[1] << 8;
	case 1: h ^= data[0];
	        h *= m;
	};
	// Do a few final mixes of the hash to ensure the last few
	// bytes are well-incorporated.
	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;
	return h;
}
int calculate_psd_vector(int *b,double *gamma,int channel,double *_p,psd_cache *cache)
{
	if (cache == NULL) {
		calculate_psd_vector(b,gamma,channel,_p);
		return 0;
	}
	pthread_mutex_lock(&(cache->_tone_locks[channel]));
	int ret = cache->check(b,gamma,channel,_p);
	switch(ret) {
		case CACHE_HIT:
			pthread_mutex_unlock(&(cache->_tone_locks[channel]));
			return CACHE_HIT;
			break;
		case CACHE_MISS:
		case CACHE_COLLISION:
			calculate_psd_vector(b,gamma,channel,_p);
			/*
			int i,j,sig;
			struct line* current;
			gsl_matrix *A = gsl_matrix_calloc(lines,lines);
			gsl_vector *B = gsl_vector_calloc(lines);
			gsl_vector *X = gsl_vector_calloc(lines);
			gsl_matrix *invA = gsl_matrix_calloc(lines,lines);
			gsl_permutation *p = gsl_permutation_calloc(lines);
			for (i=0;i<lines;i++) {
				for (j=0;j<lines;j++) {
					if (i==j)
						gsl_matrix_set(A,i,j,1);
					else {
						gsl_matrix_set(A,i,j,-1*F(gamma[i],b[i])*get_xtalk_gain(j,i,channel)/get_xtalk_gain(i,i,channel));		
					}
				}
			}
			for(i=0;i<lines;i++) {
				//gsl_vector_set(B,i,F(gamma[i],b[i])*(dbmhz_to_watts(-140)+alien_xtalk(i,channel))/get_xtalk_gain(i,i,channel));
				// forget about alien_xtalk for now!!
				//gsl_vector_set(B,i,F(gamma[i],b[i])*(dbmhz_to_watts(-140)+alien_xtalk(i,channel))/get_xtalk_gain(i,i,channel));
				gsl_vector_set(B,i,F(gamma[i],b[i])*(dbmhz_to_watts(-140))/get_xtalk_gain(i,i,channel));
			}
		#ifdef QR
			gsl_linalg_QR_decomp(A,tau);
			gsl_linalg_QR_solve(A,tau,B,X);
		#endif
		#ifdef inv
			if (gsl_linalg_LU_decomp(A,p,&sig)) {
				printf("LU decomp failed!");
				exit(1);
			}
			if (gsl_linalg_LU_invert(A,p,invA)) {
				printf("LU invert failed");
				exit(1);
			}
			if (gsl_blas_dgemv(CblasNoTrans,1.0,invA,B,0.0,X)) {           // ans = inv(D)*y
				printf("failed!");
				exit(1);
			}
		#endif
			for (i=0;i<lines;i++) {
				_p[i]=gsl_vector_get(X,i);
				if (_p[i] > 0 && _p[i] < 4.31e-14)
					_p[i]=0;
				if (_p[i] < 0 && _p[i] > -4.31e-14)
					_p[i]=0;
			}
			gsl_matrix_free(A);
			gsl_matrix_free(invA);
			gsl_vector_free(X);
			gsl_vector_free(B);
			gsl_permutation_free(p);	
			*/
			break;
	}
	if (ret == CACHE_MISS) {
		cache->insert_cache_entry(b,channel,_p);
		pthread_mutex_unlock(&(cache->_tone_locks[channel]));
		return CACHE_MISS;
	}	
	else {
		cache->replace_cache_entry(b,channel,_p);
		pthread_mutex_unlock(&(cache->_tone_locks[channel]));
		return CACHE_COLLISION;
	}
}
