#include "multiuser_load.h"
#include "psd_vector.h"
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
#define GSL //This was bugging me
#define F(g,b) pow(10,(g+3)/10)*(pow(2,b)-1)
#define TONE_CACHE_SIZE 80000
bool psd_caching;
psd_vector::psd_vector()
{
	_cache_hits=0;
	_cache_misses=0;
	_total_cache_size=0;
	_calls=0;
	_compares=0;
	_collisions=0;
	_caching_on=psd_caching;
	if (_caching_on) {	
		_cache_head = new struct cache_entry *[DMTCHANNELS];
		_cache_tail = new struct cache_entry *[DMTCHANNELS];
		_cache_size = new double [DMTCHANNELS];
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			_cache_size[tone]=0;
			_cache_head[tone]=new struct cache_entry [TONE_CACHE_SIZE];	
		}
	}
	_line_array = new struct line *[lines];
	for (int user=0;user<lines;user++) {
		_line_array[user]=get_line(user);
	} 	
	pthread_mutex_init(&_lock,NULL);	
#ifdef GSL
	A = gsl_matrix_calloc(lines,lines);
	B = gsl_vector_calloc(lines);
	X = gsl_vector_calloc(lines);
	invA = gsl_matrix_calloc(lines,lines);
	p = gsl_permutation_calloc(lines);
	tau = gsl_vector_calloc(lines);
#endif
#ifdef EIGEN
	A = MatrixXd::Identity(lines,lines);
	B = VectorXd::Zero(lines);
	x = VectorXd::Zero(lines);
#endif
}
psd_vector::~psd_vector()
{
	if (_caching_on) {
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			struct cache_entry *current=_cache_head[tone];
			struct cache_entry *next;
			while(current !=NULL) {
				next=current->next;
				delete[] current->p_vec;
				delete current;
				current=next;	
			}
		}
		delete[] _cache_head;
		delete[] _cache_tail;
		printf("This psd_vector class managed %d cache hits against %d cache misses\n",_cache_hits,_cache_misses);
		printf("calls = %d\n",_calls);
		printf("total cache size = %.0lf entries\n",_total_cache_size);
		printf("total cache in bytes = %lf kbytes\n",_total_cache_size*(double)sizeof(cache_entry)/1024);
		/*
		for (int tone=0;tone<DMTCHANNELS;tone++) {
			printf("Cache size on channel %d is %.0lf\n",tone,_cache_size[tone]);
		}
		*/
		printf("total compares = %lf\n",_compares);
		printf("total collisions = %lf\n",_collisions);
	}
#ifdef GSL
	gsl_matrix_free(A);
        gsl_matrix_free(invA);
        gsl_vector_free(X);
        gsl_vector_free(B);
        gsl_permutation_free(p);
#endif
}
int psd_vector::calc(int *b,double *gamma,int channel,double *p)
{
	//struct cache_entry* cache=_cache_head[channel];
	//bool debug=false;
	//double scale_fac = pow(2,-15);
	pthread_mutex_lock(&_lock);
	unsigned int hash=0;
	/*
	if (b[0] == 15 && b[1] == 12 && channel == 0) {
		debug=true; 
		printf("cache head on channel 0 is %p\n",cache);
		getchar();
	}
	*/
	if (_caching_on) {
		//hash = murmurhash((unsigned char *)b,lines*sizeof(int),hash);
		hash = murmurhash2((void *)b,lines*sizeof(int),hash);
	//		int index = (int)((double)hash*scale_fac);
		int index = hash % TONE_CACHE_SIZE;
		if (_cache_head[channel][index].p_vec == NULL) {
			psd_vector::calculate_psd_vector(b,gamma,channel,p);
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
		}
		else {
			if (hash != _cache_head[channel][index].hash) {	// different hashes, index collision after mod operation
				//printf("Collision!\n");
				psd_vector::calculate_psd_vector(b,gamma,channel,p);
				_cache_head[channel][index].hash=hash;
				memcpy(_cache_head[channel][index].p_vec,p,sizeof(double)*lines);
				memcpy(_cache_head[channel][index].b_vec,b,sizeof(int)*lines);
				_collisions++;
				pthread_mutex_unlock(&_lock);
				return 2;
			}
			for (int user=0;user<lines;user++) {
				if (b[user] != _cache_head[channel][index].b_vec[user]) { // proper hash collision 
					// eg (8 8 8 8 12 0 0 0) == (6 7 9 8 13 7 13 13) 
					psd_vector::calculate_psd_vector(b,gamma,channel,p);
					_cache_head[channel][index].hash=hash;
					memcpy(_cache_head[channel][index].p_vec,p,sizeof(double)*lines);
					memcpy(_cache_head[channel][index].b_vec,b,sizeof(int)*lines);
                                	_collisions++;
					pthread_mutex_unlock(&_lock);
                                	return 2;
				}
			}
			memcpy(p,_cache_head[channel][index].p_vec,sizeof(double)*lines);
			_cache_hits++;
			pthread_mutex_unlock(&_lock);
			return 1;
		}
	}
	else {
		psd_vector::calculate_psd_vector(b,gamma,channel,p);
	}
	pthread_mutex_unlock(&_lock);
	return 0;
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
	psd_vector::calculate_psd_vector(b,gamma,channel,p);
	_calls++;
	if (_caching_on) {
		insert_cache_entry(b,channel,p,hash);		// thats right! assuming gamma is the same for all entries!
		_cache_misses++;
	}	
	return 0;
	*/
}
int psd_vector::insert_cache_entry(int *b,int channel,double *p,unsigned int hash)
{
	struct cache_entry *new_ent;
	new_ent = new struct cache_entry;
	new_ent->p_vec = new double [lines];
	if (_cache_head[channel] == NULL) {
		_cache_head[channel] = new_ent;
		_cache_tail[channel] = new_ent;
	}
	else {
		_cache_tail[channel]->next = new_ent;
		_cache_tail[channel] = new_ent;
	}
	new_ent->hash=hash;
	memcpy(new_ent->p_vec,p,lines*sizeof(double));
	memcpy(new_ent->b_vec,b,lines*sizeof(int));
	new_ent->next=NULL;
	_total_cache_size++;
	_cache_size[channel]++;
	return 0;
}
#ifdef GSL
int psd_vector::calculate_psd_vector(int *b,double *gamma,int channel,double *_p)
{
	int i,j,sig;
	/*
	gsl_matrix *A = NULL;
	gsl_vector *B = NULL;
	gsl_vector *X = NULL;
	gsl_vector *tau = NULL;
	gsl_matrix *invA = NULL;
	A = gsl_matrix_calloc(lines,lines);
	B = gsl_vector_calloc(lines);
	X = gsl_vector_calloc(lines);
	tau = gsl_vector_calloc(lines);
	invA = gsl_matrix_calloc(lines,lines);
*/
	//double _A[lines][lines];
	//double _B[lines];
	//double __A[lines][lines];
	//double __B[lines];
	for (i=0;i<lines;i++) {
		for (j=0;j<lines;j++) {
			if (i==j) {
				gsl_matrix_set(A,i,j,1);
				//_A[i][j]=1;
				//__A[i][j]=1;
			}
			else {
				//printf("Matrix A element (%d,%d)\n");
				gsl_matrix_set(A,i,j,-1*F(gamma[i],b[i])*get_xtalk_gain(j,i,channel)/get_xtalk_gain(i,i,channel));		
				//_A[j][i]=-1*F(gamma[i],b[i])*get_xtalk_gain(j,i,channel)/get_xtalk_gain(i,i,channel);		
				//__A[i][j]=-1*F(gamma[i],b[i])*get_xtalk_gain(j,i,channel)/get_xtalk_gain(i,i,channel);		
			}
		}
	}
	//print_matrix(A,"A");
	for(i=0;i<lines;i++) {
		//gsl_vector_set(B,i,F(gamma[i],b[i])*(dbmhz_to_watts(-140)+alien_xtalk(i,channel))/get_xtalk_gain(i,i,channel));
		// forget about alien_xtalk for now!!
		gsl_vector_set(B,i,F(gamma[i],b[i])*(_line_array[i]->bkn+_line_array[i]->alien_xtalk[channel])/get_xtalk_gain(i,i,channel));
		//_B[i]=F(gamma[i],b[i])*(_line_array[i]->bkn+_line_array[i]->alien_xtalk[channel])/get_xtalk_gain(i,i,channel);
		//__B[i]=F(gamma[i],b[i])*(_line_array[i]->bkn+_line_array[i]->alien_xtalk[channel])/get_xtalk_gain(i,i,channel);
	}
/*
	for (i=0;i<lines;i++) {
		for (j=0;j<lines;j++) {
			printf("%g\t",_A[i][j]);	
		}
		printf("\n");
	}
	for (i=0;i<lines;i++) {
		for (j=0;j<lines;j++) {
			printf("%g\t",__A[i][j]);	
		}
		printf("\n");
	}
*/
//	getchar();
	//gsl_vector_fprintf(stdout,B,"%e");
	//
#ifdef QR
	//printf("QR method\n");
	gsl_linalg_QR_decomp(A,tau);
	gsl_linalg_QR_solve(A,tau,B,X);
#endif
	//int ipiv,info;
	//dgetrf(&lines,&lines,&_A[0][0],&lines,&ipiv,&info);
/*
	for (i=0;i<lines;i++) {
		for (j=0;j<lines;j++) {
			printf("%g\t",_A[i][j]);	
		}
		printf("\n");
	}
*/
#ifdef inv
	//printf("Direct inversion\n");
	if (gsl_linalg_LU_decomp(A,p,&sig)) {
                printf("LU decomp failed!");
                exit(1);
	}
	//print_matrix(A,"A");
	if (gsl_linalg_LU_invert(A,p,invA)) {
                printf("LU invert failed");
                exit(1);
        }
        if (gsl_blas_dgemv(CblasNoTrans,1.0,invA,B,0.0,X)) {           // ans = inv(D)*y
                printf("failed!");
                exit(1);
        }
#endif
	//char trans = 'N';
	//int n = lines;
	//int nrhs = 1;
	//int lda = lines;	
	//int ldb = lines;
	//dgetrs(&trans,&n,&nrhs,&_A[0][0],&lda,&ipiv,&_B[0],&ldb,&info);
	//dgesv(&lines,&nrhs,&__A[0][0],&lda,&ipiv,&__B[0],&ldb,&info);
	//printf("%d\n",info);
	for (i=0;i<lines;i++) {
                _p[i]=gsl_vector_get(X,i);
                //if (_p[i] > 0 && _p[i] < 4.31e-14)
                if (b[i] == 0)
                        _p[i]=0;
                //if (_p[i] < 0 && _p[i] > -4.31e-14)
                        //_p[i]=0;
        }
/*
	print_vector(b,"b");
	print_vector_psd(_p,"GSL P");
	print_vector_psd(_B,"MKL P");
	print_vector_psd(__B,"MKL P 2");
	getchar();
*/
        return 0;
}
#endif
#ifdef MKL
int psd_vector::calculate_psd_vector(int *b,double *gamma,int channel,double *_p)
{
	int i,j,sig;
	double _A[lines][lines];
	double _B[lines];
	for (i=0;i<lines;i++) {
		for (j=0;j<lines;j++) {
			if (i==j) {
				_A[i][j]=1;
			}
			else {
				//printf("Matrix A element (%d,%d)\n");
				_A[i][j]=-1*F(gamma[i],b[i])*get_xtalk_gain(j,i,channel)/get_xtalk_gain(i,i,channel);		
			}
		}
	}
	//gsl_matrix_fprintf(stdout,A,"%e");
	for(i=0;i<lines;i++) {
		//gsl_vector_set(B,i,F(gamma[i],b[i])*(dbmhz_to_watts(-140)+alien_xtalk(i,channel))/get_xtalk_gain(i,i,channel));
		// forget about alien_xtalk for now!!
		_B[i]=F(gamma[i],b[i])*(_line_array[i]->bkn+_line_array[i]->alien_xtalk[channel])/get_xtalk_gain(i,i,channel);
	}
	//gsl_vector_fprintf(stdout,B,"%e");
#ifdef QR
	//printf("QR method\n");
	gsl_linalg_QR_decomp(A,tau);
	gsl_linalg_QR_solve(A,tau,B,X);
#endif
#ifdef inv
	//printf("Direct inversion\n");
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
                //if (_p[i] > 0 && _p[i] < 4.31e-14)
                if (b[i] == 0)
                        _p[i]=0;
                //if (_p[i] < 0 && _p[i] > -4.31e-14)
                        //_p[i]=0;
        }
        return 0;
}
#endif
#ifdef EIGEN
int psd_vector::calculate_psd_vector(int *b,double *gamma,int channel,double *_p)
{
	//MatrixXd A(lines,lines);
	//VectorXd B(lines);
	//VectorXd x(lines);
	for (int i=0;i<lines;i++) {
		for (int j=0;j<lines;j++) {
			if (i==j) 
				A(i,j) =1;
			else 
				A(i,j) = -1*F(gamma[i],b[i])*get_xtalk_gain(j,i,channel)/get_xtalk_gain(i,i,channel);
		}
	}
	//std::cout << A << "\n\n";
	//
	for (int i=0;i<lines;i++) {
		B(i) = F(gamma[i],b[i])*(_line_array[i]->bkn+_line_array[i]->alien_xtalk[channel])/get_xtalk_gain(i,i,channel);
	}
	//std::cout << B << "\n\n";
	if (!(A.lu().solve(B,&x))) {
		for (int i=0;i<lines;i++) {
			_p[i] = -DBL_MAX;
		}
	}
	else {
		for (int i=0;i<lines;i++) {
			if (b[i] == 0) {
				_p[i] = 0;
			}
			else {
				_p[i] = x(i);
			}
		}
	}
	//std::cout << x << "\n\n";
	return 0;
}
#endif
unsigned int psd_vector::murmurhash(const unsigned char *data,int len,unsigned int h)
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
unsigned int psd_vector::murmurhash2( const void * key, int len, unsigned int seed )
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
#ifdef GSL
int calculate_psd_vector(double *b,double *gamma,int channel,double *_p)
{
	int i,j,sig;
	struct line* current;
	gsl_matrix *A = gsl_matrix_calloc(lines,lines);
	gsl_vector *B = gsl_vector_calloc(lines);
	gsl_vector *X = gsl_vector_calloc(lines);
	gsl_vector *tau = gsl_vector_calloc(lines);
	gsl_matrix *invA = gsl_matrix_calloc(lines,lines);
	gsl_permutation *p = gsl_permutation_calloc(lines);
	for (i=0;i<lines;i++) {
		for (j=0;j<lines;j++) {
			if (i==j)
				gsl_matrix_set(A,i,j,1);
			else {
				//printf("Matrix A element (%d,%d)\n");
				gsl_matrix_set(A,i,j,-1*F(gamma[i],b[i])*get_xtalk_gain(j,i,channel)/get_xtalk_gain(i,i,channel));		
			}
		}
	}
	//gsl_matrix_fprintf(stdout,A,"%e");
	for(i=0;i<lines;i++) {
		//gsl_vector_set(B,i,F(gamma[i],b[i])*(dbmhz_to_watts(-140)+alien_xtalk(i,channel))/get_xtalk_gain(i,i,channel));
		// forget about alien xtalk for now!!!!
		gsl_vector_set(B,i,F(gamma[i],b[i])*(dbmhz_to_watts(-140)+alien_xtalk(i,channel))/get_xtalk_gain(i,i,channel));
	}
	//gsl_vector_fprintf(stdout,B,"%e");
#ifdef QR
	//printf("QR method\n");
	gsl_linalg_QR_decomp(A,tau);
	gsl_linalg_QR_solve(A,tau,B,X);
#endif
#ifdef inv
	//printf("Direct inversion\n");
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
	//gsl_vector_fprintf(stdout,X,"%e");	
/*
	for (i=0;i<lines;i++) {
		current=get_line(i);
		current->b[channel] = b[i];
		current->psd[channel] = watts_to_dbmhz(gsl_vector_get(X,i));
		current->gamma[0] = gamma[i];
	}
*/
	//calculate_snr();
/*
	for (i=0;i<lines;i++) {	
		current=get_line(i);
		current->print_margin(channel);	
		current->print_psd(channel);	
	}
*/
	//bool ooops=false;
	for (i=0;i<lines;i++) {
		_p[i]=gsl_vector_get(X,i);
		if (_p[i] > 0 && _p[i] < 4.31e-14)
			_p[i]=0;
		if (_p[i] < 0 && _p[i] > -4.31e-14)
			_p[i]=0;
		/*if (_p[i] < 0) {
			ooops=true;				
		}*/
	}
	/*
	if (ooops) {
		printf("Ooops in calculate_psd_vector\n");
		print_vector(b,"b");	
		print_vector(gamma,"gamma");	
		print_vector(_p,"p");
		getchar();
	}
	*/
	gsl_matrix_free(A);
	gsl_matrix_free(invA);
	gsl_vector_free(X);
	gsl_vector_free(B);
	gsl_vector_free(tau);
	gsl_permutation_free(p);	
	return 0;
}
#endif
#ifdef EIGEN
int calculate_psd_vector(double *b,double *gamma,int channel,double *_p)
{
	MatrixXd A(lines,lines);
	VectorXd B(lines);
	VectorXd x(lines);
	for (int i=0;i<lines;i++) {
		for (int j=0;j<lines;j++) {
			if (i==j) 
				A(i,j) =1;
			else 
				A(i,j) = -1*F(gamma[i],b[i])*get_xtalk_gain(j,i,channel)/get_xtalk_gain(i,i,channel);
		}
	}
	//std::cout << A << "\n\n";
	//
	for (int i=0;i<lines;i++) {
		B(i) = F(gamma[i],b[i])*(line_array[i]->bkn+line_array[i]->alien_xtalk[channel])/get_xtalk_gain(i,i,channel);
	}
	//std::cout << B << "\n\n";
	if (!(A.lu().solve(B,&x))) {
		for (int i=0;i<lines;i++) {
			_p[i] = -DBL_MAX;
		}
	}
	else {
		for (int i=0;i<lines;i++) {
			if (b[i] == 0) {
				_p[i] = 0;
			}
			else {
				_p[i] = x(i);
			}
		}
	}
	//std::cout << x << "\n\n";
	return 0;
}
#endif
int calculate_psd_vector(int *b,double *gamma,int channel,double *_p)
{
	int i,j,sig;
	struct line* current;
	gsl_matrix *A = gsl_matrix_calloc(lines,lines);
	gsl_vector *B = gsl_vector_calloc(lines);
	gsl_vector *X = gsl_vector_calloc(lines);
	//gsl_vector *tau = gsl_vector_calloc(lines);
	gsl_matrix *invA = gsl_matrix_calloc(lines,lines);
	gsl_permutation *p = gsl_permutation_calloc(lines);
	for (i=0;i<lines;i++) {
		for (j=0;j<lines;j++) {
			if (i==j)
				gsl_matrix_set(A,i,j,1);
			else {
				//printf("Matrix A element (%d,%d)\n");
				gsl_matrix_set(A,i,j,-1*F(gamma[i],b[i])*get_xtalk_gain(j,i,channel)/get_xtalk_gain(i,i,channel));		
			}
		}
	}
	//gsl_matrix_fprintf(stdout,A,"%e");
	for(i=0;i<lines;i++) {
		//gsl_vector_set(B,i,F(gamma[i],b[i])*(dbmhz_to_watts(-140)+alien_xtalk(i,channel))/get_xtalk_gain(i,i,channel));
		// forget about alien_xtalk for now!!
		gsl_vector_set(B,i,F(gamma[i],b[i])*(dbmhz_to_watts(-140)+alien_xtalk(i,channel))/get_xtalk_gain(i,i,channel));
	}
	//gsl_vector_fprintf(stdout,B,"%e");
#ifdef QR
	//printf("QR method\n");
	gsl_linalg_QR_decomp(A,tau);
	gsl_linalg_QR_solve(A,tau,B,X);
#endif
#ifdef inv
	//printf("Direct inversion\n");
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
	//gsl_vector_fprintf(stdout,X,"%e");	
/*
	for (i=0;i<lines;i++) {
		current=get_line(i);
		current->b[channel] = b[i];
		current->psd[channel] = watts_to_dbmhz(gsl_vector_get(X,i));
		current->gamma[0] = gamma[i];
	}
*/
	//calculate_snr();
/*
	for (i=0;i<lines;i++) {	
		current=get_line(i);
		current->print_margin(channel);	
		current->print_psd(channel);	
	}
*/
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
	//gsl_vector_free(tau);
	gsl_permutation_free(p);	
	return 0;
}
