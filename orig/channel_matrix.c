#include "multiuser_load.h"
#include <cstdlib>
#include <cmath>

double offset[4][4] = {{0,-14.47,-4.47,-1.76},
		{-15.54,0,-13.27,-4.57},
		{-4.81,-13.85,0,-3.37},
		{-2.79,-4.04,-4.61,0}};

double offset2[8][8] = 	{{0,-17.3757,-15.1871,-16.237,-13.0119,2.1896,-12.9727,-15.5395},
			{-17.2496,0,-7.3152,-24.1604,-15.3275,-8.1194,-21.5641,-13.0353},
			{-14.7051,-8.742,0,-12.3908,-11.9613,-6.3329,-5.6299,-18.2671},
			{-15.4583,-25.0632,-13.7968,0,-15.7089,-37.6546,-2.859,-10.2667},
			{-13.6823,-15.1213,-10.6057,-16.8391,0,-14.5002,-7.9387,-17.6145},
			{3.2564,-9.1261,-6.0911,-38.9937,-13.5543,0,-27.2996,-20.5998},
			{-13.1577,-22.6951,-5.635,-3.7635,-7.1681,-26.1968,0,-18.6713},
			{-14.4681,-13.3279,-19.7019,-8.9933,-17.2753,-19.2207,-18.4048,0}};

double *freq_matrix;

bool beta_model;
int DMTCHANNELS;


double *calc_channel_matrix(const int lines)
{
	int i,j,k;
	double freq[DMTCHANNELS];
	double *matrix;
	struct line* current;
	
	freq_matrix = new double[DMTCHANNELS];

	if ((matrix = (double *)malloc(sizeof(double) * lines * lines * DMTCHANNELS)) == NULL) {
                printf("Cannot alloc memory for channel matrix\n");
                exit(1);
        }

#ifdef ADSL_DOWNSTREAM	
	for (i=0;i<DMTCHANNELS;i++) {
		//freq[i] = 2156.25 + 4312.5 * i;
		//printf("Starting from 2.156kHz\n");
		freq[i] = 140156.25 + 4312.5 * i;	// quick hack to only use downstream ADSL channels, FIXME!! use with DMTCHANNELS=224
		freq_matrix[i] = 140156.25 + 4312.5 * i;	// quick hack to only use downstream ADSL channels, FIXME!! use with DMTCHANNELS=224
		//printf("Starting from 140156.25kHz\n");
	}
#endif

#ifdef VDSL_UPSTREAM	
	for (i=0;i<LOWERCHANNELS;i++) {
		//freq[i] = 2156.25 + 4312.5 * i;
		//printf("Starting from 2.156kHz\n");
		freq[i] = 3.750e6 + 4312.5 * i;	// quick hack to only use downstream ADSL channels, FIXME!! use with DMTCHANNELS=224
		freq_matrix[i] = 3.750e6 + 4312.5 * i;	// quick hack to only use downstream ADSL channels, FIXME!! use with DMTCHANNELS=224
		//printf("Starting from 3.75MHz\n");
	}
	for (i=0;i<UPPERCHANNELS;i++) {
		//freq[i] = 2156.25 + 4312.5 * i;
		//printf("Starting from 2.156kHz\n");
		freq[i+LOWERCHANNELS] = 8.5e6 + 4312.5 * i;	// quick hack to only use downstream ADSL channels, FIXME!! use with DMTCHANNELS=224
		freq_matrix[i+LOWERCHANNELS] = 8.5e6 + 4312.5 * i;	// quick hack to only use downstream ADSL channels, FIXME!! use with DMTCHANNELS=224
		//printf("Starting from 8.5MHz\n");
	}
#endif

	
        for (k=0;k<DMTCHANNELS;k++) {
                for (j=0;j<lines;j++) {
                        for (i=0;i<lines;i++) {
                                if (i == j) { 
                                        //channel_matrix[i][j][k] = transfer_fn(i,freq[k]);       // i == line_id
					current = get_line(i);
                                        current->gain[k] = xtalk_gain(i,j,k) = transfer_fn(i,freq[k]);       // i == line_id
				}
                                else    { // victim line_id is i, crosstalker line_id is j 
                                        //channel_matrix[i][j][k]  = calc_fext_xtalk_gain(i,j,freq[k],DOWNSTREAM);
                                        xtalk_gain(j,i,k)  = calc_fext_xtalk_gain(i,j,freq[k],DOWNSTREAM);
					//if (lines>2) {
					if (beta_model) {
						if (lines>2 && lines <= 4) {
							//channel_matrix[i][j][k] *= UNdB(offset[i][j]);       // i == line_id
							xtalk_gain(j,i,k) *= UNdB(offset[i][j]);       // i == line_id
						}
						else if (lines > 4) {
							//channel_matrix[i][j][k] *= UNdB(offset2[i][j]);       // i == line_id
							xtalk_gain(j,i,k) *= UNdB(offset2[i][j]);       // i == line_id
						}
					}
				}
                        }
                }
        }

	printf("Finished the channel matrix!\n");

	return matrix;
	
}

void check_normalised_xtalk_gains()
{

	int ycount=0,ncount=0;

	for (int tone=0;tone<DMTCHANNELS;tone++) {
		for (int victim=0;victim<lines;victim++) {
			for (int xtalker=0;xtalker<lines;xtalker++) {
				if (xtalker == victim)
					continue;
				if (get_xtalk_gain(xtalker,victim,tone)/get_channel_gain(victim,tone) > 0.5)
					ycount++;
				else
					ncount++;
			}
		}
	}

	printf("Normalised xtalk gain was greater than 1/2 on %d channels versus %d channels\n",ycount,ncount);

}

double calc_fext_xtalk_gain(int v, int x, double freq, int dir)
{

	struct line* xtalker = list_head;
	struct line* victim = list_head;
	int _case;			// case is a reserved word
	double h1,h2,h3;
	double ret=-1;
	double f;


	xtalker = get_line(x);
	victim = get_line(v);

	_case = get_case(xtalker,victim,dir);


	switch(_case) {
		case 1:
			//h1 = insertion_loss(2,(victim->lt - xtalker->lt)/1000,freq);
			h1 = do_transfer_fn(2,(victim->lt - xtalker->lt)/1000,freq);
			h2 = insertion_loss(2,(xtalker->nt - victim->lt)/1000,freq);
			//h3 = insertion_loss(2,(victim->nt - xtalker->nt)/1000,freq);
			h3 = do_transfer_fn(2,(victim->nt - xtalker->nt)/1000,freq);
			//ret = h1 * fext(freq,(xtalker->nt - victim->lt)/1000,h2) * h3;
			ret = h1 * h2 * h3 * fext(freq,(xtalker->nt - victim->lt)/1000);
			break;
		case 2:	
			h2 = insertion_loss(2,xtalker->length/1000,freq);
			//h3 = insertion_loss(2,(victim->length - xtalker->length)/1000,freq);
			h3 = do_transfer_fn(2,(victim->length - xtalker->length)/1000,freq);
			//ret = h3 * fext(freq,xtalker->length/1000,h2);
			ret = h2 * h3 * fext(freq,xtalker->length/1000);
			break;
		case 3:
			h2 = insertion_loss(2,xtalker->length/1000,freq);
			//h3 = insertion_loss(2,(xtalker->lt - victim->lt)/1000,freq);
			h3 = do_transfer_fn(2,(victim->length - xtalker->length)/1000,freq);
			//ret = h3 * fext(freq,xtalker->length/1000,h2);
			ret = h2 * h3 * fext(freq,xtalker->length/1000);
			break;
		case 4:
			//h1 = insertion_loss(2,(xtalker->length - victim->length)/1000,freq);
			h1 = do_transfer_fn(2,(xtalker->length - victim->length)/1000,freq);
			h2 = insertion_loss(2,victim->length/1000,freq);
			//ret = h1 * fext(freq,victim->length/1000,h2);
			ret = h1 * h2 * fext(freq,victim->length/1000);
			break;
		case 5:
			h2 = insertion_loss(2,xtalker->length/1000,freq);
			//ret = fext(freq,xtalker->length/1000,h2);
			ret = h2 * fext(freq,xtalker->length/1000);
			break;
		case 6:
			h2 = insertion_loss(2,xtalker->length/1000,freq);
			//ret = fext(freq,xtalker->length/1000,h2);
			ret = h2 * fext(freq,xtalker->length/1000);
			break;
		case 7:
			//h1 = insertion_loss(2,(xtalker->nt - victim->nt)/1000,freq);
			h1 = do_transfer_fn(2,(xtalker->nt - victim->nt)/1000,freq);
			h2 = insertion_loss(2,victim->length/1000,freq);
			//ret = h1 * fext(freq,victim->length/1000,h2);
			ret = h1 * h2 * fext(freq,victim->length/1000);
			break;
		case 8:	
			h2 = insertion_loss(2,victim->length/1000,freq);
			//ret = fext(freq,victim->length/1000,h2);
			ret = h2 * fext(freq,victim->length/1000);
			break;
		case 9:
			h2 = insertion_loss(2,(victim->nt - xtalker->lt)/1000,freq);
			//ret = fext(freq,(victim->nt - xtalker->lt)/1000,h2);
			f=fext(freq,(victim->nt - xtalker->lt)/1000);
			ret = h2 * f;
			//printf("%g\th2 = %g\tfext = %g\t h2 * fext = %g\n",freq,10*log10(h2),10*log10(f),10*log10(h2*f));
			break;
		default:
			printf("uh oh\n");
			getchar();
			break;
	} 

	

	return ret;
	

}



double fext(double f,double l)
{

	double k_xf=0.0056;
	double l_0=1;
	double f_0=1e6;

	//double x=k_xf * f/f_0 * sqrt(l/l_0) * h;
	
	double x=-55;		// model from BT's simulation parameters document cp38-2
	double n=((double)(lines-1))/49;	// n disturbers
	x+= 20*log10(f/90e3);		// f in Hz
	x+= 6*log10(n);			// shared l in km
	x+= 10*log10(l);

	if (l < 0)	
		x=0;
	//return 0;
	//return UNdB(x)*UNdB(-10)*(1+f/2000);
	//return UNdB(x)*UNdB(-10)*(1+pow(f,2)/50e7);
	return UNdB(x);

}

int get_case(struct line* xtalker, struct line* victim, int dir)
{


	int _case;


	//if ((victim->lt == 0) && (xtalker->lt == 0) && (victim->length < xtalker->length)) {  	// case 4 upstream, case 8 downstream
	if ((victim->lt == xtalker->lt) && (victim->length < xtalker->length)) {  	// case 4 upstream, case 8 downstream
		if (dir == UPSTREAM)
			_case = 4;
		if (dir == DOWNSTREAM)
			_case = 8;
	}
	else if ((victim->nt == xtalker->nt) && (xtalker->length > victim->length)) {		// case 8 up, case 4 down
		if (dir == UPSTREAM)
			_case = 8;
		if (dir == DOWNSTREAM)
			_case = 4;
	}
	else if ((victim->nt == xtalker->nt)	&& (victim->length > xtalker->length)) {	// case 2 up, case 6 down
		if (dir == UPSTREAM)
			_case = 2;
		if (dir == DOWNSTREAM)
			_case = 6;
		
	}
	//else if ((victim->lt == 0) && (xtalker->lt == 0) && (victim->length > xtalker->length)) {	// case 6 up, case 2 down
	else if ((victim->lt == xtalker->lt) && (victim->length > xtalker->length)) {	// case 6 up, case 2 down
		if (dir == UPSTREAM)
			_case = 6;
		if (dir == DOWNSTREAM)
			_case = 2;
	}
	else if ((victim->lt > xtalker->lt) && (victim->nt < xtalker->nt)) {		// case 7
		_case=7;		
	}
	else if ((victim->lt < xtalker->lt) && (victim->nt > xtalker->nt)) {		// case 3
		_case=3;
	}
	else if ((victim->lt == xtalker->lt) && (victim->nt == xtalker->nt)) {		// case 5!
		_case=5;		
	}
	else if ((victim->lt > xtalker->lt) && (victim->nt > xtalker->nt)) {		// case 9 up, case 1 down
		if (dir == UPSTREAM)
			_case=9;
		if (dir == DOWNSTREAM)
			_case=1;
	}
	else if ((victim->lt < xtalker->lt) && (victim->nt < xtalker->nt)) {		// case 1 up, case 9 down
		if (dir == UPSTREAM)
			_case=1;
		if (dir == DOWNSTREAM)
			_case=9;
	}
	
	else if ((victim->nt <= xtalker->lt)) {		// No Fext 1  lines never overlap
		_case=0;
	}
	else if ((victim->lt >= xtalker->nt)) {		// No Fext 2  lines never overlap	
		_case=0;
	}

	else {
		printf("ooops havent implemented this case yet");
		printf("xtalker = %d victim = %d\n",xtalker->line_id,victim->line_id);
		printf("victim->nt = %lf victim->lt = %lf xtalker->nt = %lf xtalker->lt = %lf\n",victim->nt,victim->lt,xtalker->nt,xtalker->lt);
		getchar();
		_case=-1;
	}


	return _case;

}

		
/*double index(double *matrix,int lines,int i,int j,int k) 
{	
	return *(matrix + k + (DMTCHANNELS*j) + (DMTCHANNELS * lines));
}*/

/*
#ifdef OLD

double calc_fext_xtalk_gain(int v, int x, double freq)
{

	struct line* xtalker = list_head;
	struct line* victim = list_head;
	double h1,h2,h3;
	double ret=-1;

	while (xtalker->line_id != x) {
		xtalker=xtalker->next;
	}

	while (victim->line_id != v) {
		victim=victim->next;
	}

	if ((victim->lt == 0) && (xtalker->lt == 0) && (victim->length < xtalker->length)) {  	// case 4 upstream, case 8 downstream
		h1 = insertion_loss(2,(xtalker->length - victim->length)/1000,freq);
		h2 = insertion_loss(2,victim->length/1000,freq);
		ret = h1 * fext(freq,victim->length/1000,h2);
	}

	else if ((victim->nt == xtalker->nt) && (xtalker->length > victim->length)) {		// case 8 up, case 4 down
		h2 = insertion_loss(2,victim->length/1000,freq);
		ret = fext(freq,victim->length/1000,h2);
	}

	else if ((victim->nt == xtalker->nt)	&& (victim->length > xtalker->length)) {	// case 2 up, case 6 down
		h2 = insertion_loss(2,xtalker->length/1000,freq);
		h3 = insertion_loss(2,(victim->length - xtalker->length)/1000,freq);
		ret = h3 * fext(freq,xtalker->length/1000,h2);
	}

	else if ((victim->lt == 0) && (xtalker->lt == 0) && (victim->length > xtalker->length)) {	// case 6 up, case 2 down
		h2 = insertion_loss(2,xtalker->length/1000,freq);
		ret = fext(freq,xtalker->length/1000,h2);
	}

	else if ((victim->lt > xtalker->lt) && (victim->nt < xtalker->nt)) {		// case 7
		h1 = insertion_loss(2,(xtalker->nt - victim->nt)/1000,freq);
		h2 = insertion_loss(2,victim->length/1000,freq);
		ret = h1 * fext(freq,victim->length/1000,h2);
	}

	else if ((victim->lt < xtalker->lt) && (victim->nt > xtalker->nt)) {		// case 3
		h2 = insertion_loss(2,xtalker->length/1000,freq);
		h3 = insertion_loss(2,(xtalker->lt - victim->lt)/1000,freq);
		ret = h3 * fext(freq,xtalker->length/1000,h2);
	}

	else if ((victim->lt == xtalker->lt) && (victim->nt == xtalker->nt)) {		// case 5!
		h2 = insertion_loss(2,xtalker->length/1000,freq);
		ret = fext(freq,xtalker->length/1000,h2);
	}

	else if ((victim->nt <= xtalker->lt)) {		// No Fext 1  lines never overlap
		ret = 0;
	}

	else if ((victim->lt >= xtalker->nt)) {		// No Fext 2  lines never overlap	
		ret = 0;
	}

	else {
		printf("ooops havent implemented this case yet");
		getchar();
	}	

	return ret;
	

}

#endif
*/

