#include "isb_new.h"

class osb: public isb
{
	public:
	osb();
	int run();

	private:
	void optimise_p();
	void init_lines();
	double osb_lk(int *,int);
	double calc_w_distance(void);
	void update_w(void);	

	double s_w;
	double mask;
	

};


