#define P_MAT p_mat[2][MAXBITSPERTONE+1][MAXBITSPERTONE+1][DMTCHANNELS]
#define B_MAT b_mat[2][DMTCHANNELS]

struct osb_2_params *osb_2(struct osb_2_params *p);

struct osb_2_params {

	osb_2_params();
	void osb_2_params_init();

        double p_mat[2][MAXBITSPERTONE+1][MAXBITSPERTONE+1][DMTCHANNELS];
        int b_mat[2][DMTCHANNELS];
        double w;
        double w_max;
        double w_min;
        double l1;
        double l1_min;
        double l1_max;
        double l2;
        double l2_min;
        double l2_max;
        int rate0_target;
        int rate1_target;
	double p0_budget;
        double p1_budget;

	int e;
	
	void osb_init_lines();	

};
