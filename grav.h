#ifndef GRAV_HINCLUDED
#define GRAV_HINCLUDED


/*
 ** see (A1) and (A2) of TREESPH: A UNIFICATION OF SPH WITH THE 
 ** HIERARCHICAL TREE METHOD by Lars Hernquist and Neal Katz.
 ** APJ Supplemant Series 70:416-446, 1989
 ** 
 */
#define SPLINE_POT(r2,twoh,a)\
{\
	double SPLINE_r,SPLINE_u,SPLINE_dih,SPLINE_dir;\
	SPLINE_r = sqrt(r2);\
	if (SPLINE_r < twoh) {\
		SPLINE_dih = 2.0/twoh;\
		SPLINE_u = SPLINE_r*SPLINE_dih;\
		if (SPLINE_u < 1.0) {\
			a = SPLINE_dih*(7.0/5.0 - 2.0/3.0*SPLINE_u*SPLINE_u + 3.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u\
					 - 1.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u);\
			}\
		else {\
			SPLINE_dir = 1.0/SPLINE_r;\
			a = -1.0/15.0*SPLINE_dir + SPLINE_dih*(8.0/5.0 - 4.0/3.0*SPLINE_u*SPLINE_u + SPLINE_u*SPLINE_u*SPLINE_u\
			              - 3.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u + 1.0/30.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u);\
			}\
		}\
	else {\
		a = 1.0/SPLINE_r;\
		}\
	}


void kdCellPot(KD,PINIT *,int,int,double *);
void kdSubPot(KD,PINIT *,int,PINIT *,int,double *);
void kdAddScoopPot(KD,PINIT *,int,float *,float,float *,int,double *);

#endif
