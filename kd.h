#ifndef KD_HINCLUDED
#define KD_HINCLUDED

#define ROOT		1
#define LOWER(i)	(i<<1)
#define UPPER(i)	((i<<1)+1)
#define PARENT(i)	(i>>1)
#define SIBLING(i) 	((i&1)?i-1:i+1)
#define SETNEXT(i)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}

#define DARK	1
#define GAS		2
#define STAR	4

/*
 ** Softening types!
 */
#define PLUMMER		1
#define SPLINE		2

typedef struct pInitial {
	float r[3];
	float v[3];
	float fMass;
	float fSoft;
	float fTemp;
	float fBall2;
	float fDensity;
	int iOrder;
	} PINIT;

typedef struct pMoved {
	float r[3];
	float rOld[3];
	float a[3];
	int iOrder;
	} PMOVE;

typedef struct pGroup {
	float rcm[3];
	float vcm[3];
	float rCenter[3];
	float rel[3];
	float fMass;
	float fRadius;
	} PGROUP;

typedef struct bndBound {
	float fMin[3];
	float fMax[3];
	} BND;

typedef struct kdNode {
	float fSplit;
	BND bnd;
	int iDim;
	int pLower;
	int pUpper;
	} KDN;


typedef struct kdContext {
	int nBucket;
	float fPeriod[3];
	float fCenter[3];
	float G;
	float Omega0;
	float H0;
	float z;
	int nParticles;
	int nDark;
	int nGas;
	int nStar;
	int inType;
	float fTime;
	int nLevels;
	int nNodes;
	int nSplit;
	int nMove;
	int nActive;
	int nInitActive;
	PINIT *pInit;
	PMOVE *pMove;
	PGROUP *pGroup;
	KDN *kdNodes;
	int nGroup;
	int *piGroup;
	KDN *kdGroup;
	int uSecond;
	int uMicro;
	} * KD;


#define INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz)\
{\
	float dx,dy,dz,dx1,dy1,dz1,fDist2;\
	dx = c[cp].bnd.fMin[0]-x;\
	dx1 = x-c[cp].bnd.fMax[0];\
	dy = c[cp].bnd.fMin[1]-y;\
	dy1 = y-c[cp].bnd.fMax[1];\
	dz = c[cp].bnd.fMin[2]-z;\
	dz1 = z-c[cp].bnd.fMax[2];\
	if (dx > 0.0) {\
		dx1 += lx;\
		if (dx1 < dx) {\
			fDist2 = dx1*dx1;\
			sx = x+lx;\
			}\
		else {\
			fDist2 = dx*dx;\
			sx = x;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dx1 > 0.0) {\
		dx += lx;\
		if (dx < dx1) {\
			fDist2 = dx*dx;\
			sx = x-lx;\
			}\
		else {\
			fDist2 = dx1*dx1;\
			sx = x;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		fDist2 = 0.0;\
		sx = x;\
		}\
	if (dy > 0.0) {\
		dy1 += ly;\
		if (dy1 < dy) {\
			fDist2 += dy1*dy1;\
			sy = y+ly;\
			}\
		else {\
			fDist2 += dy*dy;\
			sy = y;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dy1 > 0.0) {\
		dy += ly;\
		if (dy < dy1) {\
			fDist2 += dy*dy;\
			sy = y-ly;\
			}\
		else {\
			fDist2 += dy1*dy1;\
			sy = y;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sy = y;\
		}\
	if (dz > 0.0) {\
		dz1 += lz;\
		if (dz1 < dz) {\
			fDist2 += dz1*dz1;\
			sz = z+lz;\
			}\
		else {\
			fDist2 += dz*dz;\
			sz = z;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dz1 > 0.0) {\
		dz += lz;\
		if (dz < dz1) {\
			fDist2 += dz*dz;\
			sz = z-lz;\
			}\
		else {\
			fDist2 += dz1*dz1;\
			sz = z;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sz = z;\
		}\
	}


#define INTERCONT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz)\
{\
	float dx,dy,dz,dx1,dy1,dz1,fDist2,fMax2;\
	dx = c[cp].bnd.fMin[0]-x;\
	dx1 = x-c[cp].bnd.fMax[0];\
	dy = c[cp].bnd.fMin[1]-y;\
	dy1 = y-c[cp].bnd.fMax[1];\
	dz = c[cp].bnd.fMin[2]-z;\
	dz1 = z-c[cp].bnd.fMax[2];\
	if (dx > 0.0) {\
		if (dx1+lx < dx) {\
			dx1 += lx;\
			dx -= lx;\
			sx = x+lx;\
			fDist2 = dx1*dx1;\
			fMax2 = dx*dx;\
			}\
		else {\
			sx = x;\
			fDist2 = dx*dx;\
			fMax2 = dx1*dx1;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dx1 > 0.0) {\
		if (dx+lx < dx1) {\
		    dx += lx;\
			dx1 -= lx;\
			sx = x-lx;\
			fDist2 = dx*dx;\
			fMax2 = dx1*dx1;\
			}\
		else {\
			sx = x;\
			fDist2 = dx1*dx1;\
			fMax2 = dx*dx;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sx = x;\
		fDist2 = 0.0;\
		if (dx < dx1) fMax2 = dx*dx;\
		else fMax2 = dx1*dx1;\
		}\
	if (dy > 0.0) {\
		if (dy1+ly < dy) {\
		    dy1 += ly;\
			dy -= ly;\
			sy = y+ly;\
			fDist2 += dy1*dy1;\
			fMax2 += dy*dy;\
			}\
		else {\
			sy = y;\
			fDist2 += dy*dy;\
			fMax2 += dy1*dy1;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dy1 > 0.0) {\
		if (dy+ly < dy1) {\
		    dy += ly;\
			dy1 -= ly;\
			sy = y-ly;\
			fDist2 += dy*dy;\
			fMax2 += dy1*dy1;\
			}\
		else {\
			sy = y;\
			fDist2 += dy1*dy1;\
			fMax2 += dy*dy;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sy = y;\
		if (dy < dy1) fMax2 += dy*dy;\
		else fMax2 += dy1*dy1;\
		}\
	if (dz > 0.0) {\
		if (dz1+lz < dz) {\
		    dz1 += lz;\
            dz -= lz;\
			sz = z+lz;\
			fDist2 += dz1*dz1;\
			fMax2 += dz*dz;\
			}\
		else {\
			sz = z;\
			fDist2 += dz*dz;\
			fMax2 += dz1*dz1;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dz1 > 0.0) {\
		if (dz+lz < dz1) {\
			dz += lz;\
		    dz1 -= lz;\
			sz = z-lz;\
			fDist2 += dz*dz;\
			fMax2 += dz1*dz1;\
			}\
		else {\
			sz = z;\
			fDist2 += dz1*dz1;\
			fMax2 += dz*dz;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sz = z;\
		if (dz < dz1) fMax2 += dz*dz;\
		else fMax2 += dz1*dz1;\
		}\
	if (fMax2 < fBall2) goto ContainedCell;\
	}


void kdTime(KD,int *,int *);
int kdInit(KD *,int,float *,float *);
void kdSetSoft(KD,float);
void kdSetUniverse(KD,float,float,float,float);
int kdParticleType(KD,int);
int kdReadTipsy(KD,FILE *);
int kdBuildTree(KD);
int kdBuildMoveTree(KD);
void kdInitMove(KD,float,float,float);
int kdScatterActive(KD);
int kdScatterCut(KD);
void kdMoveParticles(KD,float);
int kdPruneInactive(KD,float);
void kdReactivateMove(KD);
void kdFoF(KD,float);
void kdInGroup(KD,char *);
void kdGroupOrder(KD);
void kdTooSmall(KD,int);
void kdUnbind(KD,int,float);
void kdOutGroup(KD,char *);
void kdOutDensity(KD,char *);
void kdOutVector(KD,char *);
void kdWriteGroup(KD,char *);
void kdFinish(KD);

#endif





