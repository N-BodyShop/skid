#ifndef SMOOTH_HINCLUDED
#define SMOOTH_HINCLUDED

#include "kd.h"

#define RESMOOTH_SAFE  10


typedef struct pqNode {
	float fKey;
	struct pqNode *pqLoser;
	struct pqNode *pqFromInt;
	struct pqNode *pqFromExt;
	struct pqNode *pqWinner;	/* Only used when building initial tree */
	int p;
	float ax;
	float ay;
	float az;
	} PQ;


typedef struct nNeighborList {
	int p;
	float fDist2;
	float dx;
	float dy;
	float dz;
	} NN;


typedef struct smContext {
	KD kd;
	int nSmooth;
	PQ *pq;
	PQ *pqHead;
	char *iMark;
	} * SMX;


#define PQ_INIT(pq,n)\
{\
	int j;\
	for (j=0;j<(n);++j) {\
		if (j < 2) (pq)[j].pqFromInt = NULL;\
		else (pq)[j].pqFromInt = &(pq)[j>>1];\
		(pq)[j].pqFromExt = &(pq)[(j+(n))>>1];\
		(pq)[j].p = 0;\
		}\
	}


#define PQ_BUILD(pq,n,q)\
{\
	int i,j;\
	PQ *t,*lt;\
	for (j=(n)-1;j>0;--j) {\
		i = (j<<1);\
		if (i < (n)) t = (pq)[i].pqWinner;\
		else t = &(pq)[i-(n)];\
		++i;\
		if (i < (n)) lt = (pq)[i].pqWinner;\
		else lt = &(pq)[i-(n)];\
		if (t->fKey < lt->fKey) {\
			(pq)[j].pqLoser = t;\
			(pq)[j].pqWinner = lt;\
			}\
		else {\
			(pq)[j].pqLoser = lt;\
			(pq)[j].pqWinner = t;\
			}\
		}\
	(q) = (pq)[1].pqWinner;\
	}


#define PQ_REPLACE(q)\
{\
    PQ *t,*lt;\
	t = (q)->pqFromExt;\
	while (t) {\
		if (t->pqLoser->fKey > (q)->fKey) {\
			lt = t->pqLoser;\
			t->pqLoser = (q);\
			(q) = lt;\
			}\
		t = t->pqFromInt;\
		}\
	}



int smInit(SMX *,KD,int);
void smFinish(SMX);
void smBallSearch(SMX,float,float *);
int  smBallGather(SMX,float,float *,NN *);
void smDensityInit(SMX);
void smAccDensity(SMX);

#endif



