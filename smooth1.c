#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include "smooth1.h"


int smInit(SMX *psmx,KD kd,int nSmooth)
{
	SMX smx;

	assert(nSmooth <= kd->nInitActive);
	smx = (SMX)malloc(sizeof(struct smContext));
	assert(smx != NULL);
	smx->kd = kd;
	smx->nSmooth = nSmooth;
	smx->pq = (PQ *)malloc(nSmooth*sizeof(PQ));
	assert(smx->pq != NULL);
	PQ_INIT(smx->pq,nSmooth);
	smx->iMark = (char *)malloc(kd->nInitActive*sizeof(char));
	assert(smx->iMark);
	smx->pp = NULL;
	smx->nExtraScat = 0;
	/*
	 ** Initialize arrays for calculated quantities.
	 */
	*psmx = smx;	
	return(1);
	}


void smFinish(SMX smx)
{
	if (smx->pp) free(smx->pp);
	free(smx->iMark);
	free(smx->pq);
	free(smx);
	}


void smBallSearch(SMX smx,float fBall2,float *ri)
{
	KDN *c;
	PINIT *p;
	int cell,cp,ct,pj;
	float fDist2,dx,dy,dz,lx,ly,lz,sx,sy,sz,x,y,z;
	PQ *pq;

	c = smx->kd->kdNodes;
	if (!c) return;
	p = smx->kd->pInit;
	pq = smx->pqHead;
	x = ri[0];
	y = ri[1];
	z = ri[2];
	lx = smx->kd->fPeriod[0];
	ly = smx->kd->fPeriod[1];
	lz = smx->kd->fPeriod[2];
	cell = ROOT;
	/*
	 ** First find the "local" Bucket.
	 ** This could mearly be the closest bucket to ri[3].
	 */
	while (c[cell].iDim >= 0) {
		if (ri[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
		else cell = UPPER(cell);
		}
	/*
	 ** Now start the search from the bucket given by cell!
	 */
	for (pj=c[cell].pLower;pj<=c[cell].pUpper;++pj) {
		dx = x - p[pj].r[0];
		dy = y - p[pj].r[1];
		dz = z - p[pj].r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < fBall2) {
			if (smx->iMark[pj]) continue;
			smx->iMark[pq->p] = 0;
			smx->iMark[pj] = 1;
			pq->fKey = fDist2;
			pq->p = pj;
			pq->ax = 0.0;
			pq->ay = 0.0;
			pq->az = 0.0;
			PQ_REPLACE(pq);
			fBall2 = pq->fKey;
			}
		}
	while (cell != ROOT) {
		cp = SIBLING(cell);
		ct = cp;
		SETNEXT(ct);
		while (1) {
			INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz);
			/*
			 ** We have an intersection to test.
			 */
			if (c[cp].iDim >= 0) {
				cp = LOWER(cp);
				continue;
				}
			else {
				for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
					dx = sx - p[pj].r[0];
					dy = sy - p[pj].r[1];
					dz = sz - p[pj].r[2];
					fDist2 = dx*dx + dy*dy + dz*dz;
					if (fDist2 < fBall2) {
						if (smx->iMark[pj]) continue;
						smx->iMark[pq->p] = 0;
						smx->iMark[pj] = 1;
						pq->fKey = fDist2;
						pq->p = pj;
						pq->ax = sx - x;
						pq->ay = sy - y;
						pq->az = sz - z;
						PQ_REPLACE(pq);
						fBall2 = pq->fKey;
						}
					}
				}
		GetNextCell:
			SETNEXT(cp);
			if (cp == ct) break;
			}
		cell = PARENT(cell);
		}
	smx->pqHead = pq;
	}


void smDensityInit(SMX smx,int bPeriodic)
{
	KDN *c;
	PINIT *p;
	PQ *pq,*pqLast;
	int cell;
	int pi,pin,pj,pNext,nSmooth,bLeap;
	float dx,dy,dz,x,y,z,h2,ih2,r2,rs,fNorm,ax,ay,az;
	KDN kdRoot;
	int ix,iy,iz,j,ppi;
	float fDist2;

	if (smx->kd->bOutDiag) puts(">> smDensityInit()");
	fflush(stdout);
	p = smx->kd->pInit;
	c = smx->kd->kdNodes;
	nSmooth = smx->nSmooth;
	pqLast = &smx->pq[smx->nSmooth-1];
	for (pi=0;pi<smx->kd->nInitActive;++pi) {
		p[pi].fBall2 = -1.0;
		p[pi].fDensity = 0.0;
		}
	/*
	 ** If there was no tree then the desities can remain 0.0.
	 */
	if (!c) return;
	/*
	 ** Clear Mark array.
	 */
	for (pi=0;pi<smx->kd->nInitActive;++pi) smx->iMark[pi] = 0;
	/*
	 ** Initialize Priority Queue.
	 */
	pNext = 0;
	bLeap = 1;
	while (1) {
		if (bLeap) {
			/*
			 ** Find next particle which is not done, and load the
			 ** priority queue with nSmooth number of particles.
			 */
			while (p[pNext].fBall2 >= 0) {
				++pNext;
				/*
				 ** Check if we are really finished.
				 */
				if (pNext == smx->kd->nInitActive) goto DoneDensity;
				}
			pi = pNext;
			x = p[pi].r[0];
			y = p[pi].r[1];
			z = p[pi].r[2];
			/*
			 ** First find the "local" Bucket.
			 ** This could mearly be the closest bucket to ri[3].
			 */
			cell = ROOT;
			while (c[cell].iDim >= 0) {
				if (p[pi].r[c[cell].iDim] < c[cell].fSplit)
					cell = LOWER(cell);
				else
					cell = UPPER(cell);
				}
			/*
			 ** Remove everything from the queue.
			 */
			smx->pqHead = NULL;
			for (pq=smx->pq;pq<=pqLast;++pq) smx->iMark[pq->p] = 0;
			/*
			 ** Add everything from pj up to and including pj+nSmooth-1.
			 */
			pj = c[cell].pLower;
			if (pj > smx->kd->nInitActive - nSmooth)
				pj = smx->kd->nInitActive - nSmooth;
			for (pq=smx->pq;pq<=pqLast;++pq) {
				smx->iMark[pj] = 1;
				dx = x - p[pj].r[0];
				dy = y - p[pj].r[1];
				dz = z - p[pj].r[2];
				pq->fKey = dx*dx + dy*dy + dz*dz;
				pq->p = pj++;
				pq->ax = 0.0;
				pq->ay = 0.0;
				pq->az = 0.0;
				}
			PQ_BUILD(smx->pq,nSmooth,smx->pqHead);
			}
		else {
			/*
			 ** Calculate the priority queue using the previous particles!
			 */
			pi = pin;
			x = p[pi].r[0];
			y = p[pi].r[1];
			z = p[pi].r[2];
			smx->pqHead = NULL;
			for (pq=smx->pq;pq<=pqLast;++pq) {
				pq->ax -= ax;
				pq->ay -= ay;
				pq->az -= az;
				dx = x + pq->ax - p[pq->p].r[0];
				dy = y + pq->ay - p[pq->p].r[1];
				dz = z + pq->az - p[pq->p].r[2];
				pq->fKey = dx*dx + dy*dy + dz*dz;
				}
			PQ_BUILD(smx->pq,nSmooth,smx->pqHead);
			ax = 0.0;
			ay = 0.0;
			az = 0.0;
			}
		smBallSearch(smx,smx->pqHead->fKey,p[pi].r);
		p[pi].fBall2 = smx->pqHead->fKey;
		/*
		 ** Pick next particle, 'pin'.
		 ** Create fList and pList for function 'fncSmooth'.
		 */
		bLeap = 1;
		h2 = smx->pqHead->fKey;
		ih2 = 4.0/h2;
		fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
		for (pq=smx->pq;pq<=pqLast;++pq) {
			if (pq == smx->pqHead) continue;
			/*
			 ** Density loop
			 */
			r2 = pq->fKey*ih2;
			rs = 2.0 - sqrt(r2);
			if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
			else rs = 0.25*rs*rs*rs;
			rs *= fNorm;
			p[pi].fDensity += rs*p[pq->p].fMass;
			p[pq->p].fDensity += rs*p[pi].fMass;
			/*
			 ** Pick next particle.
			 */
			if (p[pq->p].fBall2 >= 0) continue;
			if (pq->fKey < h2) {
				pin = pq->p;
				bLeap = 0;
				h2 = pq->fKey;
				ax = pq->ax;
				ay = pq->ay;
				az = pq->az;
				}
			}
		}
 DoneDensity:
	/*
	 ** Want to make a replica scatterers if the simualtion is periodic.
	 ** Set up a cell for the entire periodic box.
	 */
	for (j=0;j<3;++j) {
		kdRoot.bnd.fMin[j] = smx->kd->fCenter[j] - 0.5*smx->kd->fPeriod[j];
		kdRoot.bnd.fMax[j] = smx->kd->fCenter[j] + 0.5*smx->kd->fPeriod[j];
		}
	smx->nExtraScat = 0;
	if (bPeriodic) {
		for (pi=0;pi<smx->kd->nInitActive;++pi) {
			for (ix=-1;ix<=1;++ix) {
				x = p[pi].r[0] + ix*smx->kd->fPeriod[0];
				for (iy=-1;iy<=1;++iy) {
					y = p[pi].r[1] + iy*smx->kd->fPeriod[1];
					for (iz=-1;iz<=1;++iz) {
						z = p[pi].r[2] + iz*smx->kd->fPeriod[2];
						if (ix || iy || iz) {
							INTERSECTNP(&kdRoot,x,y,z,fDist2);
							if (fDist2 < p[pi].fBall2) ++smx->nExtraScat;
							}
						} /* of iz */
					} /* of iy */
				} /* of ix */
			} /* of pi-loop */
		printf("nExtraScat:%d\n",smx->nExtraScat);
		/*
		 ** Now allocate the storage for the extra scatterers.
		 */
		smx->pp = malloc(smx->nExtraScat*sizeof(PINIT));
		assert(smx->pp != NULL);
		ppi = 0;
		for (pi=0;pi<smx->kd->nInitActive;++pi) {
			for (ix=-1;ix<=1;++ix) {
				x = p[pi].r[0] + ix*smx->kd->fPeriod[0];
				for (iy=-1;iy<=1;++iy) {
					y = p[pi].r[1] + iy*smx->kd->fPeriod[1];
					for (iz=-1;iz<=1;++iz) {
						z = p[pi].r[2] + iz*smx->kd->fPeriod[2];
						if (ix || iy || iz) {
							INTERSECTNP(&kdRoot,x,y,z,fDist2);
							if (fDist2 < p[pi].fBall2) {
								smx->pp[ppi] = p[pi];
								smx->pp[ppi].r[0] = x;
								smx->pp[ppi].r[1] = y;
								smx->pp[ppi].r[2] = z;
								++ppi;
								}
							}
						} /* of iz */
					} /* of iy */
				} /* of ix */
			} /* of pi-loop */
		}
	if (smx->kd->bOutDiag) puts("<< smDensityInit()");
	fflush(stdout);
	}


int smBallGather(SMX smx,float fBall2,float *ri,NN *nnList)
{
	KDN *c;
	PMOVE *p;
	int pj,nCnt,cp;
	float dx,dy,dz,x,y,z,fDist2;

	c = smx->kd->kdNodes;
	if (!c) return(0);
	p = smx->kd->pMove;
	x = ri[0];
	y = ri[1];
	z = ri[2];
	nCnt = 0;
	cp = ROOT;
	while (1) {
		INTERSECTNP(&c[cp],x,y,z,fDist2);
		if (fDist2 < fBall2) {
			/*
			 ** We have an intersection to test.
			 */
			if (c[cp].iDim >= 0) {
				cp = LOWER(cp);
				continue;
				}
			else {
				for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
					dx = x - p[pj].r[0];
					dy = y - p[pj].r[1];
					dz = z - p[pj].r[2];
					fDist2 = dx*dx + dy*dy + dz*dz;
					if (fDist2 < fBall2) {
						nnList[nCnt].p = pj;
						nnList[nCnt].fDist2 = fDist2;
						nnList[nCnt].dx = dx;
						nnList[nCnt].dy = dy;
						nnList[nCnt].dz = dz;
						++nCnt;
						}
					}
				}
			}
		SETNEXT(cp);
		if (cp == ROOT) break;
		}
	return(nCnt);
	}


int ScatterCut(PINIT *p,int n,float fScatDens)
{
	PINIT t;
	int i,j;
	
	i = 0;
	j = n-1;
	while (1) {
		while (p[i].fDensity >= fScatDens)
			if (++i > j) goto done;
		while (p[j].fDensity < fScatDens)
			if (i > --j) goto done;
		t = p[i];
		p[i] = p[j];
		p[j] = t;
		}
 done:
	return(i);
	}


int smAccDensity(SMX smx,int bInitial)
{
	PINIT *p;
	PMOVE *pm;
	NN *nnList;
	int pi,nSmooth,i,pmi,nListSize;
	float ih2,fNorm,r2,rs;
	float fScatDens;

	if (smx->kd->bOutDiag) puts(">> smAccDensity()");
	fflush(stdout);
	/*
	 ** Allocate Nearest-Neighbor list.
	 */
	nListSize = smx->kd->nActive;
	nnList = (NN *)malloc(nListSize*sizeof(NN));
	assert(nnList != NULL);
	pm = smx->kd->pMove;
	/*
	 ** Clear accelerations for the moved positions.
	 */
	for (pmi=0;pmi<smx->kd->nActive;++pmi) {
		pm[pmi].a[0] = 0.0;
		pm[pmi].a[1] = 0.0;
		pm[pmi].a[2] = 0.0;
		}
	p = smx->kd->pInit;
	fScatDens = 0.0;
	for (pi=0;pi<smx->kd->nInitActive;++pi) {
		/*
		 ** Do a Ball Gather at the radius of the most distant particle
		 ** which smDensity sets in p[pi].fBall2.
		 */
		nSmooth = smBallGather(smx,p[pi].fBall2,p[pi].r,nnList);
		if (nSmooth > 0) {
			/*
			 ** Calculate the density acceleration by scattering the 
			 ** unmoved particle's kernel to all the moved positions.
			 */
			ih2 = 4.0/p[pi].fBall2;
			fNorm = M_1_PI*ih2*ih2*sqrt(ih2)*p[pi].fMass;
			for (i=0;i<nSmooth;++i) {
				r2 = nnList[i].fDist2*ih2;
				rs = sqrt(r2);
				if (r2 < 1.0) rs = -3.0 + 2.25*rs;
				else rs = -3.0/rs + 3.0 - 0.75*rs;
				rs *= fNorm;
				pmi = nnList[i].p;
				pm[pmi].a[0] += nnList[i].dx*rs;
				pm[pmi].a[1] += nnList[i].dy*rs;
				pm[pmi].a[2] += nnList[i].dz*rs;
				}
			if (fScatDens == 0.0) fScatDens = p[pi].fDensity;
			else if (p[pi].fDensity < fScatDens) fScatDens = p[pi].fDensity;
			}
		else if (bInitial) {
			/*
			 ** On the initial iteration we want to remove all potential
			 ** high denity scatterers which did not scatter. Note: we
			 ** should only do this in the dark matter only case.
			 */
			p[pi].fDensity = 0.0;
			}
		}
	p = smx->pp;
	for (pi=0;pi<smx->nExtraScat;++pi) {
		/*
		 ** Do a Ball Gather at the radius of the most distant particle
		 ** which smDensity sets in pp[pi].fBall2.
		 */
		nSmooth = smBallGather(smx,p[pi].fBall2,p[pi].r,nnList);
		if (nSmooth > 0) {
			/*
			 ** Calculate the density acceleration by scattering the 
			 ** unmoved particle's kernel to all the moved positions.
			 */
			ih2 = 4.0/p[pi].fBall2;
			fNorm = M_1_PI*ih2*ih2*sqrt(ih2)*p[pi].fMass;
			for (i=0;i<nSmooth;++i) {
				r2 = nnList[i].fDist2*ih2;
				rs = sqrt(r2);
				if (r2 < 1.0) rs = -3.0 + 2.25*rs;
				else rs = -3.0/rs + 3.0 - 0.75*rs;
				rs *= fNorm;
				pmi = nnList[i].p;
				pm[pmi].a[0] += nnList[i].dx*rs;
				pm[pmi].a[1] += nnList[i].dy*rs;
				pm[pmi].a[2] += nnList[i].dz*rs;
				}
			if (fScatDens == 0.0) fScatDens = p[pi].fDensity;
			else if (p[pi].fDensity < fScatDens) fScatDens = p[pi].fDensity;
			}
		else if (bInitial) {
			/*
			 ** On the initial iteration we want to remove all potential
			 ** high denity scatterers which did not scatter. Note: we
			 ** should only do this in the dark matter only case.
			 */
			p[pi].fDensity = 0.0;
			}
		}
	smx->kd->nInitActive = ScatterCut(smx->kd->pInit,
									  smx->kd->nInitActive,fScatDens);
	smx->nExtraScat = ScatterCut(smx->pp,smx->nExtraScat,fScatDens);
	free(nnList);
	if (smx->kd->bOutDiag) puts("<< smAccDensity()");
	fflush(stdout);
	return(smx->kd->nInitActive + smx->nExtraScat);
 	}







