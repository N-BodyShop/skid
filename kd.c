#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
#include <limits.h>
#include "kd.h"
#include "grav.h"
#include "tipsydefs.h"


void kdTime(KD kd,int *puSecond,int *puMicro)
{
	struct rusage ru;

	getrusage(0,&ru);
	*puMicro = ru.ru_utime.tv_usec - kd->uMicro;
	*puSecond = ru.ru_utime.tv_sec - kd->uSecond;
	if (*puMicro < 0) {
		*puMicro += 1000000;
		*puSecond -= 1;
		}
	kd->uSecond = ru.ru_utime.tv_sec;
	kd->uMicro = ru.ru_utime.tv_usec;
	}


int kdInit(KD *pkd,int nBucket,float *fPeriod,float *fCenter)
{
	KD kd;
	int j;

	kd = (KD)malloc(sizeof(struct kdContext));
	assert(kd != NULL);
	kd->nBucket = nBucket;
	for (j=0;j<3;++j) {
		kd->fPeriod[j] = fPeriod[j];
		kd->fCenter[j] = fCenter[j];
		}
	kd->pMove = NULL;
	kd->pInit = NULL;
	kd->pGroup = NULL;
	kd->kdNodes = NULL;
	kd->kdGroup = NULL;
	kd->piGroup = NULL;
	*pkd = kd;
	return(1);
	}


void kdSetUniverse(KD kd,float G,float Omega0,float H0,float z)
{
	kd->G = G;
	kd->Omega0 = Omega0;
	kd->H0 = H0;
	kd->z = z;
	}


void kdSetSoft(KD kd,float fEps)
{
	int i;
	
	for (i=0;i<kd->nParticles;++i) {
		kd->pInit[i].fSoft = fEps;
		}
	}


int kdParticleType(KD kd,int iOrder)
{
	if (iOrder < kd->nGas) return(GAS);
	else if (iOrder < (kd->nGas + kd->nDark)) return(DARK);
	else if (iOrder < kd->nParticles) return(STAR);
	else return(0);
	}


int kdReadTipsy(KD kd,FILE *fp)
{
	PINIT *p;
	int i,j;
	struct dump h;
	struct gas_particle gp;
	struct dark_particle dp;
	struct star_particle sp;

	fread(&h,sizeof(struct dump),1,fp);
	kd->inType = 0;
	kd->nDark = h.ndark;
	if (kd->nDark) kd->inType |= DARK;
	kd->nGas = h.nsph;
	if (kd->nGas) kd->inType |= GAS;
	kd->nStar = h.nstar;
	if (kd->nStar) kd->inType |= STAR;
	kd->fTime = h.time;
	kd->nParticles = kd->nDark + kd->nGas + kd->nStar;
	kd->nInitActive = kd->nParticles;
	/*
	 ** Allocate arrays.
	 */
	kd->pInit = (PINIT *)malloc(kd->nParticles*sizeof(PINIT));
	assert(kd->pInit != NULL);
	printf("nDark:%d nGas:%d nStar:%d\n",kd->nDark,kd->nGas,kd->nStar);
	fflush(stdout);
	p = kd->pInit;
	/*
	 ** Read Stuff!
	 */
	for (i=0;i<kd->nParticles;++i) {
		p[i].iOrder = i;
		p[i].fDensity = 0.0;
		switch (kdParticleType(kd,i)) {
		case (GAS):
			fread(&gp,sizeof(struct gas_particle),1,fp);
			p[i].fMass = gp.mass;
			p[i].fSoft = gp.hsmooth;
			p[i].fTemp = gp.temp;
			for (j=0;j<3;++j) {
				p[i].r[j] = gp.pos[j];
				p[i].v[j] = gp.vel[j];
				}
			break;
		case (DARK):
			fread(&dp,sizeof(struct dark_particle),1,fp);
			p[i].fMass = dp.mass;
			p[i].fSoft = dp.eps;
			p[i].fTemp = 0.0;
			for (j=0;j<3;++j) {
				p[i].r[j] = dp.pos[j];
				p[i].v[j] = dp.vel[j];
				}
			break;
		case (STAR):
			fread(&sp,sizeof(struct star_particle),1,fp);
			p[i].fMass = sp.mass;
			p[i].fSoft = sp.eps;
			p[i].fTemp = 0.0;
			for (j=0;j<3;++j) {
				p[i].r[j] = sp.pos[j];
				p[i].v[j] = sp.vel[j];
				}
			break;
			}
		}
	return(kd->nParticles);
	}


void kdSelectInit(KD kd,int d,int k,int l,int r)
{
	PINIT *p,t;
	double v;
	int i,j;

	p = kd->pInit;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


void kdSelectMove(KD kd,int d,int k,int l,int r)
{
	PMOVE *p,t;
	double v;
	int i,j;

	p = kd->pMove;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


void Combine(KDN *p1,KDN *p2,KDN *pOut)
{
	int j;

	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];
		}
	}


void UpPassInit(KD kd,int iCell)
{
	KDN *c;
	int l,u,pj,j;

	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWER(iCell);
		u = UPPER(iCell);
		UpPassInit(kd,l);
		UpPassInit(kd,u);
		Combine(&c[l],&c[u],&c[iCell]);
		}
	else {
		l = c[iCell].pLower;
		u = c[iCell].pUpper;
		for (j=0;j<3;++j) {
			c[iCell].bnd.fMin[j] = kd->pInit[u].r[j];
			c[iCell].bnd.fMax[j] = kd->pInit[u].r[j];
			}
		for (pj=l;pj<u;++pj) {
			for (j=0;j<3;++j) {
				if (kd->pInit[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = kd->pInit[pj].r[j];
				if (kd->pInit[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = kd->pInit[pj].r[j];
				}
			}
		}
	}


void UpPassMove(KD kd,int iCell)
{
	KDN *c;
	int l,u,pj,j;

	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWER(iCell);
		u = UPPER(iCell);
		UpPassMove(kd,l);
		UpPassMove(kd,u);
		Combine(&c[l],&c[u],&c[iCell]);
		}
	else {
		l = c[iCell].pLower;
		u = c[iCell].pUpper;
		for (j=0;j<3;++j) {
			c[iCell].bnd.fMin[j] = kd->pMove[u].r[j];
			c[iCell].bnd.fMax[j] = kd->pMove[u].r[j];
			}
		for (pj=l;pj<u;++pj) {
			for (j=0;j<3;++j) {
				if (kd->pMove[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = kd->pMove[pj].r[j];
				if (kd->pMove[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = kd->pMove[pj].r[j];
				}
			}
		}
	}


int kdBuildTree(KD kd)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	BND bnd;

	if (kd->nInitActive == 0) {
		if (kd->kdNodes) free(kd->kdNodes);
		kd->kdNodes = NULL;
		return(1);
		}
	assert(kd->nInitActive > 0);
	n = kd->nInitActive;
	kd->nLevels = 1;
	l = 1;
	while (n > kd->nBucket) {
		n = n>>1;
		l = l<<1;
		++kd->nLevels;
		}
	kd->nSplit = l;
	kd->nNodes = l<<1;
	if (kd->kdNodes) {
		free(kd->kdNodes);
		kd->kdNodes = NULL;
		}
	kd->kdNodes = (KDN *)malloc(kd->nNodes*sizeof(KDN));
	assert(kd->kdNodes != NULL);
	/*
	 ** Calculate Bounds.
	 */
	for (j=0;j<3;++j) {
		bnd.fMin[j] = kd->pInit[0].r[j];
		bnd.fMax[j] = kd->pInit[0].r[j];
		}
	for (i=1;i<kd->nInitActive;++i) {
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->pInit[i].r[j]) 
				bnd.fMin[j] = kd->pInit[i].r[j];
			else if (bnd.fMax[j] < kd->pInit[i].r[j])
				bnd.fMax[j] = kd->pInit[i].r[j];
			}
		}
	/*
	 ** Set up ROOT node
	 */
	c = kd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = kd->nInitActive-1;
	c[ROOT].bnd = bnd;
	i = ROOT;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			kdSelectInit(kd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = kd->pInit[m].r[d];
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m-1;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXT(i);
			if (i == ROOT) break;
			}
		}
	UpPassInit(kd,ROOT);
	return(1);
	}


int kdBuildMoveTree(KD kd)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	BND bnd;

	if (kd->nActive == 0) {
		if (kd->kdNodes) free(kd->kdNodes);
		kd->kdNodes = NULL;
		return(1);
		}
	assert(kd->nActive > 0);
	n = kd->nActive;
	kd->nLevels = 1;
	l = 1;
	while (n > kd->nBucket) {
		n = n>>1;
		l = l<<1;
		++kd->nLevels;
		}
	kd->nSplit = l;
	kd->nNodes = l<<1;
	if (kd->kdNodes) {
		free(kd->kdNodes);
		kd->kdNodes = NULL;
		}
	kd->kdNodes = (KDN *)malloc(kd->nNodes*sizeof(KDN));
	assert(kd->kdNodes != NULL);
	/*
	 ** Calculate Bounds.
	 */
	for (j=0;j<3;++j) {
		bnd.fMin[j] = kd->pMove[0].r[j];
		bnd.fMax[j] = kd->pMove[0].r[j];
		}
	for (i=1;i<kd->nMove;++i) {
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->pMove[i].r[j]) 
				bnd.fMin[j] = kd->pMove[i].r[j];
			else if (bnd.fMax[j] < kd->pMove[i].r[j])
				bnd.fMax[j] = kd->pMove[i].r[j];
			}
		}
	/*
	 ** Set up ROOT node
	 */
	c = kd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = kd->nActive-1;
	c[ROOT].bnd = bnd;
	i = ROOT;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			kdSelectMove(kd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = kd->pMove[m].r[d];
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m-1;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXT(i);
			if (i == ROOT) break;
			}
		}
	UpPassMove(kd,ROOT);
	return(1);
	}


int CutCriterion(KD kd,int pi,float fDensMin,float fTempMax,
				 float fMassMax,int bGasAndDark)
{
	if (kd->pInit[pi].fMass > fMassMax) return(0);
	switch (kd->inType) {
	case (DARK):
		if (kd->pInit[pi].fDensity >= fDensMin) 
			return(1);
		break;
	case (GAS):
	case (DARK|GAS):
		if (bGasAndDark) {
			if (kdParticleType(kd,kd->pInit[pi].iOrder) == DARK) {
				if (kd->pInit[pi].fDensity >= fDensMin) return(1);
				}
			}
		if (kdParticleType(kd,kd->pInit[pi].iOrder) == GAS) {
			if (kd->pInit[pi].fDensity >= fDensMin && 
				kd->pInit[pi].fTemp <= fTempMax) return(1);
			}
		break;
	case (STAR):
	case (DARK|STAR):
		if (kdParticleType(kd,kd->pInit[pi].iOrder) == STAR) return(1);
		break;
	case (GAS|STAR):
	case (DARK|GAS|STAR):
		if (kdParticleType(kd,kd->pInit[pi].iOrder) == GAS) {
			if (kd->pInit[pi].fDensity >= fDensMin && 
				kd->pInit[pi].fTemp <= fTempMax) return(1);
			}
		else if (kdParticleType(kd,kd->pInit[pi].iOrder) == STAR) {
			return(1);
			}
		}
	return(0);
	}


int ScatterCriterion(KD kd,int pi,int bGasAndDark)
{
	switch (kd->inType) {
	case (DARK):
		return(1);
	case (GAS):
	case (DARK|GAS):
		if (bGasAndDark) return(1);
		else if (kdParticleType(kd,kd->pInit[pi].iOrder) == GAS) {
			return(1);
			}
		break;
	case (STAR):
	case (DARK|STAR):
		if (kdParticleType(kd,kd->pInit[pi].iOrder) == STAR) return(1);
		break;
	case (GAS|STAR):
	case (DARK|GAS|STAR):
		if (kdParticleType(kd,kd->pInit[pi].iOrder) == GAS) {
			return(1);
			}
		else if (kdParticleType(kd,kd->pInit[pi].iOrder) == STAR) {
			return(1);
			}
		}
	return(0);
	}


void kdInitMove(KD kd,float fDensMin,float fTempMax,float fMassMax,
				float fCvg,int bGasAndDark)
{
	int pi,nCnt,j;
	float fCvg2;

	/*
	 ** First count number of particles meating the criterion.
	 */
	kd->nMove = 0;
	for (pi=0;pi<kd->nParticles;++pi) {
		kd->nMove += CutCriterion(kd,pi,fDensMin,fTempMax,
								  fMassMax,bGasAndDark);
		}
	kd->nActive = kd->nMove;
	printf("Number of Moving particles:%d\n",kd->nMove);
	fflush(stdout);
	/*
	 ** Allocate moving particles and initialize.
	 */
	kd->pMove = (PMOVE *)malloc(kd->nMove*sizeof(PMOVE));
	assert(kd->pMove != NULL);
	nCnt = 0;
	for (pi=0;pi<kd->nParticles;++pi) {
		if (CutCriterion(kd,pi,fDensMin,fTempMax,fMassMax,bGasAndDark)) {
			for (j=0;j<3;++j) {
				kd->pMove[nCnt].r[j] = kd->pInit[pi].r[j];
				kd->pMove[nCnt].rOld[j] = kd->pInit[pi].r[j];
				}
			kd->pMove[nCnt].iOrder = kd->pInit[pi].iOrder;
			++nCnt;
			}
		}
	/*
	 ** Change the fBall2 of all the unmoved particles to have a minimum
	 ** hSmooth of fCvg.
	 */
	fCvg2 = fCvg*fCvg;
	for (pi=0;pi<kd->nParticles;++pi) {
		if (kd->pInit[pi].fBall2 < fCvg2) kd->pInit[pi].fBall2 = fCvg2;
		}
	}


int kdScatterActive(KD kd,int bGasAndDark)
{
	PINIT *p,t;
	int i,j;
	
	p = kd->pInit;
	i = 0;
	j = kd->nParticles-1;
	while (1) {
		while (ScatterCriterion(kd,i,bGasAndDark))
			if (++i > j) goto done;
		while (!ScatterCriterion(kd,j,bGasAndDark))
			if (i > --j) goto done;
		t = p[i];
		p[i] = p[j];
		p[j] = t;
		}
 done:
	kd->nInitActive = i;
	return(i);
	}


int kdScatterCut(KD kd)
{
	PINIT *p,t;
	int i,j;
	
	if(kd->inType == DARK) {
		p = kd->pInit;
		i = 0;
		j = kd->nInitActive-1;
		while (1) {
			while (p[i].fDensity > 0.0)
				if (++i > j) goto done;
			while (p[j].fDensity < 0.0)
				if (i > --j) goto done;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			}
	done:
		kd->nInitActive = i;
		}
	return(kd->nInitActive);
	}


void kdMoveParticles(KD kd,float fStep)
{
	PMOVE *p;
	int i,j;
	float ax,ay,az,ai;

	p = kd->pMove;
	for (i=0;i<kd->nActive;++i) {
		ax = p[i].a[0];
		ay = p[i].a[1];
		az = p[i].a[2];
		ai = sqrt(ax*ax + ay*ay + az*az);
		if(ai > 0.0)
		    ai = fStep/sqrt(ax*ax + ay*ay + az*az);
		else
		    ai = 0.0;
		p[i].r[0] -= ai*ax;
		p[i].r[1] -= ai*ay;
		p[i].r[2] -= ai*az;
		for (j=0;j<3;++j) {
			if (p[i].r[j] > kd->fCenter[j]+0.5*kd->fPeriod[j])
				p[i].r[j] -= kd->fPeriod[j];
			if (p[i].r[j] <= kd->fCenter[j]-0.5*kd->fPeriod[j])
				p[i].r[j] += kd->fPeriod[j];
			}
		}
	}


int kdPruneInactive(KD kd,float fCvg)
{
	PMOVE *p,t;
	int i,j;
	float dx,dy,dz,dr2,fCvg2,hx,hy,hz;

	p = kd->pMove;
	hx = 0.5*kd->fPeriod[0];
	hy = 0.5*kd->fPeriod[1];
	hz = 0.5*kd->fPeriod[2];
	i = 0;
	j = kd->nActive-1;
	fCvg2 = fCvg*fCvg;
	while (1) {
		while (1) {
			dx = p[i].r[0] - p[i].rOld[0];
			dy = p[i].r[1] - p[i].rOld[1];
			dz = p[i].r[2] - p[i].rOld[2];
			if (dx > hx) dx -= 2*hx;
			if (dx <= -hx) dx += 2*hx;
			if (dy > hy) dy -= 2*hy;
			if (dy <= -hy) dy += 2*hy;
			if (dz > hz) dz -= 2*hz;
			if (dz <= -hz) dz += 2*hz;
			dr2 = dx*dx + dy*dy + dz*dz;
			if (dr2 < fCvg2) break;
			if (++i > j) goto done;
			}
		while (1) {
			dx = p[j].r[0] - p[j].rOld[0];
			dy = p[j].r[1] - p[j].rOld[1];
			dz = p[j].r[2] - p[j].rOld[2];
			if (dx > hx) dx -= 2*hx;
			if (dx <= -hx) dx += 2*hx;
			if (dy > hy) dy -= 2*hy;
			if (dy <= -hy) dy += 2*hy;
			if (dz > hz) dz -= 2*hz;
			if (dz <= -hz) dz += 2*hz;
			dr2 = dx*dx + dy*dy + dz*dz;
			if (dr2 >= fCvg2) break;
			if (i > --j) goto done;
			}
		t = p[i];
		p[i] = p[j];
		p[j] = t;
		}
 done:
	kd->nActive = i;
	for (i=0;i<kd->nActive;++i) {
		for (j=0;j<3;++j) {
			p[i].rOld[j] = p[i].r[j];
			}
		}
	return(i);
	}


void kdReactivateMove(KD kd)
{
    kd->nActive = kd->nMove;
	}


void kdFoF(KD kd,float fTau)
{
	PMOVE *p;
	KDN *c;
	int pi,pj,pn,cp;
	int *Group,iGroup;
	int *Fifo,iHead,iTail,nFifo;
	float fTau2;
	float dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;

    kd->nActive = kd->nMove;
	if (kd->nActive == 0) return;
	kdBuildMoveTree(kd);
	p = kd->pMove;
	c = kd->kdNodes;
	lx = kd->fPeriod[0];
	ly = kd->fPeriod[1];
	lz = kd->fPeriod[2];
	fTau2 = fTau*fTau;
	Group = (int *)malloc(kd->nActive*sizeof(int));
	assert(Group != NULL);
	for (pn=0;pn<kd->nActive;++pn) Group[pn] = 0;
	nFifo = kd->nActive;
	Fifo = (int *)malloc(nFifo*sizeof(int));
	assert(Fifo != NULL);
	iHead = 0;
	iTail = 0;
	iGroup = 0;
	for (pn=0;pn<kd->nActive;++pn) {
		if (Group[pn]) continue;
		++iGroup;
		/*
		 ** Mark it and add to the do-fifo.
		 */
		Group[pn] = iGroup;
		Fifo[iTail++] = pn;
		if (iTail == nFifo) iTail = 0;
		while (iHead != iTail) {
			pi = Fifo[iHead++];
			if (iHead == nFifo) iHead = 0;
			/*
			 ** Now do an fEps-Ball Gather!
			 */
			x = p[pi].r[0];
			y = p[pi].r[1];
			z = p[pi].r[2];
			cp = ROOT;
			while (1) {
				INTERCONT(c,cp,fTau2,lx,ly,lz,x,y,z,sx,sy,sz);
				/*
				 ** We have an intersection to test.
				 */
				if (c[cp].iDim >= 0) {
					cp = LOWER(cp);
					continue;
					}
				else {
					for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
						if (Group[pj]) continue;
						dx = sx - p[pj].r[0];
						dy = sy - p[pj].r[1];
						dz = sz - p[pj].r[2];
						fDist2 = dx*dx + dy*dy + dz*dz;
						if (fDist2 < fTau2) {
							/*
							 ** Mark it and add to the do-fifo.
							 */
							Group[pj] = iGroup;
							Fifo[iTail++] = pj;
							if (iTail == nFifo) iTail = 0;
							}
						}
					SETNEXT(cp);
					if (cp == ROOT) break;
					continue;
					}
			ContainedCell:
				for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
					if (Group[pj]) continue;
					/*
					 ** Mark it and add to the do-fifo.
					 */
					Group[pj] = iGroup;
					Fifo[iTail++] = pj;
					if (iTail == nFifo) iTail = 0;
					}
			GetNextCell:
				SETNEXT(cp);
				if (cp == ROOT) break;
				}
			}
		}
	free(Fifo);
	kd->nGroup = iGroup+1;
	kd->piGroup = (int *)malloc(kd->nParticles*sizeof(int));
	assert(kd->piGroup);
	for (pi=0;pi<kd->nParticles;++pi) {
		kd->piGroup[pi] = 0;
		}
	for (pi=0;pi<kd->nActive;++pi) {
		kd->piGroup[p[pi].iOrder] = Group[pi];
		}
	free(Group);
	}


void kdInGroup(KD kd,char *pszIn)
{
	FILE *fp;
	int nGroup,iGroup,n,pi;

	kd->piGroup = (int *)malloc(kd->nParticles*sizeof(int));
	assert(kd->piGroup);
	for (pi=0;pi<kd->nParticles;++pi) {
		kd->piGroup[pi] = 0;
		}
	if (pszIn) {
		fp = fopen(pszIn,"r");
		if (!fp) {
			fprintf(stderr,"ERROR: Could not open file:%s\n",pszIn);
			exit(1);
			}
		}
	else {
		fp = stdin;
		}
	fscanf(fp,"%d",&n);
	/*
	 ** Sanity Check
	 */
	if (n != kd->nParticles) {
		fprintf(stderr,"ERROR: Mismatched number of particles\n");
		fprintf(stderr,"Number in Group file %s: %d\n",pszIn,n);
		fprintf(stderr,"Number in TIPSY BINARY input file: %d\n",
				kd->nParticles);
		exit(1);
		}
	nGroup = 0;
	for (pi=0;pi<kd->nParticles;++pi) {
		fscanf(fp,"%d",&iGroup);
		kd->piGroup[pi] = iGroup;
		if (iGroup > nGroup) nGroup = iGroup;
		}
	kd->nGroup = nGroup+1;
	}


/*
 ** This function sets the relative coordinate for each group. 
 ** The relative coordinate can be the coordinates of any group member.
 */
void kdInitpGroup(KD kd)
{
	PINIT *p;
	int i,iGroup,j;

	kd->pGroup = (PGROUP *)malloc(kd->nGroup*sizeof(PGROUP));
	assert(kd->pGroup != NULL);
	for (i=0;i<kd->nGroup;++i) {
		kd->pGroup[i].nMembers = 0;
		kd->pGroup[i].fMass = 0.0;
		kd->pGroup[i].fRadius = 0.0;
		for (j=0;j<3;++j) {
			kd->pGroup[i].vcm[j] = 0.0;
			kd->pGroup[i].rCenter[j] = 0.0;
			kd->pGroup[i].rel[j] = 0.0;
			}
		}
	p = kd->pInit;
	for (i=0;i<kd->nParticles;++i) {
		iGroup = kd->piGroup[p[i].iOrder];
		if (!iGroup) continue;
		kd->pGroup[iGroup].fMass += p[i].fMass;
		kd->pGroup[iGroup].nMembers += 1;
		for (j=0;j<3;++j) {
			kd->pGroup[iGroup].rel[j] = p[i].r[j];
			}
		}
	}


/*
 ** Calculates the density "center" for all the groups.
 ** It uses the moving particle mean position to determine this.
 ** Must call kdInitpGroup() before calling this function.
 */
void kdCalcCenter(KD kd)
{
	PMOVE *p;
	PINIT *q;
	int i,j,iGroup;
	float del;

	for (i=0;i<kd->nGroup;++i) {
		for (j=0;j<3;++j) {
			kd->pGroup[i].rCenter[j] = 0.0;
			kd->pGroup[i].vcm[j] = 0.0;
			}
		}
	/*
	 ** Also need to calculate kd->pGroup[i].vcm[j] and 
	 ** to completely match the kdReadCenter() function.
	 */
	p = kd->pMove;
	for (i=0;i<kd->nMove;++i) {
		iGroup = kd->piGroup[p[i].iOrder];
		if (!iGroup) continue;
		for (j=0;j<3;++j) {
			del = p[i].r[j] - kd->pGroup[iGroup].rel[j];
			if (del > 0.5*kd->fPeriod[j]) del -= kd->fPeriod[j];
			if (del <= -0.5*kd->fPeriod[j]) del += kd->fPeriod[j];
			kd->pGroup[iGroup].rCenter[j] += del;
			}
		}
	q = kd->pInit;
	for (i=0;i<kd->nParticles;++i) {
		iGroup = kd->piGroup[q[i].iOrder];
		if (!iGroup) continue;
		for (j=0;j<3;++j) {
			kd->pGroup[iGroup].vcm[j] += q[i].fMass*q[i].v[j];
			}
		}
	for (i=1;i<kd->nGroup;++i) {
		for (j=0;j<3;++j) {
			kd->pGroup[i].rCenter[j] /= kd->pGroup[i].nMembers;
			kd->pGroup[i].rCenter[j] += kd->pGroup[i].rel[j];
			if (kd->pGroup[i].rCenter[j] > kd->fCenter[j]+0.5*kd->fPeriod[j])
				kd->pGroup[i].rCenter[j] -= kd->fPeriod[j];
			if (kd->pGroup[i].rCenter[j] <= kd->fCenter[j]-0.5*kd->fPeriod[j])
				kd->pGroup[i].rCenter[j] += kd->fPeriod[j];
			kd->pGroup[i].vcm[j] /= kd->pGroup[i].fMass;
			}
		}
	/*
	 ** Can free up the moving particles and move tree here.
	 */
	free(kd->pMove);
	kd->pMove = NULL;
	if (kd->kdNodes) {
		free(kd->kdNodes);
		kd->kdNodes = NULL;
		}
	}


/*
 ** This function reads in the density centers from the .gtp file if it 
 ** exists, otherwise it sets the center to be the center-of-mass.
 ** Requires that kdInitpGroup() has been called first!
 */
void kdReadCenter(KD kd,char *pszGtp)
{
	FILE *fp;
	PINIT *p;
	struct dump h;
	struct star_particle sp;
	int i,j,iGroup,bMismatch;
	float del;

	fp = fopen(pszGtp,"r");
	if (fp != NULL) {
		/*
		 ** Read from the file.
		 */
		fread(&h,sizeof(struct dump),1,fp);
		/*
		 ** Make sure that the number of Group-particles in this file
		 ** agrees with that which we set before in kdInGroup!
		 */
		if (h.nstar != kd->nGroup-1) {
			fprintf(stderr,"ERROR: Could not open file:%s\n",pszGtp);
			exit(1);
			}
		/*
		 ** Skip over any gas or dark matter particles that may be
		 ** in the file as well. There shouldn't be any though unless
		 ** future software uses a file interpretation of this sort.
		 */
		fseek(fp,h.nsph*sizeof(struct gas_particle),SEEK_CUR);
		fseek(fp,h.ndark*sizeof(struct dark_particle),SEEK_CUR);
		bMismatch = 0;
		for (i=0;i<h.nstar;++i) {
			fread(&sp,sizeof(struct star_particle),1,fp);
			/*
			 ** Verify that the mass of each group in the .gtp file
			 ** matches that calculated from the .grp file as a sort
			 ** of sanity check.
			 */
			if (!bMismatch) {
				if (sp.mass != kd->pGroup[i+1].fMass) {
					fprintf(stderr,"WARNING: Mismatch in group masses between\n");
					fprintf(stderr,"         the .grp and .gtp files.\n");
					}
				bMismatch = 1;
				}
			for (j=0;j<3;++j) {
				kd->pGroup[i+1].rCenter[j] = sp.pos[j];
				kd->pGroup[i+1].vcm[j] = sp.vel[j];
				}
			}
		}
	else {
		/*
		 ** Set the center to be the center-of-mass.
		 */
		for (i=0;i<kd->nGroup;++i) {
			for (j=0;j<3;++j) {
				kd->pGroup[i].rCenter[j] = 0.0;
				kd->pGroup[i].vcm[j] = 0.0;
				}
			}
		p = kd->pInit;
		for (i=0;i<kd->nParticles;++i) {
			iGroup = kd->piGroup[p[i].iOrder];
			if (!iGroup) continue;
			for (j=0;j<3;++j) {
				del = p[i].r[j] - kd->pGroup[iGroup].rel[j];
				if (del > 0.5*kd->fPeriod[j]) del -= kd->fPeriod[j];
				if (del <= -0.5*kd->fPeriod[j]) del += kd->fPeriod[j];
				kd->pGroup[iGroup].rCenter[j] += p[i].fMass*del;
				kd->pGroup[iGroup].vcm[j] += p[i].fMass*p[i].v[j];
				}
			}
		for (i=1;i<kd->nGroup;++i) {
			for (j=0;j<3;++j) {
				kd->pGroup[i].rCenter[j] /= kd->pGroup[i].fMass;
				kd->pGroup[i].rCenter[j] += kd->pGroup[i].rel[j];
				if (kd->pGroup[i].rCenter[j] > kd->fCenter[j]+0.5*kd->fPeriod[j])
					kd->pGroup[i].rCenter[j] -= kd->fPeriod[j];
				if (kd->pGroup[i].rCenter[j] <= kd->fCenter[j]-0.5*kd->fPeriod[j])
					kd->pGroup[i].rCenter[j] += kd->fPeriod[j];
				kd->pGroup[i].vcm[j] /= kd->pGroup[i].fMass;
				}
			}
		}
	}


void kdTooSmall(KD kd,int nMembers)
{
	int *pnMembers,*pMap;
	int i,pi,nGroup;

	pnMembers = (int *)malloc(kd->nGroup*sizeof(int));
	assert(pnMembers != NULL);
	pMap = (int *)malloc(kd->nGroup*sizeof(int));
	assert(pMap != NULL);
	for (i=0;i<kd->nGroup;++i) pnMembers[i] = 0;
	for (pi=0;pi<kd->nParticles;++pi) {
		++pnMembers[kd->piGroup[pi]];
		}
	for (i=1;i<kd->nGroup;++i) {
		if (pnMembers[i] < nMembers) {
			pnMembers[i] = 0;
			}
		}
	/*
	 ** Create a remapping!
	 */
	pMap[0] = 0;
	nGroup = 1;
	for (i=1;i<kd->nGroup;++i) {
		pMap[i] = nGroup;
		kd->pGroup[nGroup] = kd->pGroup[i];
		kd->kdGroup[nGroup] = kd->kdGroup[i];
		if (pnMembers[i] == 0) {
			pMap[i] = 0;
			}
		else {
			++nGroup;
			}
		}
	/*
	 ** Remap the groups.
	 */
	for (pi=0;pi<kd->nParticles;++pi) {
		kd->piGroup[pi] = pMap[kd->piGroup[pi]];
		}
	free(pMap);
	free(pnMembers);
	kd->nGroup = nGroup;
	}


void kdGroupOrder(KD kd)
{
	PINIT *p,t;
	int i,j,iGroup,pi;
	
	/*
	 ** First split off the "non-group" particles.
	 */	
	kd->kdGroup = (KDN *)malloc(kd->nGroup*sizeof(KDN));
	assert(kd->kdGroup != NULL);
	p = kd->pInit;
	pi = 0;
	for (iGroup=0;iGroup<kd->nGroup;++iGroup) {
		i = pi;
		j = kd->nParticles-1;
		while (1) {
			while (kd->piGroup[p[i].iOrder] == iGroup)
				if (++i > j) goto done;
			while (kd->piGroup[p[j].iOrder] != iGroup)
				if (i > --j) goto done;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			}
	done:
		kd->kdGroup[iGroup].pLower = pi;
		kd->kdGroup[iGroup].pUpper = i-1;
		pi = i;
		}
	}


void kdUnbind(KD kd,int iSoftType,float fScoop,int bGasAndDark)
{
	PINIT *q,t;
	KDN *pkdn;
	int iGroup,n,i,j,iBig,pi;
	float hx,hy,hz,dx,dy,dz,dv,dv2,fShift,fCosmo,fTot,fTotBig;
	double dMass,rcm[3],vcm[3],*pdPot,dPot;
	int nUnbind = 0;

	hx = 0.5*kd->fPeriod[0];
	hy = 0.5*kd->fPeriod[1];
	hz = 0.5*kd->fPeriod[2];
	fShift = 1.0/(1.0+kd->z);
	fCosmo = kd->H0*sqrt(1.0+kd->Omega0*kd->z);
	/*
	 ** Now build a tree for the non-grouped particles!
	 ** They are the first "group" in the pInit array after the group order.
	 ** The tree is for scooping particles!
	 */
	kd->nInitActive = kd->kdGroup[0].pUpper+1;
	kdBuildTree(kd);
	printf("Groups before Unbind:%d\n",kd->nGroup-1);
	fflush(stdout);
	for (iGroup=1;iGroup<kd->nGroup;++iGroup) {
		pkdn = &kd->kdGroup[iGroup];
		n = pkdn->pUpper-pkdn->pLower+1;
		/*
		 ** Make temporary particles to hold new positions
		 */
		q = malloc(n*sizeof(PINIT));
		assert(q != NULL);
		/*
		 ** Make all group particles have coordinates relative to 
		 ** a group reference point, given by pGroup[i].rel.
		 */
		for (i=0,pi=pkdn->pLower;i<n;++i,++pi) {
			dx = kd->pInit[pi].r[0] - kd->pGroup[iGroup].rel[0];
			dy = kd->pInit[pi].r[1] - kd->pGroup[iGroup].rel[1];
			dz = kd->pInit[pi].r[2] - kd->pGroup[iGroup].rel[2];
			if (dx > hx) dx -= 2*hx;
			if (dx <= -hx) dx += 2*hx;
			if (dy > hy) dy -= 2*hy;
			if (dy <= -hy) dy += 2*hy;
			if (dz > hz) dz -= 2*hz;
			if (dz <= -hz) dz += 2*hz;
			q[i] = kd->pInit[pi];
			q[i].r[0] = dx;
			q[i].r[1] = dy;
			q[i].r[2] = dz;
			}
		/*
		 ** Calculate center of mass and center of mass 
		 ** velocity for the group.
		 */
		dMass = 0.0;
		for (j=0;j<3;++j) {
			rcm[j] = 0.0;
			vcm[j] = 0.0;
			}
		for (i=0;i<n;++i) {
			dMass += q[i].fMass;
			for (j=0;j<3;++j) {
				rcm[j] += q[i].fMass*q[i].r[j];
				vcm[j] += q[i].fMass*q[i].v[j];
				}
			}
		for (j=0;j<3;++j) {
			rcm[j] /= dMass;
			vcm[j] /= dMass;
			}
		pdPot = (double *)malloc(n*sizeof(double));
		assert(pdPot != NULL);
		kdCellPot(kd,q,n,iSoftType,pdPot);
		if ((kd->inType == DARK|GAS && !bGasAndDark) || 
			kd->inType == DARK|GAS|STAR || kd->inType == DARK|STAR) {
			kdAddScoopPot(kd,q,n,kd->pGroup[iGroup].rCenter,fScoop,
						  kd->pGroup[iGroup].rel,iSoftType,pdPot);
			}
		while (1) {
			/*
			 ** Find the largest total Energy.
			 */
			iBig = 0;
			fTotBig = -1.0;
			for (i=0;i<n;++i) {
				dv2 = 0.0;
				for (j=0;j<3;++j) {
					dv = fShift*(q[i].v[j]-vcm[j]) + fCosmo*(q[i].r[j]-rcm[j]);
					dv2 += dv*dv;
					}
				fTot = 0.5*dv2 - pdPot[i]*(1.0+kd->z);
				if (fTot > fTotBig) {
					fTotBig = fTot;
					iBig = i;
					}
				}
			if (fTotBig < 0) break;
			/*
			 ** Unbind particle iBig!
			 */
			kd->piGroup[q[iBig].iOrder] = 0;
			/*
			 ** Adjust rcm and vcm.
			 */
			dMass -= q[iBig].fMass;
			for (j=0;j<3;++j) {
				rcm[j] += q[iBig].fMass/dMass*(rcm[j]-q[iBig].r[j]);
				vcm[j] += q[iBig].fMass/dMass*(vcm[j]-q[iBig].v[j]);
				}
			++nUnbind;
			--n;
			--kd->kdGroup[iGroup].pUpper;
			/*
			 ** Swap particle iBig and last, swap also corresp. potential E.
			 */
			t = q[iBig];
			q[iBig] = q[n];
			q[n] = t;
			dPot = pdPot[iBig];
			pdPot[iBig] = pdPot[n];
			pdPot[n] = dPot;
			if (kd->inType == DARK || kd->inType == STAR) {
				/*
				 ** Modify the potential energies.
				 */
				kdSubPot(kd,q,n,&q[n],iSoftType,pdPot);
				}
			}
		free(pdPot);		
		free(q);
		kd->pGroup[iGroup].fMass = dMass;
		for (j=0;j<3;++j) {
			kd->pGroup[iGroup].vcm[j] = vcm[j];
			}
		}
	printf("Number of particles Unbound:%d\n",nUnbind);
	fflush(stdout);
	}


int CmpInit(const void *p1,const void *p2)
{
	PINIT *a = (PINIT *)p1;
	PINIT *b = (PINIT *)p2;

	return(a->iOrder - b->iOrder);
	}


int CmpMove(const void *p1,const void *p2)			
{
	PMOVE *a = (PMOVE *)p1;
	PMOVE *b = (PMOVE *)p2;

	return(a->iOrder - b->iOrder);
	}


void Order(KD kd)
{
	qsort(kd->pInit,kd->nParticles,sizeof(PINIT),CmpInit);
	qsort(kd->pMove,kd->nMove,sizeof(PMOVE),CmpMove);
	}


void kdOutGroup(KD kd,char *pszFile)
{	
	FILE *fp;
	int i;

	printf("Number of Groups:%d\n",kd->nGroup-1);
	fflush(stdout);
	if (pszFile) {
		fp = fopen(pszFile,"w");
		assert(fp != NULL);
		}
	else {
		fp = stdout;
		}
	fprintf(fp,"%d\n",kd->nParticles);
	for (i=0;i<kd->nParticles;++i) fprintf(fp,"%d\n",kd->piGroup[i]);
	fclose(fp);
	}


void kdOutDensity(KD kd,char *pszFile)
{	
	FILE *fp;
	int i;

	if (pszFile) {
		fp = fopen(pszFile,"w");
		assert(fp != NULL);
		}
	else {
		fp = stdout;
		}
	/*
	 ** Make sure the particles are ordered before outputing the densities.
	 */
	Order(kd);
	fprintf(fp,"%d\n",kd->nParticles);
	for (i=0;i<kd->nParticles;++i) {
		fprintf(fp,"%.10g\n",kd->pInit[i].fDensity);
		}
	fclose(fp);
	}


void kdOutVector(KD kd,char *pszFile)
{
	PINIT *p;
	PMOVE *pm;
	FILE *fp;
	int pi,pmi;
	float hx,hy,hz,dx,dy,dz;

	/*
	 ** Make sure the particles are ordered before outputing the vectors.
	 */
	Order(kd);

	p = kd->pInit;
	pm = kd->pMove;
	hx = 0.5*kd->fPeriod[0];
	hy = 0.5*kd->fPeriod[1];
	hz = 0.5*kd->fPeriod[2];
	fp = fopen(pszFile,"w");
	assert(fp != NULL);
	fprintf(fp,"%d\n",kd->nParticles);
	for (pmi=0,pi=0;pmi<kd->nMove;++pmi,++pi) {
		while (p[pi].iOrder < pm[pmi].iOrder) {
			fprintf(fp,"0\n");
			++pi;
			}
		assert(p[pi].iOrder == pm[pmi].iOrder);
		dx = pm[pmi].r[0] - p[pi].r[0];
		if (dx > hx) dx -= 2*hx;
		if (dx <= -hx) dx += 2*hx;
		fprintf(fp,"%g\n",dx);
		}
	for (;pi<kd->nParticles;++pi) fprintf(fp,"0\n");
	for (pmi=0,pi=0;pmi<kd->nMove;++pmi,++pi) {
		while (p[pi].iOrder < pm[pmi].iOrder) {
			fprintf(fp,"0\n");
			++pi;
			}
		assert(p[pi].iOrder == pm[pmi].iOrder);
		dy = pm[pmi].r[1] - p[pi].r[1];
		if (dy > hy) dy -= 2*hy;
		if (dy <= -hy) dy += 2*hy;
		fprintf(fp,"%g\n",dy);
		}
	for (;pi<kd->nParticles;++pi) fprintf(fp,"0\n");
	for (pmi=0,pi=0;pmi<kd->nMove;++pmi,++pi) {
		while (p[pi].iOrder < pm[pmi].iOrder) {
			fprintf(fp,"0\n");
			++pi;
			}
		assert(p[pi].iOrder == pm[pmi].iOrder);
		dz = pm[pmi].r[2] - p[pi].r[2];
		if (dz > hz) dz -= 2*hz;
		if (dz <= -hz) dz += 2*hz;
		fprintf(fp,"%g\n",dz);
		}
	for (;pi<kd->nParticles;++pi) fprintf(fp,"0\n");
	fclose(fp);
	}


void kdWriteGroup(KD kd,char *pszFile)
{
	PINIT *p;
	int iGroup;
	float del,fRad;
	FILE *fp;
	int i,j;
	struct dump h;
	struct star_particle sp;

	/*
	 ** Calculate the radius of each group.
	 */
	for (i=1;i<kd->nGroup;++i) {
		kd->pGroup[i].fRadius = 0.0;
		}
	p = kd->pInit;
	for (i=0;i<kd->nParticles;++i) {
		iGroup = kd->piGroup[p[i].iOrder];
		if (!iGroup) continue;
		fRad = 0.0;
		for (j=0;j<3;++j) {
			del = p[i].r[j] - kd->pGroup[iGroup].rCenter[j];
			if (del > 0.5*kd->fPeriod[j]) del -= kd->fPeriod[j];
			if (del <= -0.5*kd->fPeriod[j]) del += kd->fPeriod[j];
			fRad += del*del;
			}
		fRad = sqrt(fRad);
		if (fRad > kd->pGroup[iGroup].fRadius) 
			kd->pGroup[iGroup].fRadius = fRad;
		}
	/*
	 ** Now output the GTP file.
	 */
	fp = fopen(pszFile,"wb");
	assert(fp != NULL);
	h.nbodies = kd->nGroup-1;
	h.nstar = kd->nGroup-1;
	h.ndark = 0;
	h.nsph = 0;
	h.ndim = 3;
	h.time = kd->fTime;
	fwrite(&h,sizeof(struct dump),1,fp);
	for (i=1;i<kd->nGroup;++i) {
		sp.mass = kd->pGroup[i].fMass;
		for (j=0;j<3;++j) {
			sp.pos[j] = kd->pGroup[i].rCenter[j];
			sp.vel[j] = kd->pGroup[i].vcm[j];
			}
		sp.eps = kd->pGroup[i].fRadius;
		sp.metals = 0.0;
		sp.tform = kd->fTime;
		sp.phi = 0.0;
		fwrite(&sp,sizeof(struct star_particle),1,fp);
		}
	fclose(fp);
	}


int CmpRadius(const void *p1,const void *p2)
{
	PINIT *a = (PINIT *)p1;
	PINIT *b = (PINIT *)p2;

	if(a->fBall2 > b->fBall2)
	    return 1;
	if(a->fBall2 < b->fBall2)
	    return -1;
	return 0;
	}


void kdOutStats(KD kd,char *pszFile, float fDensMin, float fTempMax)
{
	FILE *fp;
	int iGroup,i,pi,j,k,n;
	PINIT *q;
	KDN *pkdn;
	float hx,hy,hz,dx,dy,dz;
	float fTotMass, fGasMass, fStarMass, fVcirc;
	float flVcirc;
	float fmVcirc;

	fp = fopen(pszFile,"w");
	assert(fp != NULL);
	hx = 0.5*kd->fPeriod[0];
	hy = 0.5*kd->fPeriod[1];
	hz = 0.5*kd->fPeriod[2];
	/*
	 * Particles should be in Group Order.
	 */
	for (iGroup=1;iGroup<kd->nGroup;++iGroup) {
		pkdn = &kd->kdGroup[iGroup];
		n = pkdn->pUpper-pkdn->pLower+1;
		/*
		 ** Make temporary particles to hold new positions
		 */
		q = malloc(n*sizeof(PINIT));
		assert(q != NULL);
		/*
		 ** Make all group particles have coordinates relative to 
		 ** a group reference point, given by pGroup[i].rel.
		 */
		for (i=0,pi=pkdn->pLower;i<n;++i,++pi) {
			dx = kd->pInit[pi].r[0] - kd->pGroup[iGroup].rCenter[0];
			dy = kd->pInit[pi].r[1] - kd->pGroup[iGroup].rCenter[1];
			dz = kd->pInit[pi].r[2] - kd->pGroup[iGroup].rCenter[2];
			if (dx > hx) dx -= 2*hx;
			if (dx <= -hx) dx += 2*hx;
			if (dy > hy) dy -= 2*hy;
			if (dy <= -hy) dy += 2*hy;
			if (dz > hz) dz -= 2*hz;
			if (dz <= -hz) dz += 2*hz;
			q[i] = kd->pInit[pi];
			q[i].r[0] = dx;
			q[i].r[1] = dy;
			q[i].r[2] = dz;
			}
		/*
		 * Sort particles by distance from center.
		 */
	    for(j = 0; j < n; ++j) {
			float radius2 = 0.0;
			for(k = 0; k < 3; ++k) {
				radius2 += q[j].r[k]*q[j].r[k];
				}
			q[j].fBall2 = radius2;
			}
	    qsort(q,n,sizeof(PINIT),CmpRadius);
	    fGasMass = fStarMass = fTotMass = fVcirc = 0.0;
	    fmVcirc = 0.0;
	    for(j = 0; j < n; j++) {
			fTotMass += q[j].fMass;
			if(kd->G*fTotMass/sqrt(q[j].fBall2) > fVcirc)
				fVcirc = kd->G*fTotMass/sqrt(q[j].fBall2);
			if (kdParticleType(kd,q[j].iOrder) == GAS &&
				q[j].fDensity >= fDensMin && q[j].fTemp <= fTempMax)
				fGasMass += q[j].fMass;
			if (kdParticleType(kd,q[j].iOrder) == STAR)
				fStarMass += q[j].fMass;
			if(j == n/2)
				fmVcirc = kd->G*fTotMass/sqrt(q[j].fBall2); 
			}
	    flVcirc = sqrt(kd->G*fTotMass/sqrt(q[n-1].fBall2));
	    fprintf(fp, "%d %d %g %g %g %g %g %g %g %g %g %g %g %g %g\n", i, n,
				fTotMass, fGasMass, fStarMass, sqrt(fVcirc),
				fmVcirc, flVcirc, sqrt(q[n-1].fBall2),
				kd->pGroup[i].rCenter[0],
				kd->pGroup[i].rCenter[1],
				kd->pGroup[i].rCenter[2], 
				kd->pGroup[i].vcm[0],
				kd->pGroup[i].vcm[1],
				kd->pGroup[i].vcm[2]);
		}
	fclose(fp);
	}


void kdFinish(KD kd)
{
	if (kd->pMove) free(kd->pMove);
	if (kd->pInit) free(kd->pInit);
	if (kd->pGroup) free(kd->pGroup);
	if (kd->kdNodes) free(kd->kdNodes);
	if (kd->kdGroup) free(kd->kdGroup);
	if (kd->piGroup) free(kd->piGroup);
	free(kd);
	}












