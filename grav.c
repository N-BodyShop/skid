#include <stdio.h>
#include <math.h>
#include "kd.h"
#include "grav.h"



void kdCellPot(KD kd,PINIT *p,int n,int iSoftType,double *pdPot)
{
	int i,j;
	float x,y,z,dx,dy,dz,d2,twoh,dir;

	for (i=0;i<n;++i) pdPot[i] = 0.0;
	for (i=0;i<(n-1);++i) {
		x = p[i].r[0];
		y = p[i].r[1];
		z = p[i].r[2];
		for (j=i+1;j<n;++j) {
			dx = x - p[j].r[0];
			dy = y - p[j].r[1];
			dz = z - p[j].r[2];
			d2 = dx*dx + dy*dy + dz*dz;
			twoh = p[i].fSoft + p[j].fSoft;
			switch (iSoftType) {
			case PLUMMER:
				dir = 1.0/sqrt(d2 + 0.25*twoh*twoh);
				break;
			default:
			case SPLINE:
				SPLINE_POT(d2,twoh,dir);
				}
			pdPot[i] += kd->G*p[j].fMass*dir;
			pdPot[j] += kd->G*p[i].fMass*dir;
			}
		}
	}


void kdSubPot(KD kd,PINIT *p,int n,PINIT *ps,int iSoftType,double *pdPot)
{
	int i;
	float dx,dy,dz,d2,twoh,dir;

	for (i=0;i<n;++i) {
		dx = ps->r[0] - p[i].r[0];
		dy = ps->r[1] - p[i].r[1];
		dz = ps->r[2] - p[i].r[2];
		d2 = dx*dx + dy*dy + dz*dz;
		twoh = ps->fSoft + p[i].fSoft;
		switch (iSoftType) {
		case PLUMMER:
			dir = 1.0/sqrt(d2 + 0.25*twoh*twoh);
			break;
		default:
		case SPLINE:
			SPLINE_POT(d2,twoh,dir);
			}
		pdPot[i] -= kd->G*ps->fMass*dir;
		}
	}


void kdAddScoopPot(KD kd,PINIT *pGroup,int n,float *ri,float fScoop,float *rel,
				   int iSoftType,double *pdPot)
{
	KDN *c;
	PINIT *p;
	int pj,cp,i;
	float dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2,fBall2,nx,ny,nz;
	float d2,twoh,dir;

	c = kd->kdNodes;
	if (!c) return;
	p = kd->pInit;
	lx = kd->fPeriod[0];
	ly = kd->fPeriod[1];
	lz = kd->fPeriod[2];
	x = ri[0];
	y = ri[1];
	z = ri[2];
	fBall2 = fScoop*fScoop;
	cp = ROOT;
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
					/*
					 ** Add this particle's potential to all the particles in
					 ** pGroup.
					 */
					nx = p[pj].r[0] - rel[0];
					ny = p[pj].r[1] - rel[1];
					nz = p[pj].r[2] - rel[2];
					if (nx > 0.5*lx) nx -= lx;
					if (nx <= -0.5*lx) nx += lx;
					if (ny > 0.5*ly) ny -= ly;
					if (ny <= -0.5*ly) ny += ly;
					if (nz > 0.5*lz) nz -= lz;
					if (nz <= -0.5*lz) nz += lz;					
					for (i=0;i<n;++i) {
						dx = nx - pGroup[i].r[0];
						dy = ny - pGroup[i].r[1];
						dz = nz - pGroup[i].r[2];
						d2 = dx*dx + dy*dy + dz*dz;
						twoh = p[pj].fSoft + pGroup[i].fSoft;
						switch (iSoftType) {
						case PLUMMER:
							dir = 1.0/sqrt(d2 + 0.25*twoh*twoh);
							break;
						default:
						case SPLINE:
							SPLINE_POT(d2,twoh,dir);
							}
						pdPot[i] += kd->G*p[pj].fMass*dir;
						}
					}
				}
			}
	GetNextCell:
		SETNEXT(cp);
		if (cp == ROOT) break;
		}
	}







