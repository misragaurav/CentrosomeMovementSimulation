#include "proto.h"
#define boundary_depoly

void polymerize(struct mt *mt, struct params *params, real dt)
{
  int n, k, i=0;
  int G;

  G=0;
  for(n=0; n<Nmt; n++)
    G += mt[n].N;

for(n=0; n<Nmt; n++)
{
  // switching between poly and depoly
  if (mt[n].status == 0 && gsl_rng_uniform(rng) < params->k_catast*dt)
    mt[n].status = 1;
  if (mt[n].status == 1 && gsl_rng_uniform(rng) < params->k_recov*dt)
    mt[n].status = 0;

  // growing and shrinking
  if (mt[n].status == 0) // growth mode
    {
#ifdef stallforceON
      if(mt[n].F3tip > -params->F3stall)   // polymerization dependent on stall force
#endif
	if(G<Gmax)
	  {mt[n].dL = mt[n].dL + params->Vpoly*dt;}

      if (mt[n].dL > ds)
	{
	  for(i=0; i<(int)(mt[n].dL/ds); i++)
	    {
	      if (mt[n].N >= Nmax)
		{
		  mt[n].dL = 0;  // reset dL if N=Nmax and add no more
		}
	      else
		{
		  addseg(&mt[n], params);
		}
	    }
	}
    }

  else  //shrink
    {
      mt[n].dL = mt[n].dL - params->Vdepoly*dt;
      if (mt[n].dL <0)
	{
	  for(i=0; i>(int)(mt[n].dL/ds); i--)
	    {
	      if (mt[n].N <= 1)
		{
		  randomorient(&mt[n]);
		  mt[n].dL = 0;   // reset dL if N=1
		}
	      else
		{
		  delseg(&mt[n], params);
		}
	    }
	}
    }
}
}

void addseg(struct mt *mt, struct params *params)
{
  struct seg *ps1, *ps2;
  real rx, ry, rz, drx, dry, drz, ds1, ds2, ds3;
  real cosx, cosy, cosz;
  int proximity[3], normal[3];
  int boundaryflag, m;
  real prob[3], randnum;

  ds1=0;
  ds2=0;
  ds3=ds;

  ps1 = &mt->seg[mt->N-1];
  rotmatrix(ps1);

  drx = ps1->R[0][0]*ds1 + ps1->R[1][0]*ds2 + ps1->R[2][0]*ds3;
  dry = ps1->R[0][1]*ds1 + ps1->R[1][1]*ds2 + ps1->R[2][1]*ds3;
  drz = ps1->R[0][2]*ds1 + ps1->R[1][2]*ds2 + ps1->R[2][2]*ds3;

  rx = ps1->rx + drx;
  ry = ps1->ry + dry;
  rz = ps1->rz + drz;

  mt->dL = mt->dL - ds;
  mt->N  = mt->N + 1;
  //      printf("Nmax: %i \t N: %i \n",Nmax, mt->N);
  ps2     = &mt->seg[mt->N-1];
  ps2->rx = rx;
  ps2->ry = ry;
  ps2->rz = rz;
  ps2->q0 = ps1->q0;
  ps2->qx = ps1->qx;
  ps2->qy = ps1->qy;
  ps2->qz = ps1->qz;

  /*
  //calculate if there is extra friction at the tip
  for(m=0; m<mt->N; m++)
    {
      ps2 = &mt->seg[m];

      ps2->zt1 = params->fric_perp;
      ps2->zt2 = params->fric_perp;
      ps2->zt3 = params->fric_para;

      if(m>mt->N-4)
	{
	  checkboundarysquare(&proximity[0], ps2->rx, ps2->ry, ps2->rz, Dia-2*cutoffdist, thickness-2*cutoffdist);

	  if(proximity[0] != 0 || proximity[1] != 0) // close to the x-z or y-z planes
	    {
	      ps2->zt1 = frictionattip*params->fric_perp;
	      ps2->zt2 = frictionattip*params->fric_perp;
	    }
	}
    }
  */

#ifdef pinning
  // Square boundary
  rotmatrix(ps2);
  checkboundarysquare(&proximity[0], ps2->rx, ps2->ry, ps2->rz, Dia-2*cutoffdist/1000.0, thickness-2*cutoffdist/1000.0);

  prob[0] = cosx = proximity[0]*ps2->R[2][0]; // cosine of the angle between normal to x and d3
  prob[1] = cosy = proximity[1]*ps2->R[2][1]; // cosine of the angle between normal to y and d3
  prob[2] = cosz = proximity[2]*ps2->R[2][2]; // cosine of the angle between normal to z and d3

  randnum = gsl_rng_uniform(rng);

  if(randnum < prob[0] || randnum < prob[1] || randnum < prob[2])
    {
      mt->pinned = 1; //pinned
    }
  else
    {
      mt->pinned = 0; //free
    }

#endif

  /*
  // Circular boundary
  checkboundary(ps2->rx, ps2->ry, ps2->rz, Dia-2*cutoffdist, thickness-2*cutoffdist);

  if(boundaryflag != 0)
    if(gsl_rng_uniform(rng) < pinnedMTfrac)
      mt->pinned = 1; //pinned
    else
      mt->pinned = 0; //free
  else
    mt->pinned = 0;
  */

}

void delseg(struct mt *mt, struct params *params)
{
  struct seg *ps2;
  real cosx, cosy, cosz;
  int proximity[3], normal[3];
  int boundaryflag, m;
  real prob[3], randnum;

  ps2   = &mt->seg[mt->N-1];

  ps2->rx = 0;
  ps2->ry = 0;
  ps2->rz = 0;

  ps2->q0 = 0;
  ps2->qx = 0;
  ps2->qy = 0;
  ps2->qz = 0;
  rotmatrix(ps2);

  mt->N =  mt->N - 1;
  mt->dL = mt->dL + ds;

  /*

  //calculate if there is extra friction at the tip
  for(m=0; m<mt->N; m++)
    {
      ps2 = &mt->seg[m];

      ps2->zt1 = params->fric_perp;
      ps2->zt2 = params->fric_perp;
      ps2->zt3 = params->fric_para;

      if(m>mt->N-4)
	{
	  checkboundarysquare(&proximity[0], ps2->rx, ps2->ry, ps2->rz, Dia-2*cutoffdist, thickness-2*cutoffdist);

	  if(proximity[0] != 0 || proximity[1] != 0) // close to the x-z or y-z planes
	    {
	      ps2->zt1 = frictionattip*params->fric_perp;
	      ps2->zt2 = frictionattip*params->fric_perp;
	    }
	}
    }
  */

#ifdef pinning
  // Square boundary
  rotmatrix(ps2);
  checkboundarysquare(&proximity[0], ps2->rx, ps2->ry, ps2->rz, Dia-2*cutoffdist/1000.0, thickness-2*cutoffdist/1000.0);

  prob[0] = cosx = proximity[0]*ps2->R[2][0]; // cosine of the angle between normal to x and d3
  prob[1] = cosy = proximity[1]*ps2->R[2][1]; // cosine of the angle between normal to y and d3
  prob[2] = cosz = proximity[2]*ps2->R[2][2]; // cosine of the angle between normal to z and d3

  randnum = gsl_rng_uniform(rng);

  if(randnum < prob[0] || randnum < prob[1] || randnum < prob[2])
    {
      mt->pinned = 1; //pinned
    }
  else
    {
      mt->pinned = 0; //free
    }

#endif


  /*
  // Circular boundary
  boundaryflag = checkboundary(ps2->rx, ps2->ry, ps2->rz, Dia-2*cutoffdist, thickness-2*cutoffdist);

  if(boundaryflag != 0)
    if(gsl_rng_uniform(rng) < pinnedMTfrac)
      mt->pinned = 1; //pinned
    else
      mt->pinned = 0; //free
  else
    mt->pinned = 0;
  */

}
