#include "proto.h"

void Qinv (real Q[4][8])  /* Invert Q matrix */
{
  real ratio;
  int    i, j, k;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      Q[i][j+4] = 0.0;
  for (i = 0; i < 4; i++)
    Q[i][i+4] = 1.0;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      {
        if (i == j)  continue;
        ratio = Q[j][i]/Q[i][i];
        for (k = 0; k < 8; k++)
          Q[j][k] -= ratio*Q[i][k];
      }
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      Q[i][j+4] /= Q[i][i];
}


void initparams(struct params *params, struct mt *mt)
{
  struct seg *ps2;
  int proximity[3], m;

  // Init params
  params->YI          = 25; // in pN micron^2
  params->Vpoly       = 1e-4; //micron per ms
  params->Vdepoly     = 3e-4; //micron per ms

  params->k_catast    = 5e-5; // per ms
  params->k_recov     = 2e-4; // per ms
  params->Vmax        = 1e-3; // micron per ms
  params->Fmax 	      = 10;   // pN
  params->koff 	      = 1e-3; // per ms
  params->kappa       = 1000; // pN/micron

  params->rho_motor   = 5;   // per micron

  params->F0          = params->rho_motor*params->Fmax/( params->Fmax*params->koff/(params->kappa*params->Vmax) + 1 );
  params->fric_para   = params->F0/params->Vmax;

#ifdef motorforceON
  params->fric_perp   =      params->rho_motor*params->kappa/params->koff;
#else
  params->fric_perp   = 0.01*params->rho_motor*params->kappa/params->koff;
#endif

  //printf("F0: %f \t fric_para: %f \t fric_perp: %f \n", params->F0, params->fric_para, params->fric_perp );

  params->fric_centro = params->fric_para*10;

  params->k           = 5e4; // pN/micron - centrosome MT spring

  params->F3stall     = 100; // MT polymerization stall force in pN

  //calculate if there is extra friction at the tip
  for(m=0; m<mt->N; m++)
    {
      ps2 = &mt->seg[m];

      ps2->zt1 = params->fric_perp;
      ps2->zt2 = params->fric_perp;
      ps2->zt3 = params->fric_para;

      if(m > mt->N-4)
	{
	  checkboundarysquare(&proximity[0], ps2->rx, ps2->ry, ps2->rz, Dia-2*cutoffdist, thickness-2*cutoffdist);

	  if(proximity[0] != 0 || proximity[1] != 0) // close to the x-z or y-z planes
	    {
	      ps2->zt1 = frictionattip*params->fric_perp;
	      ps2->zt2 = frictionattip*params->fric_perp;
	    }
	}
    }
}


void buildmt(struct mt *mt)
{
  struct seg *ps;

  real Q[4][8], Qplus[4][4];
  real rx, ry, rz;
  real q0, qx, qy, qz, q;
  real dq0, dqx, dqy, dqz;
  real b1 = 0.0, b2 = 0.0, b3 = 0;
  real s1 = 0, s2 = 0, s3 = 1;
  int k, l, m;

  ps = &mt->seg[0];
  
  q0 = ps->q0;
  qx = ps->qx;
  qy = ps->qy;
  qz = ps->qz;

  for (m = 1; m < mt->N; m++)
    {
      ps = &mt->seg[m];

      /*
      if((mt->id)%2==0)
	{	
	  b2 = 0.5;
	}
      else
	{
	  if(m<=mt->N/4)
	    b2=0;
	  if(m>mt->N/4 && m<=3*mt->N/4)
	    b2= 0.04*cos(2*Pi*(m-mt->N/4)/(2*mt->N/4));
	  if(m>3*mt->N/4)
	    b2=0;
	}
      */

      Q[0][0] =  0.0; Q[0][1] = -b1;  Q[0][2] = -b2;  Q[0][3] = -b3;
      Q[1][0] =  b1;  Q[1][1] = 0.0;  Q[1][2] =  b3;  Q[1][3] = -b2;
      Q[2][0] =  b2;  Q[2][1] = -b3;  Q[2][2] = 0.0;  Q[2][3] = -b1;
      Q[3][0] =  b3;  Q[3][1] =  b2;  Q[3][2] = -b1;  Q[3][3] = 0.0;

      for (k = 0; k < 4; k++)
	for (l = 0; l < 4; l++)
	  {
	    Q[k][l] *= 0.25*ds;
	    Qplus[k][l] = Q[k][l];
	  }

      for (k = 0; k < 4; k++)
	for (l = 0; l < 4; l++)
	  Q[k][l] *= -1.0;
      for (k = 0; k < 4; k++)
	{
	  Q[k][k]     += 1.0;
	  Qplus[k][k] += 1.0;
	}
      Qinv (Q);

      dq0 = Qplus[0][0]*q0 + Qplus[0][1]*qx + Qplus[0][2]*qy + Qplus[0][3]*qz;
      dqx = Qplus[1][0]*q0 + Qplus[1][1]*qx + Qplus[1][2]*qy + Qplus[1][3]*qz;
      dqy = Qplus[2][0]*q0 + Qplus[2][1]*qx + Qplus[2][2]*qy + Qplus[2][3]*qz;
      dqz = Qplus[3][0]*q0 + Qplus[3][1]*qx + Qplus[3][2]*qy + Qplus[3][3]*qz;	  

      q0 =  Q[0][4]*dq0 + Q[0][5]*dqx + Q[0][6]*dqy + Q[0][7]*dqz;
      qx =  Q[1][4]*dq0 + Q[1][5]*dqx + Q[1][6]*dqy + Q[1][7]*dqz;
      qy =  Q[2][4]*dq0 + Q[2][5]*dqx + Q[2][6]*dqy + Q[2][7]*dqz;
      qz =  Q[3][4]*dq0 + Q[3][5]*dqx + Q[3][6]*dqy + Q[3][7]*dqz;

      q  = sqrt( q0*q0 + qx*qx + qy*qy + qz*qz);
      q0 /= q; qx /= q; qy /= q;  qz /= q; 
      ps->q0 = q0; ps->qx = qx; ps->qy = qy; ps->qz = qz;

      rotmatrix(ps);
    }

  ps = &mt->seg[0];

  rx = ps->rx;  
  ry = ps->ry;  
  rz = ps->rz;
  
  for (m = 1; m < mt->N; m++)
    {
      ps = &mt->seg[m];

      rx += ds*(ps->R[0][0]*s1 + ps->R[1][0]*s2 + ps->R[2][0]*s3);
      ry += ds*(ps->R[0][1]*s1 + ps->R[1][1]*s2 + ps->R[2][1]*s3);
      rz += ds*(ps->R[0][2]*s1 + ps->R[1][2]*s2 + ps->R[2][2]*s3);

      ps->rx = rx;  ps->ry = ry;  ps->rz = rz;
    }
}

void initmtcentro(struct mt *mt, struct centro *centro)
{
  struct seg *ps;
  int n = 0;

  for(n=0; n<Nmt; n++)
    {
      mt[n].seg = (struct seg *) calloc(Nmax, sizeof(struct seg) ); //Nmax used

      mt[n].id     = n;
      mt[n].status = 0; // Initially growing status for all MTs
      mt[n].dL     = 0;

      mt[n].N      = 2;

      //      mt[n].N      = (int)(1.5*sqrt((Dia/4.0)*(Dia/4.0) + (thickness/2.0)*(thickness/2.0))/ds); // + (int)(gsl_rng_uniform(rng)*thickness/(4*ds)); // start with 1 at least

      //      randomorientcone(&mt[n]); //orient node 0 - populates rotmatrix for node 0

      randomorient(&mt[n]); //orient node 0 - populates rotmatrix for node 0

      // Init Centro 
      centro->rx = 0; centro->ry = 0; centro->rz = 0;

      ps = &mt[n].seg[0];  // placing all MTs at centrosome
      ps->rx = centro->rx;
      ps->ry = centro->ry;
      ps->rz = centro->rz;      

      buildmt(&mt[n]);
    }
}


void popmtcentro(char *checkpointfile, struct mt *mt, struct centro *centro, int *t)
{
  struct seg *ps;  
  FILE  *file_ptr=0;
  int n=0, k=0;

  for(n=0; n<Nmt; n++) mt[n].seg = (struct seg *) calloc(Nmax, sizeof(struct seg) ); //Nmax used

  file_ptr = fopen(checkpointfile, "rb");

  if(file_ptr==0) sprintf("file %s missing \n", checkpointfile);

  gsl_rng_fread(file_ptr, rng);

  fread(t, sizeof(int), 1, file_ptr);

  fread(&centro->rx, sizeof(real), 1, file_ptr);
  fread(&centro->ry, sizeof(real), 1, file_ptr);
  fread(&centro->rz, sizeof(real), 1, file_ptr);

  for(n=0; n<Nmt; n++)
    {
      fread(&mt[n].id,      sizeof(int),  1, file_ptr);
      fread(&mt[n].N,       sizeof(int),  1, file_ptr);
      fread(&mt[n].status,  sizeof(int),  1, file_ptr);
      fread(&mt[n].pinned,  sizeof(int),  1, file_ptr);
      fread(&mt[n].dL,      sizeof(real), 1, file_ptr);
    }

  for(n=0; n<Nmt; n++)
    {
      for (k = 0; k < mt[n].N; k++)
	{
	  ps = &(mt[n].seg[k]);

	  fread(&ps->rx, sizeof(real), 1, file_ptr);
	  fread(&ps->ry, sizeof(real), 1, file_ptr);
	  fread(&ps->rz, sizeof(real), 1, file_ptr);
	  fread(&ps->q0, sizeof(real), 1, file_ptr);
	  fread(&ps->qx, sizeof(real), 1, file_ptr);
	  fread(&ps->qy, sizeof(real), 1, file_ptr);
	  fread(&ps->qz, sizeof(real), 1, file_ptr);
	}
    } 
  fclose(file_ptr);
}
