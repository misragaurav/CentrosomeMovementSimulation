#include "proto.h"

#define convertdist  1e3   //micro to nm
#define convertt     1e3   //ms to micro s
#define convertYI    1e6   //pN micron^2 to fN nm^2
#define convertk     1e-3  //spring constant - pN/micron to pN/nm
#define convertfric  1e-3  //pN.ms/micron^2 to N.s/m^2
#define convertforce 1     //pN to pN

int checkboundary(real rx, real ry, real rz, real d, real thick)
{
  int x,z;
  int status;
  
  if(rx*rx + ry*ry > d*d/4.0)
    x=1;
  else
    x=0;
  
  if(rz > thick/2.0 || rz < -thick/2.0)
    z=1;
  else
    z=0;
  
  if(x==0 && z==0) status = 0;
  if(x==1 && z==0) status = 1;
  if(x==0 && z==1) status = 2;
  if(x==1 && z==1) status = 3;
  
  return status;
}

void checkboundarysquare(int *proximity, real rx, real ry, real rz, real side, real thick)
{
  if(rx > side/2.0)        proximity[0]= 1;
  else if(rx < -side/2.0)  proximity[0]=-1;
  else                     proximity[0]= 0;

  if(ry > side/2.0)        proximity[1]= 1;
  else if(ry < -side/2.0)  proximity[1]=-1;
  else                     proximity[1]= 0;
 
  if(rz > thick/2.0)       proximity[2]= 1;
  else if(rz < -thick/2.0) proximity[2]=-1;
  else                     proximity[2]= 0;
}

void initrodparams(struct params *params, struct rodparams *rparams, struct rod *rod, real dt)
{
  int n;
  real ext_mod_correct = 0.5; //extensional modulus is softer by a factor of 2 due to the MT being hollow.

  // time step for position update
  rparams->dt  = convertt*dt;

  //parameters from the main code
  rparams->le  = convertdist*ds;                // microns to nm
  rparams->YI  = convertYI*params->YI;
  rparams->k   = convertk*params->k;           // force per unit length - pN/micron to fN/nm 
  rparams->F0  = convertk*params->F0;          // force per unit length - pN/micron to fN/nm
  
  rparams->zt1 = convertfric*params->fric_perp; // 
  rparams->zt2 = convertfric*params->fric_perp;	// fric/length converted to fric/mass in momenta.c
  rparams->zt3 = convertfric*params->fric_para;					       

  // Independent parameters - in units of rod code
  rparams->dia   = 25;    // nm
  rparams->rho   = 1e-9;  // fg/nm^3
  rparams->sigma = 0.33;  // dimensionless

  rparams->s01 = 0;       // Reference strains
  rparams->s02 = 0;
  rparams->s03 = 1;
  rparams->b01 = 0;
  rparams->b02 = 0;
  rparams->b03 = 0;

  // Derived parameters
  rparams->A = Pi*rparams->dia*rparams->dia/4.0; //Area
  rparams->I = rparams->A*rparams->A/(4.0*Pi);   //Area moment
  rparams->Y = rparams->YI/rparams->I;
  rparams->G = rparams->Y/(2+2*rparams->sigma);

  rparams->xm =   rparams->A*rparams->rho; // mass per unit length
            
  rparams->i1 =   rparams->I*rparams->rho; // Moment of inertia per unit length
  rparams->i2 =   rparams->I*rparams->rho;
  rparams->i3 = 2*rparams->I*rparams->rho;

  rparams->cs1 =  ext_mod_correct*rparams->A*rparams->G;  // Moments of inertia per unit length
  rparams->cs2 =  ext_mod_correct*rparams->A*rparams->G;       
  rparams->cs3 =  ext_mod_correct*rparams->A*rparams->Y;
  rparams->cb1 =   rparams->I*rparams->Y;
  rparams->cb2 =   rparams->I*rparams->Y;
  rparams->cb3 = 2*rparams->I*rparams->G;

  // Cell size
  rparams->d     = convertdist*Dia;
  rparams->thick = convertdist*thickness;

  //Allocate memory to nodes and segments
  for(n=0; n<Nmt; n++)
    {
      rod[n].nodes = (struct node *) calloc(Nmax, sizeof(struct node) ); //Nmax used
      rod[n].elems = (struct elem *) calloc(Nmax, sizeof(struct elem) ); //Nmax used
    }
}

void interfacein(struct mt *mt, struct rod *rod, struct centro *centro)
{
  struct seg *ps;
  struct node *pn;
  int i;

  rod->n_nodes = mt->N;
  rod->n_elems = mt->N+1;

  rod->pinned = mt->pinned;
  rod->dL = convertdist*mt->dL;

  // Populate the quaternions in rod nodes
  for(i=0;i<rod->n_nodes;i++)
    {  
      ps = &mt->seg[i];
      pn = &rod->nodes[i];

      pn->q0 = ps->q0;
      pn->qx = ps->qx;
      pn->qy = ps->qy;
      pn->qz = ps->qz;
      
      pn->rx = convertdist*ps->rx;
      pn->ry = convertdist*ps->ry;
      pn->rz = convertdist*ps->rz;

      pn->zt1 = convertfric*ps->zt1;
      pn->zt2 = convertfric*ps->zt2;
      pn->zt3 = convertfric*ps->zt3;

    }

  rod->centrorx = convertdist*centro->rx;
  rod->centrory = convertdist*centro->ry;
  rod->centrorz = convertdist*centro->rz;

}

void interfaceout(struct mt *mt, struct rod *rod)
{
  struct seg  *ps;
  struct node *pn;
  int i;

  for(i=0;i<rod->n_nodes;i++)
    {  
      pn = &rod->nodes[i];
      ps = &mt->seg[i];

      ps->q0 = pn->q0;
      ps->qx = pn->qx;
      ps->qy = pn->qy;
      ps->qz = pn->qz;

      ps->rx = pn->rx/convertdist;
      ps->ry = pn->ry/convertdist;
      ps->rz = pn->rz/convertdist;
      
      rotmatrix(ps);
    }

  mt->F3tip = convertforce*rod->F3tip;

}
