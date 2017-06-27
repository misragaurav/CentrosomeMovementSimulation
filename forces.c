#include "rod_proto.h"

void motorforce(struct rod *rod, struct rodparams *rparams)
{
  struct node *pn;
  int m;
  
  real f1 = 0;
  real f2 = 0;
  real f3 = rparams->F0; // motor force per unit length
  real t1 = 0; 
  real t2 = 0; 
  real t3 = 0; 

  for (m=0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];

      pn->f1 += f1;	//force per unit length
      pn->f2 += f2;
      pn->f3 += f3;
      pn->t1 += t1;
      pn->t2 += t2;
      pn->t3 += t3;
    }
}

void centroforce(struct rod *rod, struct rodparams *rparams)
{
  struct node *pn;
  real Fx, Fy, Fz;
  
  pn = &rod->nodes[0];

  Fx = -rparams->k*(pn->rx - rod->centrorx)/rparams->le; //force per unit length
  Fy = -rparams->k*(pn->ry - rod->centrory)/rparams->le;
  Fz = -rparams->k*(pn->rz - rod->centrorz)/rparams->le;

  pn->f1 += pn->d[0][0]*Fx + pn->d[0][1]*Fy + pn->d[0][2]*Fz;
  pn->f2 += pn->d[1][0]*Fx + pn->d[1][1]*Fy + pn->d[1][2]*Fz;
  pn->f3 += pn->d[2][0]*Fx + pn->d[2][1]*Fy + pn->d[2][2]*Fz;
  
}

void exvolcircle(struct rod *rod, struct rodparams *rparams)
{
  struct node *pn;
  real Fx, Fy, Fz;
  real F1, F2, F3;
  real dist;
  int boundaryflag, m;
  
  for (m = 0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];
      
      boundaryflag = checkboundary(pn->rx, pn->ry, pn->rz, rparams->d-2*cutoffdist, rparams->thick-2*cutoffdist);

      dist = sqrt(pn->rx*pn->rx + pn->ry*pn->ry);
      
      if(boundaryflag==0)
	{
	  Fx = 0;
	  Fy = 0;
	  Fz = 0;
	}
      if(boundaryflag==1)
	{
	  Fx = -U0*2*pn->rx*exp(-(rparams->d/2.0 - dist)/lambda)/(lambda*dist);
	  Fy = -U0*2*pn->ry*exp(-(rparams->d/2.0 - dist)/lambda)/(lambda*dist);
	  Fz =   0;
	}
      if(boundaryflag==2)
	{
	  Fx = 0;
	  Fy = 0;
	  Fz = -U0*pn->rz*exp(-(rparams->thick/2.0 - fabs(pn->rz))/lambda)/(fabs(pn->rz)*lambda);	  	  
	}
      if(boundaryflag==3)
	{
	  Fx = -U0*2*pn->rx*exp(-(rparams->d/2.0 - dist)/lambda)/(lambda*dist);
	  Fy = -U0*2*pn->ry*exp(-(rparams->d/2.0 - dist)/lambda)/(lambda*dist);
	  Fz = -U0*pn->rz*exp(-(rparams->thick/2.0 - fabs(pn->rz))/lambda)/(fabs(pn->rz)*lambda);	  	  
	}

      F1 = pn->d[0][0]*Fx + pn->d[0][1]*Fy + pn->d[0][2]*Fz;
      F2 = pn->d[1][0]*Fx + pn->d[1][1]*Fy + pn->d[1][2]*Fz;
      F3 = pn->d[2][0]*Fx + pn->d[2][1]*Fy + pn->d[2][2]*Fz;
      
      pn->f1 += F1/rparams->le;
      pn->f2 += F2/rparams->le;
      pn->f3 += F3/rparams->le;
      
    }
  
  rod->F3tip = F3; // For comparison with stall force

      pn  = &rod->nodes[rod->n_nodes-1]; 	 
      boundaryflag = checkboundary(pn->rx, pn->ry, pn->rz, rparams->d-2*cutoffdist, rparams->thick-2*cutoffdist);
      if(boundaryflag != 0 && rod->pinned == 1)
	{
	  pn->f1 = 0;			    	 
	  pn->f2 = 0;			    	 
	  /*      
	  pn->f3 = 0;        		    	 
	  pn->t1 = 0;			    	 
	  pn->t2 = 0;			    	 
	  pn->t3 = 0;                        	 
	  */
	}
}

void exvolsquare(struct rod *rod, struct rodparams *rparams)
{
  struct node *pn;
  real drx, dry, drz, rx, ry, rz;
  real Fx, Fy, Fz;
  real F1, F2, F3;
  real dist;
  int proximity[3], m;

  for (m = 0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];

      if(m==rod->n_nodes-1) // for the last node add dL
	{
	  drx = pn->d[2][0]*rod->dL;
	  dry = pn->d[2][1]*rod->dL;
	  drz = pn->d[2][2]*rod->dL;

	  rx = pn->rx + drx;
	  ry = pn->ry + dry;
	  rz = pn->rz + drz;
	}
      else
	{
	  rx = pn->rx;
	  ry = pn->ry;
	  rz = pn->rz;
	}
      
      checkboundarysquare(&proximity[0], rx, ry, rz, rparams->d-2*cutoffdist, rparams->thick-2*cutoffdist);

      Fx = (-U0*exp(-(rparams->d/2.0     - fabs(rx))/lambda)/lambda)*proximity[0];
      Fy = (-U0*exp(-(rparams->d/2.0     - fabs(ry))/lambda)/lambda)*proximity[1];
      Fz = (-U0*exp(-(rparams->thick/2.0 - fabs(rz))/lambda)/lambda)*proximity[2];  	  

      F1 = pn->d[0][0]*Fx + pn->d[0][1]*Fy + pn->d[0][2]*Fz;
      F2 = pn->d[1][0]*Fx + pn->d[1][1]*Fy + pn->d[1][2]*Fz;
      F3 = pn->d[2][0]*Fx + pn->d[2][1]*Fy + pn->d[2][2]*Fz;
      
      pn->f1 += F1/rparams->le;
      pn->f2 += F2/rparams->le;
      pn->f3 += F3/rparams->le;

      /*
      if(m==rod->n_nodes-1) // additional torque on last node
	{
	  pn->t1 +=  -rod->dL*F2/rparams->le;
	  pn->t2 +=   rod->dL*F1/rparams->le;
	  pn->t3 += 0*rod->dL*F3/rparams->le;
	}
      */

    }
  
  rod->F3tip = F3; // For comparison with stall force

#ifdef pinning

  if( rod->pinned == 1 )
    {
      pn  = &rod->nodes[rod->n_nodes-1]; 	 

      pn->f1 = 0;			    	 
      pn->f2 = 0;			    	 
            
      // pn->f3 = 0;        		    	 
      /*
      pn->t1 = 0;			    	 
      pn->t2 = 0;			    	 
      pn->t3 = 0;
      */
      
    }
#endif
}

void force(struct rod *rod, struct rodparams *rparams)  /* Accelerations in BFF */
{
  struct elem *pe;
  struct node *pn, *pn1, *pn2;
  real ep[3][4];
  real drx, dry, drz, dq0, dqx, dqy, dqz;
  real s1, s2, s3, b1, b2, b3;
  real f1, f2, f3, t1, t2, t3;
  real fx, fy, fz, t0, tx, ty, tz;
  real tx0, tx1, tx2, tx3, tdq;
  real le;
  int  m;

  le = rparams->le;

  // Setting nodal forces to zero at the beginning of each time step
  for (m = 0; m < rod->n_nodes; m++)
  {
      pn = &rod->nodes[m];
      pn->fx = pn->fy = pn->fz = 0.0;
      pn->t0 = pn->tx = pn->ty = pn->tz = 0.0;
  }

  /*  Calculate nodal forces and torques  */
  for (m = 1; m < rod->n_elems-1; m++)
    {
      pe  = &rod->elems[m];              
      pn1 = &rod->nodes[m-1];
      pn2 = &rod->nodes[m];

      drx = (pn2->rx - pn1->rx)/le;        /* Derivatives */
      dry = (pn2->ry - pn1->ry)/le;
      drz = (pn2->rz - pn1->rz)/le;
      dq0 = (pn2->q0 - pn1->q0)/le;
      dqx = (pn2->qx - pn1->qx)/le;
      dqy = (pn2->qy - pn1->qy)/le;
      dqz = (pn2->qz - pn1->qz)/le;

      ep[0][0] = -dqx;  ep[0][1] =  dq0; ep[0][2] =  dqz;  ep[0][3] = -dqy;
      ep[1][0] = -dqy;  ep[1][1] = -dqz; ep[1][2] =  dq0;  ep[1][3] =  dqx;
      ep[2][0] = -dqz;  ep[2][1] =  dqy; ep[2][2] = -dqx;  ep[2][3] =  dq0;

      ep[0][0] /= pe->q;  ep[0][1] /= pe->q; ep[0][2] /= pe->q;  ep[0][3] /= pe->q;
      ep[1][0] /= pe->q;  ep[1][1] /= pe->q; ep[1][2] /= pe->q;  ep[1][3] /= pe->q;
      ep[2][0] /= pe->q;  ep[2][1] /= pe->q; ep[2][2] /= pe->q;  ep[2][3] /= pe->q;
      

      s1  = pe->d[0][0]*drx + pe->d[0][1]*dry + pe->d[0][2]*drz;  
      s2  = pe->d[1][0]*drx + pe->d[1][1]*dry + pe->d[1][2]*drz;
      s3  = pe->d[2][0]*drx + pe->d[2][1]*dry + pe->d[2][2]*drz;                      
      b1  = 2*(pe->e[0][0]*dq0 + pe->e[0][1]*dqx + pe->e[0][2]*dqy + pe->e[0][3]*dqz);
      b2  = 2*(pe->e[1][0]*dq0 + pe->e[1][1]*dqx + pe->e[1][2]*dqy + pe->e[1][3]*dqz);
      b3  = 2*(pe->e[2][0]*dq0 + pe->e[2][1]*dqx + pe->e[2][2]*dqy + pe->e[2][3]*dqz);

      /* Direct shear and bending forces */
      f1  = rparams->cs1*(s1-rparams->s01);
      f2  = rparams->cs2*(s2-rparams->s02);
      f3  = rparams->cs3*(s3-rparams->s03);
      t1  = rparams->cb1*(b1-rparams->b01);
      t2  = rparams->cb2*(b2-rparams->b02);
      t3  = rparams->cb3*(b3-rparams->b03);
      
/* Gamma X F */
      tx0 = s1*f1 + s2*f2 + s3*f3;                /* Include dot product for projection */
      tx1 = s2*f3 - s3*f2;
      tx2 = s3*f1 - s1*f3;
      tx3 = s1*f2 - s2*f1;

/* Rotate forces to natural coordinate systems */
      fx  = pe->d[0][0]*f1 + pe->d[1][0]*f2 + pe->d[2][0]*f3;                  /* dF/ds */
      fy  = pe->d[0][1]*f1 + pe->d[1][1]*f2 + pe->d[2][1]*f3;
      fz  = pe->d[0][2]*f1 + pe->d[1][2]*f2 + pe->d[2][2]*f3;

      pn1->fx += fx/le;
      pn1->fy += fy/le;
      pn1->fz += fz/le;
      pn2->fx -= fx/le;
      pn2->fy -= fy/le;
      pn2->fz -= fz/le;

      t0  = 2*(pe->e[0][0]*t1 + pe->e[1][0]*t2 + pe->e[2][0]*t3)/le;           /* dT/ds */
      tx  = 2*(pe->e[0][1]*t1 + pe->e[1][1]*t2 + pe->e[2][1]*t3)/le;
      ty  = 2*(pe->e[0][2]*t1 + pe->e[1][2]*t2 + pe->e[2][2]*t3)/le;
      tz  = 2*(pe->e[0][3]*t1 + pe->e[1][3]*t2 + pe->e[2][3]*t3)/le;

      pn1->t0 += t0;
      pn1->tx += tx;
      pn1->ty += ty;
      pn1->tz += tz;
      pn2->t0 -= t0;
      pn2->tx -= tx;
      pn2->ty -= ty;
      pn2->tz -= tz;

      t0  = (pe->e[0][0]*tx1 + pe->e[1][0]*tx2 + pe->e[2][0]*tx3)/pe->q; // Gamma X F
      tx  = (pe->e[0][1]*tx1 + pe->e[1][1]*tx2 + pe->e[2][1]*tx3)/pe->q;
      ty  = (pe->e[0][2]*tx1 + pe->e[1][2]*tx2 + pe->e[2][2]*tx3)/pe->q;
      tz  = (pe->e[0][3]*tx1 + pe->e[1][3]*tx2 + pe->e[2][3]*tx3)/pe->q;

      pn1->t0 += t0;
      pn1->tx += tx;
      pn1->ty += ty;
      pn1->tz += tz;
      pn2->t0 += t0;
      pn2->tx += tx;
      pn2->ty += ty;
      pn2->tz += tz;

      t0  = ep[0][0]*t1 + ep[1][0]*t2 + ep[2][0]*t3;               /* Omega X T */
      tx  = ep[0][1]*t1 + ep[1][1]*t2 + ep[2][1]*t3;
      ty  = ep[0][2]*t1 + ep[1][2]*t2 + ep[2][2]*t3;
      tz  = ep[0][3]*t1 + ep[1][3]*t2 + ep[2][3]*t3;

      tdq = t0*pe->q0 + tx*pe->qx + ty*pe->qy + tz*pe->qz;
      t0  -= pe->q0*tdq;
      tx  -= pe->qx*tdq;
      ty  -= pe->qy*tdq;
      tz  -= pe->qz*tdq;

      pn1->t0 += t0;
      pn1->tx += tx;
      pn1->ty += ty;
      pn1->tz += tz;
      pn2->t0 += t0;
      pn2->tx += tx;
      pn2->ty += ty;
      pn2->tz += tz;
    }
  
/* Rotate forces and torques to body frame */

  for (m = 0; m < rod->n_nodes; m++)
    {
      pn  = &rod->nodes[m];
      fx = pn->fx; 
      fy = pn->fy; 
      fz = pn->fz;
      t0 = pn->t0; 
      tx = pn->tx; 
      ty = pn->ty; 
      tz = pn->tz;
      
      pn->f1 =  pn->d[0][0]*fx + pn->d[0][1]*fy + pn->d[0][2]*fz;
      pn->f2 =  pn->d[1][0]*fx + pn->d[1][1]*fy + pn->d[1][2]*fz;
      pn->f3 =  pn->d[2][0]*fx + pn->d[2][1]*fy + pn->d[2][2]*fz;
      pn->t1 = (pn->e[0][0]*t0 + pn->e[0][1]*tx + pn->e[0][2]*ty + pn->e[0][3]*tz)*0.5;
      pn->t2 = (pn->e[1][0]*t0 + pn->e[1][1]*tx + pn->e[1][2]*ty + pn->e[1][3]*tz)*0.5;
      pn->t3 = (pn->e[2][0]*t0 + pn->e[2][1]*tx + pn->e[2][2]*ty + pn->e[2][3]*tz)*0.5;
    }

  centroforce(rod, rparams);

#ifdef motorforceON
  motorforce(rod, rparams);
#endif

#ifdef excludedvolON
  exvolsquare(rod, rparams);
#endif
  
  // BC at 0
  /*
    pn  = &rod->nodes[0];
      
    pn->f1 = 0;
    pn->f2 = 0;
    pn->f3 = 0;        
    pn->t1 = 0;
    pn->t2 = 0;
    pn->t3 = 0;
  */
  // BC at  L

}
