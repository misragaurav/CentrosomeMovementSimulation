#include "rod_proto.h"

void q(struct rod *rod, struct rodparams *rparams, real dt, int dir)
{
  struct node  *pn;
  real l[3], q[4];
  real i[3];
  real eq[4], c, s, c2, s2, halfangle;
  real l1_temp, l2_temp;           
  int  m, n, dir1, dir2;

  i[0] = rparams->i1;  i[1] = rparams->i2;  i[2] = rparams->i3;

  for (m = 0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];      
      q[0] = pn->q0;  q[1] = pn->qx;  q[2] = pn->qy; q[3] = pn->qz;
      l[0] = pn->l1;  l[1] = pn->l2;  l[2] = pn->l3;

      switch (dir)
	{
	case 0: eq[0] =-q[1]; eq[1] = q[0]; eq[2] = q[3]; eq[3] = -q[2]; break;
	case 1: eq[0] =-q[2]; eq[1] =-q[3]; eq[2] = q[0]; eq[3] =  q[1]; break;
	case 2: eq[0] =-q[3]; eq[1] = q[2]; eq[2] =-q[1]; eq[3] =  q[0]; break;
	}

      halfangle = 0.5*( (l[dir])/(i[dir]) )*dt;

      c = cos(halfangle);  s = sin(halfangle);

      for (n = 0; n < 4; n++)
	{ q[n] = c*q[n] + s*eq[n]; }
      pn->q0 = q[0];
      pn->qx = q[1];
      pn->qy = q[2];
      pn->qz = q[3];

      c2 = cos(2*halfangle); s2 = sin(2*halfangle);

      dir1 = (dir+1)%3;
      dir2 = (dir+2)%3;
      
      l1_temp = l[dir1];                   // stored old values of momenta
      l2_temp = l[dir2];
      l[dir2] = l2_temp*c2 - l1_temp*s2;   // Rotating the frames
      l[dir1] = l2_temp*s2 + l1_temp*c2;
      
      pn->l1 = l[0];
      pn->l2 = l[1];
      pn->l3 = l[2];
    }
}

void r(struct rod *rod, struct rodparams *rparams, real dt)
{
  struct node *pn;
  real xm;
  int m;

  xm = rparams->xm;
  for (m = 0; m < rod->n_nodes; m++)
    {      
      pn = &rod->nodes[m];
      pn->rx += pn->px*dt/xm;
      pn->ry += pn->py*dt/xm;
      pn->rz += pn->pz*dt/xm;
    }
}

void momenta(struct rod *rod, struct rodparams *rparams, real dt)
{
  struct node *pn;
  int  m;
  
  for (m = 0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];

      pn->p1 = pn->d[0][0]*pn->px + pn->d[0][1]*pn->py + pn->d[0][2]*pn->pz;
      pn->p2 = pn->d[1][0]*pn->px + pn->d[1][1]*pn->py + pn->d[1][2]*pn->pz;
      pn->p3 = pn->d[2][0]*pn->px + pn->d[2][1]*pn->py + pn->d[2][2]*pn->pz;

      /*
      pn->p1 = (pn->p1*(1 - (dt/2.0)*pn->zt1/rparams->xm) + pn->f1*dt)/(1 + (dt/2.0)*pn->zt1/rparams->xm);
      pn->p2 = (pn->p2*(1 - (dt/2.0)*pn->zt2/rparams->xm) + pn->f2*dt)/(1 + (dt/2.0)*pn->zt2/rparams->xm);
      pn->p3 = (pn->p3*(1 - (dt/2.0)*pn->zt3/rparams->xm) + pn->f3*dt)/(1 + (dt/2.0)*pn->zt3/rparams->xm);
      */

      pn->p1 = (pn->p1*(1 - (dt/2.0)*rparams->zt1/rparams->xm) + pn->f1*dt)/(1 + (dt/2.0)*rparams->zt1/rparams->xm);
      pn->p2 = (pn->p2*(1 - (dt/2.0)*rparams->zt2/rparams->xm) + pn->f2*dt)/(1 + (dt/2.0)*rparams->zt2/rparams->xm);
      pn->p3 = (pn->p3*(1 - (dt/2.0)*rparams->zt3/rparams->xm) + pn->f3*dt)/(1 + (dt/2.0)*rparams->zt3/rparams->xm);

      pn->l1 = (pn->l1*(1 - (dt/2.0)*rparams->zr1/rparams->i1) + pn->t1*dt)/(1 + (dt/2.0)*rparams->zr1/rparams->i1);
      pn->l2 = (pn->l2*(1 - (dt/2.0)*rparams->zr2/rparams->i2) + pn->t2*dt)/(1 + (dt/2.0)*rparams->zr2/rparams->i2);
      pn->l3 = (pn->l3*(1 - (dt/2.0)*rparams->zr3/rparams->i3) + pn->t3*dt)/(1 + (dt/2.0)*rparams->zr3/rparams->i3);

      pn->px = pn->d[0][0]*pn->p1 + pn->d[1][0]*pn->p2 + pn->d[2][0]*pn->p3;
      pn->py = pn->d[0][1]*pn->p1 + pn->d[1][1]*pn->p2 + pn->d[2][1]*pn->p3;
      pn->pz = pn->d[0][2]*pn->p1 + pn->d[1][2]*pn->p2 + pn->d[2][2]*pn->p3;
    }
}


void propagator(struct rod *rod, struct rodparams *rparams, real dt)
{
  q(rod, rparams, 0.5*dt/(2), 2);  // Free rotation (L3)
  q(rod, rparams, 0.5*dt/(2), 1);  // Free rotation (L2)
  q(rod, rparams, 0.5*dt/(2), 0);  // Free rotation (L1)
  q(rod, rparams, 0.5*dt/(2), 0);  // Free rotation (L1)
  q(rod, rparams, 0.5*dt/(2), 1);  // Free rotation (L2)
  q(rod, rparams, 0.5*dt/(2), 2);  // Free rotation (L3)

  rotmatrix_node(rod);
  rotmatrix_elem(rod);
  
  r(rod, rparams, dt/2);

  force(rod, rparams);
  momenta(rod, rparams, dt);
  
  r(rod, rparams, dt/2);

  q(rod, rparams, 0.5*dt/(2), 2);  // Free rotation (L3) //
  q(rod, rparams, 0.5*dt/(2), 1);  // Free rotation (L2) // 
  q(rod, rparams, 0.5*dt/(2), 0);  // Free rotation (L1) //
  q(rod, rparams, 0.5*dt/(2), 0);  // Free rotation (L1) // 
  q(rod, rparams, 0.5*dt/(2), 1);  // Free rotation (L2) // 
  q(rod, rparams, 0.5*dt/(2), 2);  // Free rotation (L3) // 

  rotmatrix_node(rod);
  rotmatrix_elem(rod);

}
