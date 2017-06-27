#include "rod_proto.h"

void euler(struct rod *rod, struct rodparams *rparams, real dt)
{
  struct node *pn;
  int  m;
  
  for (m = 0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];
      pn->p1 = rparams->xm*pn->f1/rparams->zt1;
      pn->p2 = rparams->xm*pn->f2/rparams->zt2;
      pn->p3 = rparams->xm*pn->f3/rparams->zt3;
    }
  
  for (m = 0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];
      pn->px = pn->d[0][0]*pn->p1 + pn->d[1][0]*pn->p2 + pn->d[2][0]*pn->p3;
      pn->py = pn->d[0][1]*pn->p1 + pn->d[1][1]*pn->p2 + pn->d[2][1]*pn->p3;
      pn->pz = pn->d[0][2]*pn->p1 + pn->d[1][2]*pn->p2 + pn->d[2][2]*pn->p3;
    }

  for (m = 0; m < rod->n_nodes; m++)
    {      
      pn = &rod->nodes[m];
      pn->rx += pn->px*dt/rparams->xm;
      pn->ry += pn->py*dt/rparams->xm;
      pn->rz += pn->pz*dt/rparams->xm;
    }
}
