#include "rod_proto.h"

//#define debug

void rodflex(struct rodparams *rparams, struct rod *rod, int t)
{
  struct node *pn;
  real suml1 = 0, suml2 = 0, suml3 = 0, sump1 = 0, sump2 = 0, sump3 = 0;
  real err = 1, tol = 1e-7;
  int i, counter=0;
  real delt = 1e-5;
  real factor;

#ifdef debug
  FILE  *fileptr=0;
  char iters[64];
  strcpy(iters, datadir);
  strcat(iters, "/iters.ods");
  fileptr = fopen(iters, "a");
#endif

  rotmatrix_node(rod);
  rotmatrix_elem(rod);

  // Rolling fric
#ifdef motorforceON
  factor = 5e-5;
#else
  factor = 5e-3;
#endif
  rparams->zr1 = factor*rparams->zt1; //need critical damping for the higher rotational modes 
  rparams->zr2 = factor*rparams->zt1;                                        
  rparams->zr3 = factor*rparams->zt1;

  counter=0;
  while(err>tol)
    {
      counter++;
      
      propagator(rod, rparams, delt);

      suml1 = 0; suml2 = 0; suml3 = 0; sump1 = 0; sump2 = 0; sump3 = 0;
      for(i=0;i<rod->n_nodes;i++)
	{
	  pn = &rod->nodes[i];
	  
	  suml1 = suml1 + fabs(pn->l1);
	  suml2 = suml2 + fabs(pn->l2);
	  suml3 = suml3 + fabs(pn->l3);
	  sump1 = sump1 + fabs(pn->p1);
	  sump2 = sump2 + fabs(pn->p2);
	  sump3 = sump3 + fabs(pn->p3);
	}

      err = (suml1 + suml2 + suml3 + sump1 + sump2 + sump3)/rod->n_nodes;

#ifdef debug
      //      if (t%(iprint) == 0)
	{
	  //	  if(err<tol)
	    {
	      fprintf(fileptr, "%i ", counter);
	      
	      pn = &rod->nodes[rod->n_nodes-1];
	      //fprintf(fileptr, "%.16e %.16e \n", pn->q0 ,pn->qx );
	      //fprintf(fileptr, "%.16e %.16e ", pn->ry ,pn->rz );
	      fprintf(fileptr, "%.16e %.16e %.16e \n", pn->l1, pn->p2, pn->p3);
	      
	      //printf(fileptr, "%.16e \n", err);
	      //printf(fileptr, "\n");
	    }
	}
#endif

    }
  
  force(rod, rparams);
  euler(rod, rparams, rparams->dt);

  //  if (t%(iprint) == 0) energy(rod, rparams, t);

#ifdef debug
  fclose (fileptr);
#endif

}
