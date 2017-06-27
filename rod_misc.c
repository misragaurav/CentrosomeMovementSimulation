#include "rod_proto.h"

void file_name (char *name, char *work_dir, int task_number)
{
  static char   name0[] = {'0','1','2','3','4','5','6','7','8','9'};
  int    index;

  if (work_dir != '\0')
    {
      for (index = strlen(name); index > 0; index--)
	name[index+strlen(work_dir)]=name[index-1];
      name[strlen(work_dir)] = '/';
      for (index = 0; index < strlen(work_dir); index++)
	name[index] = work_dir[index];
    }

  index = strlen(name);
  name[index]   = name0[task_number/100];
  task_number %= 100;
  name[index+1] = name0[task_number/10];
  task_number %= 10;
  name[index+2] = name0[task_number];
  name[index+3] = '.';
  name[index+4] = 'o';
  name[index+5] = 'd';
  name[index+6] = 's';
  name[index+7] = '\0';
}

void rotmatrix_node(struct rod *rod)
{
  struct node *pn;
  real q0, qx, qy, qz;
  int m;

  for (m = 0; m < rod->n_nodes; m++)
    {
      pn  = &rod->nodes[m];
      q0 = pn->q0; qx = pn->qx; qy = pn->qy; qz = pn->qz;
      
      pn->d[0][0] = -qy*qy+qx*qx-qz*qz+q0*q0; pn->d[0][1] = 2*qy*qx+2*qz*q0;          pn->d[0][2] = -2*qy*q0+2*qx*qz;
      pn->d[1][0] = 2*qy*qx-2*qz*q0;          pn->d[1][1] = qy*qy-qx*qx-qz*qz+q0*q0;  pn->d[1][2] = 2*qy*qz+2*qx*q0;
      pn->d[2][0] = 2*qy*q0+2*qx*qz;          pn->d[2][1] = 2*qy*qz-2*qx*q0;          pn->d[2][2] = -qy*qy-qx*qx+qz*qz+q0*q0;
      pn->e[0][0] = -qx;  pn->e[0][1] =  q0; pn->e[0][2] =  qz;  pn->e[0][3] = -qy;
      pn->e[1][0] = -qy;  pn->e[1][1] = -qz; pn->e[1][2] =  q0;  pn->e[1][3] =  qx;
      pn->e[2][0] = -qz;  pn->e[2][1] =  qy; pn->e[2][2] = -qx;  pn->e[2][3] =  q0;
    }
}

void rotmatrix_elem(struct rod *rod)
{
  struct node *pn1, *pn2;
  struct elem *pe;
  real q0, qx, qy, qz, q2, q;
  int m;

  for (m = 1; m < rod->n_elems-1; m++)
    {
      pe  = &rod->elems[m];              
      pn1 = &rod->nodes[m-1];
      pn2 = &rod->nodes[m];

      q0 = 0.5*(pn1->q0 + pn2->q0);
      qx = 0.5*(pn1->qx + pn2->qx);
      qy = 0.5*(pn1->qy + pn2->qy);
      qz = 0.5*(pn1->qz + pn2->qz);
      q2 = q0*q0 + qx*qx + qy*qy + qz*qz;

      q   = sqrt(q2);
      q0 /= q; qx /= q; qy /= q; qz /= q;

      pe->q0=q0; pe->qx=qx; pe->qy=qy; pe->qz=qz;
      pe->q =q ; pe->q2=q2;

      pe->d[0][0] = -qy*qy+qx*qx-qz*qz+q0*q0; pe->d[0][1] = 2*qy*qx+2*qz*q0;          pe->d[0][2] = -2*qy*q0+2*qx*qz;
      pe->d[1][0] = 2*qy*qx-2*qz*q0;          pe->d[1][1] = qy*qy-qx*qx-qz*qz+q0*q0;  pe->d[1][2] = 2*qy*qz+2*qx*q0;
      pe->d[2][0] = 2*qy*q0+2*qx*qz;          pe->d[2][1] = 2*qy*qz-2*qx*q0;          pe->d[2][2] = -qy*qy-qx*qx+qz*qz+q0*q0;
      pe->e[0][0] = -qx;  pe->e[0][1] =  q0; pe->e[0][2] =  qz;  pe->e[0][3] = -qy;
      pe->e[1][0] = -qy;  pe->e[1][1] = -qz; pe->e[1][2] =  q0;  pe->e[1][3] =  qx;
      pe->e[2][0] = -qz;  pe->e[2][1] =  qy; pe->e[2][2] = -qx;  pe->e[2][3] =  q0;
    }
}

void com(struct rod *rod, struct rodparams *rparams)
{
  struct node *pn;
  real Rx=0, Ry=0, Rz=0;
  int k;

  rod->xcm  = 0.0;
  rod->ycm  = 0.0;
  rod->zcm  = 0.0;
  rod->pxcm = 0.0;
  rod->pycm = 0.0;
  rod->pzcm = 0.0;
  rod->lxcm = 0.0;
  rod->lycm = 0.0;
  rod->lzcm = 0.0;

  for (k = 0; k < rod->n_nodes; k++)
    {
      pn = &rod[0].nodes[k];
      rod->xcm += pn->rx/(rod->n_nodes);
      rod->ycm += pn->ry/(rod->n_nodes);
      rod->zcm += pn->rz/(rod->n_nodes);
    }
    
  for (k = 0; k < rod->n_nodes; k++)
    {
      pn = &rod->nodes[k];

      rod->pxcm += rparams->le*pn->px;
      rod->pycm += rparams->le*pn->py;
      rod->pzcm += rparams->le*pn->pz;
    }

  for (k = 0; k < rod->n_nodes; k++)
    {
      pn = &rod->nodes[k];
      pn->lx = pn->d[0][0]*pn->l1 + pn->d[1][0]*pn->l2 + pn->d[2][0]*pn->l3;
      pn->ly = pn->d[0][1]*pn->l1 + pn->d[1][1]*pn->l2 + pn->d[2][1]*pn->l3;
      pn->lz = pn->d[0][2]*pn->l1 + pn->d[1][2]*pn->l2 + pn->d[2][2]*pn->l3;
    }

  for (k = 0; k < rod->n_nodes; k++)
    {
      pn = &rod->nodes[k];
	  
      Rx = pn->rx - rod->xcm;
      Ry = pn->ry - rod->ycm;
      Rz = pn->rz - rod->zcm;

      rod->lxcm += rparams->le*(pn->lx + Ry*pn->pz - Rz*pn->py);
      rod->lycm += rparams->le*(pn->ly + Rz*pn->px - Rx*pn->pz);
      rod->lzcm += rparams->le*(pn->lz + Rx*pn->py - Ry*pn->px);
    }
}


void energy(struct rod *rod, struct rodparams *rparams, int t)
{
  struct node *pn1, *pn2;
  struct elem *pe;
  real drx, dry, drz, dq0, dqx, dqy, dqz;
  real s1, s2, s3, b1, b2, b3;
  real p1, p2, p3, l1, l2, l3;
  real xm, i1, i2, i3;
  int m;
  char   filename[64];
  FILE   *fileptr, *fileptr1, *fileptr2;

  strcpy(filename, datadir);
  strcat(filename, "/e.ods");
  fileptr = fopen (filename, "a");

  strcpy(filename,datadir);
  strcat(filename, "/ke.ods");
  fileptr1 = fopen (filename, "a");

  strcpy(filename,datadir);
  strcat(filename, "/pe.ods");
  fileptr2 = fopen (filename, "a");
  
  xm = rparams->xm;  i1 = rparams->i1;  i2 = rparams->i2;  i3 = rparams->i3;

  for (m = 0; m < rod->n_nodes; m++)
    {
      pn1 = &rod->nodes[m];

      p1 = pn1->d[0][0]*pn1->px + pn1->d[0][1]*pn1->py + pn1->d[0][2]*pn1->pz;
      p2 = pn1->d[1][0]*pn1->px + pn1->d[1][1]*pn1->py + pn1->d[1][2]*pn1->pz;
      p3 = pn1->d[2][0]*pn1->px + pn1->d[2][1]*pn1->py + pn1->d[2][2]*pn1->pz;
      l1 = pn1->l1;  
      l2 = pn1->l2;  
      l3 = pn1->l3;

      pn1->keT1 = p1*p1/(2*xm);
      pn1->keT2 = p2*p2/(2*xm);
      pn1->keT3 = p3*p3/(2*xm);	       
      pn1->keR1 = l1*l1/(2*i1);
      pn1->keR2 = l2*l2/(2*i2);
      pn1->keR3 = l3*l3/(2*i3);

      pn1->ke = pn1->keT1 + pn1->keT2 + pn1->keT3 + pn1->keR1 + pn1->keR2 + pn1->keR3 ;
    }

  for (m = 1; m < rod->n_elems-1; m++)
    {
      pe  = &rod->elems[m];              
      pn1 = &rod->nodes[m-1];
      pn2 = &rod->nodes[m];

      drx = (pn2->rx - pn1->rx)/rparams->le;
      dry = (pn2->ry - pn1->ry)/rparams->le;
      drz = (pn2->rz - pn1->rz)/rparams->le;
      dq0 = (pn2->q0 - pn1->q0)/rparams->le;
      dqx = (pn2->qx - pn1->qx)/rparams->le;
      dqy = (pn2->qy - pn1->qy)/rparams->le;
      dqz = (pn2->qz - pn1->qz)/rparams->le;

      s1  = pe->d[0][0]*drx + pe->d[0][1]*dry + pe->d[0][2]*drz;
      s2  = pe->d[1][0]*drx + pe->d[1][1]*dry + pe->d[1][2]*drz;
      s3  = pe->d[2][0]*drx + pe->d[2][1]*dry + pe->d[2][2]*drz;                      
      b1  = 2*(pe->e[0][0]*dq0 + pe->e[0][1]*dqx + pe->e[0][2]*dqy + pe->e[0][3]*dqz);
      b2  = 2*(pe->e[1][0]*dq0 + pe->e[1][1]*dqx + pe->e[1][2]*dqy + pe->e[1][3]*dqz);
      b3  = 2*(pe->e[2][0]*dq0 + pe->e[2][1]*dqx + pe->e[2][2]*dqy + pe->e[2][3]*dqz);

      s1 -= rparams->s01;
      s2 -= rparams->s02;
      s3 -= rparams->s03;
      b1 -= rparams->b01;
      b2 -= rparams->b02;
      b3 -= rparams->b03;

      pe->peT1 = 0.5*s1*s1*rparams->cs1;
      pe->peT2 = 0.5*s2*s2*rparams->cs2;
      pe->peT3 = 0.5*s3*s3*rparams->cs3;
      pe->peR1 = 0.5*b1*b1*rparams->cb1;
      pe->peR2 = 0.5*b2*b2*rparams->cb2;
      pe->peR3 = 0.5*b3*b3*rparams->cb3;
      
      pe->pe = pe->peT1+pe->peT2+pe->peT3+pe->peR1+pe->peR2+pe->peR3;
    }

  for (m=1, rod->pe=0, rod->peT1=0, rod->peT2=0, rod->peT3=0, rod->peR1=0, rod->peR2=0, rod->peR3=0; m < rod->n_elems-1; m++)
    {
      pe = &rod->elems[m];
      rod->pe   += pe->pe;
      rod->peT1 += pe->peT1 ;
      rod->peT2 += pe->peT2 ;
      rod->peT3 += pe->peT3 ;
      rod->peR1 += pe->peR1 ;
      rod->peR2 += pe->peR2 ;
      rod->peR3 += pe->peR3 ;
    }

  for (m=0, rod->ke=0, rod->keT1=0, rod->keT2=0, rod->keT3=0, rod->keR1=0, rod->keR2=0, rod->keR3=0; m < rod->n_nodes; m++)
    {
      pn1 = &rod->nodes[m];
      rod->ke   += pn1->ke;
      rod->keT1 += pn1->keT1 ;
      rod->keT2 += pn1->keT2 ;
      rod->keT3 += pn1->keT3 ;
      rod->keR1 += pn1->keR1 ;
      rod->keR2 += pn1->keR2 ;
      rod->keR3 += pn1->keR3 ;
    }
  
  fprintf (fileptr,  "%f % .5e % .5e % .5e\n", t*rparams->dt/1000, rod->pe, rod->ke, rod->pe + rod->ke);
  fprintf (fileptr1, "%f % .5e % .5e % .5e % .5e % .5e % .5e\n", t*rparams->dt/1000, rod->keT1, rod->keT2, rod->keT3, rod->keR1, rod->keR2, rod->keR3);
  fprintf (fileptr2, "%f % .5e % .5e % .5e % .5e % .5e % .5e\n", t*rparams->dt/1000, rod->peT1, rod->peT2, rod->peT3, rod->peR1, rod->peR2, rod->peR3);

  fclose (fileptr);
  fclose (fileptr1);
  fclose (fileptr2);

}
