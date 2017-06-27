#include "proto.h"

#define animation
#define animation_small

void makedatadir(char **argv)
{
  strcpy(datadir,"data");

  strcat(datadir,"-");
  strcat(datadir,argv[1]);
  strcat(datadir,"-");
  strcat(datadir,argv[2]);
  strcat(datadir,"-");
  strcat(datadir,argv[3]);
  strcat(datadir,"-");
  strcat(datadir,argv[4]);
  strcat(datadir,"-");
  strcat(datadir,argv[5]);

  char command[64]={"mkdir "};
  strcat(command,datadir);
  if (system(command) != 0) // check if the command was successful
    {printf("mkdir not successful, bailing out \n");exit(1);}
}

void opendatadir(char **argv)
{
  strcpy(datadir,"data");

  strcat(datadir,"-");
  strcat(datadir,argv[1]);
  strcat(datadir,"-");
  strcat(datadir,argv[2]);
  strcat(datadir,"-");
  strcat(datadir,argv[3]);
  strcat(datadir,"-");
  strcat(datadir,argv[4]);
  strcat(datadir,"-");
  strcat(datadir,argv[5]);

  char command[64]={"ls "};
  strcat(command,datadir);
  if (system(command) != 0) // check if the command was successful
    {printf("datadir does not exist, bailing out \n");exit(1);}

}

void print(struct mt *mt, struct centro *centro, int Nmt, int t, real dt)
{
  struct seg *ps; 
  FILE  *file_ptr=0;
  int n, k;
  int div; 
  real angle;
  int totalsegs;

#ifdef animation
  char animationxyz[64];
  strcpy(animationxyz, datadir);
  strcat(animationxyz, "/a.xyz");
  file_ptr  = fopen(animationxyz, "a");
 
#ifdef animation_small		
  if (t%(20*iprint) == 0)	
#endif
    {
      totalsegs = 0;
      for(n=0; n<Nmt; n++)
	totalsegs = totalsegs + mt[n].N;

      fprintf (file_ptr,"%3i \nAtom_name \n", totalsegs);

      totalsegs = 0; //reset
      for(n=0; n<Nmt; n++)
	{
	  for (k = 0; k < mt[n].N; k++)
	    {
	      ps = &(mt[n].seg[k]);
	      totalsegs = totalsegs + 1;
	      fprintf (file_ptr, "AtomNum%i \t %f \t %f \t %f \n", totalsegs, ps->rx, ps->ry, ps->rz);
	    }
	}
    }
  fclose (file_ptr);
#endif

  char centrofile[64];
  strcpy(centrofile, datadir);
  strcat(centrofile, "/c.xyz");
  file_ptr = fopen(centrofile, "a");
  fprintf (file_ptr,"%i \nAtom_name \n", 1);
  fprintf (file_ptr,"AtomNum%i %f \t %f \t %f \n", 1, centro->rx, centro->ry, centro->rz);
  fclose (file_ptr);

  char centrotime[64];
  strcpy(centrotime, datadir);
  strcat(centrotime, "/c.ods");
  file_ptr = fopen(centrotime, "a");
  fprintf (file_ptr,"%f %.16e \t %.16e \t %.16e \n", t*dt, centro->rx, centro->ry, centro->rz);
  fclose (file_ptr);

  //print circular boundary
  /*
  div = 16; angle = 0;

  if (t%(20*iprint) == 0)
    {
      char boundary[64];
      strcpy(boundary, datadir);
      strcat(boundary, "/b.xyz");
      file_ptr = fopen(boundary, "a");
      fprintf (file_ptr,"%i \nAtom_name \n", 2*div);
      for(k=0; k<div; k++)
	{       
	  angle = angle + 2*Pi/div;
	  fprintf (file_ptr,"AtomNum%i %f \t %f \t %f \n", k+1, (Dia/2.0)*cos(angle), (Dia/2.0)*sin(angle), thickness/2.0);
	  fprintf (file_ptr,"AtomNum%i %f \t %f \t %f \n", k+1+div, (Dia/2.0)*cos(angle), (Dia/2.0)*sin(angle), -thickness/2.0);
	}

      fclose(file_ptr);
    }
  */


  //print square boundary

  if (t%(20*iprint) == 0)
    {
      char boundary[64];
      strcpy(boundary, datadir);
      strcat(boundary, "/b.xyz");
      file_ptr = fopen(boundary, "a");
      fprintf (file_ptr,"%i \nAtom_name \n", 8);
      fprintf (file_ptr,"AtomNum%i %f \t %f \t %f \n", 1, +Dia/2.0, +Dia/2.0, thickness/2.0);
      fprintf (file_ptr,"AtomNum%i %f \t %f \t %f \n", 2, -Dia/2.0, +Dia/2.0, thickness/2.0);
      fprintf (file_ptr,"AtomNum%i %f \t %f \t %f \n", 3, +Dia/2.0, -Dia/2.0, thickness/2.0);
      fprintf (file_ptr,"AtomNum%i %f \t %f \t %f \n", 4, -Dia/2.0, -Dia/2.0, thickness/2.0);

      fprintf (file_ptr,"AtomNum%i %f \t %f \t %f \n", 5, +Dia/2.0, +Dia/2.0, -thickness/2.0);
      fprintf (file_ptr,"AtomNum%i %f \t %f \t %f \n", 6, -Dia/2.0, +Dia/2.0, -thickness/2.0);
      fprintf (file_ptr,"AtomNum%i %f \t %f \t %f \n", 7, +Dia/2.0, -Dia/2.0, -thickness/2.0);
      fprintf (file_ptr,"AtomNum%i %f \t %f \t %f \n", 8, -Dia/2.0, -Dia/2.0, -thickness/2.0);

      fclose(file_ptr);
    }

}

void checkpoint(char *checkpointfile, struct mt *mt, struct centro *centro, int Nmt, int t, real dt)
{
  struct seg *ps; 
  FILE  *file_ptr=0;
  int n, k;

  file_ptr = fopen(checkpointfile, "wb");

  gsl_rng_fwrite(file_ptr, rng);

  fwrite(&t, sizeof(int), 1, file_ptr);

  fwrite(&centro->rx, sizeof(real), 1, file_ptr);
  fwrite(&centro->ry, sizeof(real), 1, file_ptr);
  fwrite(&centro->rz, sizeof(real), 1, file_ptr);

  for(n=0; n<Nmt; n++)
    {
      fwrite(&mt[n].id,     sizeof(int),  1, file_ptr);
      fwrite(&mt[n].N,      sizeof(int),  1, file_ptr);
      fwrite(&mt[n].status, sizeof(int),  1, file_ptr);
      fwrite(&mt[n].pinned, sizeof(int),  1, file_ptr);
      fwrite(&mt[n].dL,     sizeof(real), 1, file_ptr);
    }

  for(n=0; n<Nmt; n++)
    {
      for (k = 0; k < mt[n].N; k++)
	{
	  ps = &(mt[n].seg[k]);
	  fwrite(&ps->rx, sizeof(real), 1, file_ptr);
	  fwrite(&ps->ry, sizeof(real), 1, file_ptr);
	  fwrite(&ps->rz, sizeof(real), 1, file_ptr);

	  fwrite(&ps->q0, sizeof(real), 1, file_ptr);
	  fwrite(&ps->qx, sizeof(real), 1, file_ptr);
	  fwrite(&ps->qy, sizeof(real), 1, file_ptr);
	  fwrite(&ps->qz, sizeof(real), 1, file_ptr);
	}
    } 
  fclose(file_ptr);
}

void randomorient(struct mt *mt)
{
  struct seg *ps;
  real th=0, ph=0, si=0;
  real q, q2;
  real x=1, y=1, z=1;

  ps =&mt->seg[0];  //// Always 0th node because randomorient() is called only for the first node.

  x=1; y=1; z=1;
  while(x*x+y*y+z*z>=1)
    {
      // generate points in a cube
      x = 2*gsl_rng_uniform(rng)-1; // [-1 +1)
      y = 2*gsl_rng_uniform(rng)-1; // [-1 +1)
      z = 2*gsl_rng_uniform(rng)-1; // [-1 +1)
      // accept those which lie in the enclosed sphere
    }
      
  th = acos(z);
  
  if(atan2(y, x) < 0)
    ph = atan2(y, x) + 2*Pi; // to convert from (-pi pi) to (0 2pi)
  else 
    ph = atan2(y, x);

  si = 0; // only two angles are needed for defining the orientation of a vector      

  //////////
  /*
    if(2*gsl_rng_uniform(rng)-1 > 0)
    ph = 3*Pi/2.0;
    else
    ph = Pi/2.0;

    printf("%f \n", ph);

    th = Pi/2.0;  
    si = 0;
  */
  //////////

  ps->q0 = cos(th/2.0)*cos((ph+si)/2.0);
  ps->qx = sin(th/2.0)*cos((ph-si)/2.0);
  ps->qy = sin(th/2.0)*sin((ph-si)/2.0);
  ps->qz = cos(th/2.0)*sin((ph+si)/2.0);
      
  q2 = ps->q0*ps->q0 + ps->qx*ps->qx + ps->qy*ps->qy + ps->qz*ps->qz;
  q  = sqrt(q2);

  ps->q0 = ps->q0/q;
  ps->qx = ps->qx/q;
  ps->qy = ps->qy/q;
  ps->qz = ps->qz/q;
}

void randomorientcone(struct mt *mt)
{
  struct seg *ps;
  real th=0, ph=0, si=0;
  real q, q2;
  real x=1, y=1, z=1;

  ps =&mt->seg[0];  //// Always 0th node because randomorient() is called only for the first node.

  while(th < atan2(Dia, 2*thickness) || th > Pi - atan2(Dia, 2*thickness))
    {
      x=1; y=1; z=1;
      while(x*x+y*y+z*z>=1)
	{
	  // generate points in a cube
	  x = 2*gsl_rng_uniform(rng)-1; // [-1 +1)
	  y = 2*gsl_rng_uniform(rng)-1; // [-1 +1)
	  z = 2*gsl_rng_uniform(rng)-1; // [-1 +1)
	  // accept those which lie in the enclosed sphere
	}
      th = acos(z);
      //      printf("th: %f \t %f \t %f \n", th, atan2(Dia, 2*thickness), Pi - atan2(Dia, 2*thickness));      
    }
  
  if(atan2(y, x) < 0)
    ph = atan2(y, x) + 2*Pi; // to convert from (-pi pi) to (0 2pi)
  else 
    ph = atan2(y, x);

  si = 0; // only two angles are needed for defining the orientation of a vector    

  ps->q0 = cos(th/2.0)*cos((ph+si)/2.0);
  ps->qx = sin(th/2.0)*cos((ph-si)/2.0);
  ps->qy = sin(th/2.0)*sin((ph-si)/2.0);
  ps->qz = cos(th/2.0)*sin((ph+si)/2.0);
      
  q2 = ps->q0*ps->q0 + ps->qx*ps->qx + ps->qy*ps->qy + ps->qz*ps->qz;
  q  = sqrt(q2);

  ps->q0 = ps->q0/q;
  ps->qx = ps->qx/q;
  ps->qy = ps->qy/q;
  ps->qz = ps->qz/q;
}

void rotmatrix(struct seg *ps)
{
  real q0, qx, qy, qz;

  q0 = ps->q0; qx = ps->qx; qy = ps->qy; qz = ps->qz;

  ps->R[0][0] =  -qy*qy+qx*qx-qz*qz+q0*q0; ps->R[0][1] = 2*qy*qx+2*qz*q0;         ps->R[0][2] =  -2*qy*q0+2*qx*qz;
  ps->R[1][0] =  2*qy*qx-2*qz*q0;          ps->R[1][1] = qy*qy-qx*qx-qz*qz+q0*q0; ps->R[1][2] =  2*qy*qz+2*qx*q0;
  ps->R[2][0] =  2*qy*q0+2*qx*qz;          ps->R[2][1] = 2*qy*qz-2*qx*q0;         ps->R[2][2] =  -qy*qy-qx*qx+qz*qz+q0*q0;
}

void movecentro(struct params *params, struct centro *centro, struct mt *mt, real dt)
{
  struct seg *ps;
  real rx, ry, rz; 
  int n;

  centro->fx = 0; centro->fy = 0; centro->fz = 0;

  for(n=0; n<Nmt; n++)
    { 
      ps = &mt[n].seg[0];
      mt[n].fcx = params->k*(ps->rx - centro->rx); // individual MT force
      mt[n].fcy = params->k*(ps->ry - centro->ry);
      mt[n].fcz = params->k*(ps->rz - centro->rz);
    }

  for(n=0; n<Nmt; n++)
    {
      centro->fx = centro->fx + mt[n].fcx;  //sum of all MT forces
      centro->fy = centro->fy + mt[n].fcy;
      centro->fz = centro->fz + mt[n].fcz;
    }

  //  printf("%f \t %f \t %f \n", centro->fx, centro->fy, centro->fz);

  rx = centro->rx + (centro->fx/params->fric_centro)*dt;
  ry = centro->ry + (centro->fy/params->fric_centro)*dt;
  rz = centro->rz + (centro->fz/params->fric_centro)*dt;

  if(checkboundary(rx, ry, rz, Dia, thickness) == 0)
    {
      centro->rx = rx;
      centro->ry = ry;
      centro->rz = rz;
    }

}
