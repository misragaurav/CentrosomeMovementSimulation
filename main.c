#include "proto.h"

int main( argc, argv )
int argc;
char **argv;
{
  int rank, size; // size = number of processors in the hosts file
  real tstart, tstop;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  tstart = MPI_Wtime();
  fprintf (stdout, "Begin process %d of %d\n", rank+1, size);

  struct rodparams *rparams; struct params *params;
  struct centro *centro; struct mt *mt; struct rod *rod;
  struct seg *ps;

  real *buffer, *bufsend, *bufrecv;
  real *ptr2buf;

  char checkpointfile[64], ext[3];
  strcpy(checkpointfile,"checkpoint");
  sprintf(ext, "%03i", rank);
  strcat(checkpointfile, ext);
  strcat(checkpointfile, ".bin");

  int t, n, k, i, iprint_checkpoint;
  real dt;

  if (argc < 6)
  {printf("\n Too few arguments to the C code, bailing out. \n"); exit(0); }

  Nmt  	 = atoi(argv[1]);
  dt 	 = atof(argv[2]);
  ds     = atof(argv[3]);
  iter 	 = atoi(argv[4]);	  //number of iterations
  iprint = atoi(argv[5]);	  //iterations to save after

  Gmax = Nmt*Nmax/2;
  iprint_checkpoint = 100*iprint;

  #ifndef resume
  if (rank == 0) makedatadir(argv);
  #else
  if (rank == 0) opendatadir(argv);
  #endif

  // GSL random number generator setup
  gsl_rng_default_seed=1;
  rng = gsl_rng_alloc(gsl_rng_taus);

  //  printf ("generator type: %s \nseed = %lu \nrandom number: %f \n", gsl_rng_name(rng), gsl_rng_default_seed, gsl_rng_uniform(rng));

  // memory allocation
  params  = (struct params    *) calloc(1,   sizeof(struct params   ));
  centro  = (struct centro    *) calloc(1,   sizeof(struct centro   ));
  mt      = (struct mt        *) calloc(Nmt, sizeof(struct mt       ));
  rod     = (struct rod       *) calloc(Nmt, sizeof(struct rod      ));
  rparams = (struct rodparams *) calloc(1,   sizeof(struct rodparams));

  buffer  = calloc(Nmt*Nmax*7, sizeof(real));
  bufsend = calloc(Nmt*3, sizeof(real));
  bufrecv = calloc(Nmt*3, sizeof(real));

  #ifndef resume
  initmtcentro(mt, centro);           // initialize MT and centro positions
  t=0;
  #else
  popmtcentro(checkpointfile, mt, centro, &t);
  #endif

  initparams(params, mt);                   // initialize MT parameters
  initrodparams(params, rparams, rod, dt);  // initialize rod parameters

  // the main loop

  for (; t<=iter; t++)
  {
    ///////// printing /////////
    if (t%iprint == 0)
    {
      /*
      // Bcast N of each MT
      for(n=0; n<Nmt; n++)
      MPI_Bcast (&mt[n].N, 1, MPI_INT, (n%size), MPI_COMM_WORLD);
      */
      // Bcast positions

      for(n=0; n<Nmt; n++)
      {
        for(k=0; k<Nmax; k++)
        {
          ps = &mt[n].seg[k];
          buffer[(n*Nmax+k)*7+0]= ps->rx; buffer[(n*Nmax+k)*7+1]= ps->ry; buffer[(n*Nmax+k)*7+2]= ps->rz;
          buffer[(n*Nmax+k)*7+3]= ps->q0; buffer[(n*Nmax+k)*7+4]= ps->qx; buffer[(n*Nmax+k)*7+5]= ps->qy; buffer[(n*Nmax+k)*7+6]= ps->qz;
        }
      }

      for(n=0; n<Nmt; n++)
      {
        ptr2buf = &buffer[0] + n*Nmax*7;
        MPI_Bcast (ptr2buf, Nmax*7, MPIreal, (n%size), MPI_COMM_WORLD);
      }

      for(n=0; n<Nmt; n++)
      {
        for(k=0; k<Nmax; k++)
        {
          ps = &mt[n].seg[k];
          ps->rx=buffer[(n*Nmax+k)*7+0]; ps->ry=buffer[(n*Nmax+k)*7+1]; ps->rz=buffer[(n*Nmax+k)*7+2];
          ps->q0=buffer[(n*Nmax+k)*7+3]; ps->qx=buffer[(n*Nmax+k)*7+4]; ps->qy=buffer[(n*Nmax+k)*7+5]; ps->qz=buffer[(n*Nmax+k)*7+6];
        }
      }

      if (rank == 0) print(mt, centro, Nmt, t, dt);

      if (t%iprint_checkpoint == 0)
      checkpoint(checkpointfile, mt, centro, Nmt, t, dt);

    }

    polymerize(mt, params, dt);

    for(n=rank; n<Nmt; n=n+size)
    {
      interfacein(&mt[n], &rod[n], centro);
      rodflex(rparams, &rod[n], t);
      interfaceout(&mt[n], &rod[n]);
    }

    ////// move centro ///////
    if (t%centrofactor == 0)
    {
      for(n=rank; n<Nmt; n=n+size)
      {
        ps = &mt[n].seg[0];
        bufsend[n*3+0]= ps->rx; bufsend[n*3+1]= ps->ry; bufsend[n*3+2]= ps->rz;
      }

      MPI_Allreduce(bufsend, bufrecv, 3*Nmt, MPIreal, MPI_SUM, MPI_COMM_WORLD);

      for(n=0; n<Nmt; n++)
      {
        ps = &mt[n].seg[0];
        ps->rx=bufrecv[n*3+0]; ps->ry=bufrecv[n*3+1]; ps->rz=bufrecv[n*3+2];
      }
      movecentro(params, centro, mt, centrofactor*dt);
    }
  }

  tstop = MPI_Wtime();

  fprintf (stdout, "Elapsed time on proc %3d: %.6e \n", rank+1, tstop-tstart);
  fflush (stdout);

  MPI_Finalize( );

  // free all memory
  gsl_rng_free(rng);
  free(buffer);
  free(bufsend);
  free(bufrecv);

  for (n=0; n<Nmt; n++)
  {
    free(mt[n].seg);
    free(rod[n].elems);
    free(rod[n].nodes);
  }

  free(mt);
  free(rod);
  free(params);
  free(rparams);
  free(centro);

  return 0;
}
