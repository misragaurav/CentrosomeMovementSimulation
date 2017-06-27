// MT data structures
#include "params.h"

struct seg {				     
  real rx, ry, rz;             // Coordinates
  real q0, qx, qy, qz;         // quaternions
  real R[3][3];                // rotation matrix - space to body
  real zt1, zt2, zt3;          // translational friction
};

struct mt {
  struct seg *seg;           // Nodal data structure
  int id;                    // MT identification number
  int N, status;             // Number of segments
  int pinned;
  real dL;
  real fcx, fcy, fcz;        // force between mt and centro in space frame
  real F3tip;
};

struct centro {
  real rx, ry, rz;
  real fx, fy, fz;
};

// Global variables
gsl_rng * rng, * MPI_rng;

real ds;   // segment size
int Nmt;   // number of microtubules
int Gmax; // Maximum segments globally
