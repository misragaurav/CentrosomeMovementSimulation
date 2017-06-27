#include "rod_params.h"

struct elem {
  real q0, qx, qy, qz, q, q2;
  real d[3][3], e[3][4];
  real f1, f2, f3, t1, t2, t3;             // Force and couple
  real pe;                                 // Pot Energy
  real peT1, peT2, peT3, peR1, peR2, peR3; // pot energies of the elements
}; 					     
					     
struct node {
  real q0, qx, qy, qz;		           // Quaternions		     
  real d[3][3], e[3][4];		   
  real rx, ry, rz;                         // Coordinates
  real px, py, pz, lx, ly, lz;             // Momenta in Space frame 
  real p1, p2, p3, l1, l2, l3;             // Momenta in Body frame
  real f1, f2, f3, t1, t2, t3;             // Force and couple
  real fx, fy, fz, t0, tx, ty, tz;         // Forces in SFF
  real zt1, zt2, zt3;                      // translational friction
  real ke; 				     
  real keT1, keT2, keT3, keR1, keR2, keR3; // energies of the DoF */
};

struct rod {
  struct node *nodes;                        // Nodal data structure
  struct elem *elems;                        // Element data structure
  int  n_nodes, n_elems;                     // Number of nodes, elements
  int pinned;
  real F3tip;
  real dL;
  real centrorx, centrory, centrorz;         // centrosome position
  real pe, ke;
  real keT1, keT2, keT3, keR1, keR2, keR3;     
  real peT1, peT2, peT3, peR1, peR2, peR3;     
  real xcm, ycm, zcm;			     // Center of mass coordinates
  real pxcm, pycm, pzcm, lxcm, lycm, lzcm;   // Center of mass momenta
};
