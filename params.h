#include "switches.h"

struct params {
  real YI;                       // bending stiffness of the MTs
  real Vpoly, Vdepoly;           // speeds of polymerization and depolymerization
  real k_catast, k_recov;        // rates of catastrophe and recovery
  real Vmax, Fmax;               // motor max velocity and max force
  real koff, kappa;              // motor koff and stiffness
  real rho_motor;                // number of motors per unit length
  real F0, fric_para, fric_perp; // motor force/length in tangent direction, fric in tangent and normal directions
  real k;                        // spring constant between centrosome and MTs
  real fric_centro;              // friction on centrosome
  real F3stall;                  // stall force for MT polymerization
};
