#define real double

struct rodparams {
  real dt;
  // Parameters inherited from the main code
  real YI;                                  
  real F0;                                   // motor force per unit length
  real k;
  real le;
  real zt1, zt2, zt3, zr1, zr2, zr3;         //Frictions per unit length

  // Independent parameters
  real dia, rho, sigma;                      // diameter, volume density, Poisson's ratio - Independent parameters

  // Derived parameters
  real A, I, Y, G;                           // Area, area moment, Young's modulus, shear modulus
  real i1, i2, i3, xm;               
  real cs1, cs2, cs3, cb1, cb2, cb3;         // Modulii
  real s01, s02, s03, b01, b02, b03;         // Reference strains

  // Cell size
  real d, thick;
};
