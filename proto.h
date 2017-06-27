// Function prototypes
#include "struct.h"
#include "rod_struct.h"

void initparams(struct params *params, struct mt *mt);
void initmtcentro(struct mt *mt, struct centro *centro);
void popmtcentro(char *checkpointfile, struct mt *mt, struct centro *centro, int *t);
void initrodparams(struct params *params, struct rodparams *rparams, struct rod *rod, real dt);
void polymerize(struct mt *mt, struct params *params, real dt);
void randomorient(struct mt *mt);
void randomorientcone(struct mt *mt);
void addseg(struct mt *mt, struct params *params);
void delseg(struct mt *mt, struct params *params);
void rotmatrix(struct seg *ps);
int  checkboundary(real rx, real ry, real rz, real dia, real thick);
void checkboundarysquare(int *proximity, real rx, real ry, real rz, real side, real thick);

void interfacein (struct mt *mt, struct rod *rod, struct centro *centro);
void rodflex(struct rodparams *rparams, struct rod *rod, int t);
void interfaceout(struct mt *mt, struct rod *rod);

void tempinterfacein (struct mt *mt, struct rod *rod, struct centro *centro, struct rodparams *rparams);
void tempinterfaceout(struct mt *mt, struct rod *rod);

void eulerupdate(struct mt *mt, real dt);
void movecentro(struct params *params, struct centro *centro, struct mt *mt, real dt);

void makedatadir(char **argv);
void opendatadir(char **argv);

void print(struct mt *mt, struct centro *centro, int Nmt, int t, real dt);
void checkpoint(char *checkpointfile, struct mt *mt, struct centro *centro, int Nmt, int t, real dt);
void Qinv(real Q[][8]);
