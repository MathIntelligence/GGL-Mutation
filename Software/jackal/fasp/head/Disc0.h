#ifndef _Disc
#define _Disc

class Disc
{
public:
Disc();
~Disc();
void setupall(Pdb *,float);
void setrange();
void setdipoleascharge();
void setdipoleascharge(Dipole *c);
void setdipoleascharge(Res *,int);
void ready();
float clash(Res *r,int n);
float clash(Lattice *,Res *r,int n);
float clash(Lattice *,Atm *);
void write(FILE *);
Dipole *finddipole(Atm *a);
Charge *findcharge(Atm *a);
void setup();
void setup(Res *,int n);
void calcenglist(FILE *fp,char *file,char *org);
Dipole  **dipole;
Charge  **charge;
Dipole  **range;
char      force[100];
float     energy[200]; 
float cutoff;
int logg;
Pdb *pdb;
};

#endif
