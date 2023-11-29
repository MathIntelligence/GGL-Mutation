#ifndef _Disc
#define _Disc

class Disc
{
public:
Disc();
~Disc();
void clear();
void setupall(Pdb *,float,int);
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
void setup(int);
void setup(Res *,int n,int);
void calcenglist(FILE *fp,char *file,char *org);
int       smoothclash;
Dipole  **dipole;
Charge  **charge;
Dipole  **range;
float  *change;
char      force[100];
float     energy[200]; 
float cutoff;
float eps;
Pdb  *pdb;
};

#endif
