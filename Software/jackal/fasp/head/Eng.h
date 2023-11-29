#ifndef _Eng
#define _Eng

class Tatm;
class Eng
{
public: 
Eng(Tatm *);
~Eng();
void tank();//create dimensions for parameters
void charm(char *,char *);//read charmm parameter
void amber(char *,char *);//read amber parameter
Tatm *tatm;

//bonds ---   V(bond) = Kb(b - b0)**2
float *bond0;
float *kbond;
int mbond;
//angles --   V(angle) = Ktheta(Theta - Theta0)**2
float *angle0;
float *kangle;
int mangle;
//dihedral -- V(dihedral) = Kchi(1 + cos(n(chi) - delta))
float *chi0;
float *kchi;
float *nchi;
float *chi02;
float *kchi2; 
float *nchi2;
int mchi;
//nonbonded--V(Lennard-Jones) = Eps[(Rmin/r)**12 -2(Rmin/r)**6] 
float dipole;
float charge;
float radius;
float epslon;
float charge0;
};
#endif
