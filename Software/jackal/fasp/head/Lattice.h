#ifndef _Lattice
#define _Lattice

class Pdb; 
class Cell;
class Lattice
{
public:
Lattice();
~Lattice();
void putonhetatm(Pdb *);
void initial();
int  indx(float *s);
int  indx(float s);
void clean();
void printoutLattice();
void addifnotexist(Atm *);
void ready(Pdb *s);
void putonhetatm();
void puton();
void puton(Chn *,int,int,int);
void puton(Pdb *,int,int,int);
void puton(Res *,int,int,int);
void puton(Chn *,int begin,int end);
void puton(Pdb *,int begin,int end);
void puton(Res *,int begin,int end);
void puton(Res *,int);
void puton(Pdb *);
void puton(Chn *);
void puton(Res *);
void puton(Atm *);
int  putoff(Chn *,int,int,int);
int  putoff(Pdb *,int,int,int);
int  putoff(Res *,int,int,int);
int  putoff(Chn *,int begin,int end);
int  putoff(Pdb *,int begin,int end);
int  putoff(Res *,int begin,int end);
int  putoff(Res *,int);
int  putoff(Res *);
int  putoff(Atm *);
int  putoff();
int  getcell(float,float,float,float);
int  getcell(float *,float);
int  getcell(Atm *,float);
int  cutoff(float *,float);
int  cutoff(Atm *,float);
int  cutoff(float *);
int  cutoff(Atm *);
int  resonly();
int  getrescutcharge(Atm *);
int  getrescut(Atm *);
int  cover(float *,float);
int  exist(Atm *);
Pdb   *pdb;
Cell **busket;
Cell **atom;
Cell **obtain;
int    hash; 
int    nget;
float  grdsiz;
float  radall; //maximum of include all
float  radmax; //maximum of radius
int  flag; //0 not include this residue
int  hydrogen; //1 including hydrogen, 0 no hydrogen
int  offset;   //offset for atom sequential id
int  ishet;
//Atm *atm;
};
#endif
