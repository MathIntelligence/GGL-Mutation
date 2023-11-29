#ifndef _DistBound
#define _DistBound

class ModPdb;
class HashTable;
class DistSecond;
class DistBound
{
public:
DistBound();
~DistBound();
void setmodel(Pdb *);
int findpair(HashTable **,int,int);
int isnextres(Res *,Res *);
int quasiringdist(Res *r);
int modpdbquasiringdist(Res *r);
int setphipsidist();
int setmodpdbphipsidist();
int set14dist();
int setmodpdb14dist();
int set14dist(Res *r);
int setmodpdb14dist(Res *r);
int setn14dist();
int setmodpdbn14dist();
int setn14dist(Res *r0,Res *r);
int setmodpdbn14dist(Res *r0,Res *r);
int setsidechainrestriction();
int setmodpdbsidechainrestriction();
int setsidechainrestriction(Res *r);
int setmodpdbsidechainrestriction(Res *r);
int setquasiringdist();
int setmodpdbquasiringdist();
int setmodpdbdistance();
void setstartsize(int n);
void defaultaction();
void boundsmooth(int);
void adjust();
void searchstructure(int);
void initialstructure();
void printoutimproper();
void settag();
void boundsquare(float);
float evaluate();
void writeoutdistbound(char *);
void writeoutdistbound(char *,int,int);
void mirror(Atm *,Atm *,Atm *,Atm *);
void reorder();
void deletefixedbound(int);
void deletebound(int);
int  atmexist(Atm*);
void swapatm(Atm *,Atm *);
void swapatm(int,int *);
float getta(Atm *,Atm *,Atm *,Atm *);
void setimproper();
void ready();
void addnewatm(Atm *s,int);
void getsequential(Atm **,Res *);
void getring(Atm **,Res *);
void clearredundancy(); 
void addbounds(Atm *,Atm *,float,float);
void addbounds(int,int,float,float,float,float);
int  setinternaldist(); //1-2,1-3,ring
int  setmodpdbinternaldist();
void setsize(int);
int readdist(char *filename);
int  optimize(int);
int  optimize();
void createnewbound();

//data member
//atom list
Atm  **atoms;     //all atom list
int   *predt;	  //flag for all atoms whose position to be predicted, default is 1, to be predicted
int   *tag;       //the atm list start and stop position, the first atoms
int    natom;     //number of atoms
//improper dihedral list
Atm   **improper;
int    nimp;
//dist constraint list
float *weight;
float *dist;      //constrained distance
int   *atmpair;   //atom pair id
int    ndist;	  //number of constraint;
int    arraysize; //totol array size;
Pdb   *model;
int    seed;
int    flag;
int    chiropt;
ModPdb *modpdb;
DistSecond *distsecond;
DistBound *next;  
};
#endif
