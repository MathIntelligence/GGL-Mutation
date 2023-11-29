#ifndef _DistMutateBound
#define _DistMutateBound

class ModPdb;
class HashTable;
class DistSecond;
class DistMutateBound
{
public:
DistMutateBound();
~DistMutateBound();
int setcloseresidueold();
int setsecondarybounds();
void setallbackbone();
void setdssp();
void detectbadbound(FILE *fp);
void detectbadtriangle();
void writeatmbound(FILE *fp,int nn);
void setmodel(Pdb *);
int sethbondbounds();
int setexternalbounds();
int setalldist();
int setcloseresidue();
int setclosedistance();
int setnearnextbounds();
int findpair(HashTable **,int,int);
int isnextres(Res *,Res *);
int quasiringdist(Res *r);
int setphipsidist();
int set14dist();
int set14dist(Res *r);
int setn14dist();
int setn14dist(Res *r0,Res *r);
int setsidechainrestriction();
int setsidechainrestriction(Res *r);
int setquasiringdist();
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
void writeatmbound(char *);
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
Atm   **improper; //improper dihedral list
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
int   *dssp;    
ModPdb *modpdb;
DistSecond *distsecond;
DistMutateBound *next;  
};
#endif
