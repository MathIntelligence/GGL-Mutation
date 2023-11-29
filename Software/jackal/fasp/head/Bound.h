#ifndef _Bound
#define _Bound

class ModTop;
class HashTable;
class DistSecond;
class Bound
{
public:
Bound();
~Bound();
int  setpdbfixbound(int flg);
void getnewsequential(Atm **,Res *);
void setpdbfixbound();
void setatmat(Chn *);
void setsmallbound();
void easyready();
void setbond();
void setbondold();
void delsidechainbound();
void setcurve(float);
int myoptimize(int);
int bondoptimize(int,int);
int smoothoptimize(int);
int calcpredt(Res *r);
int  getnumberofpredt();
int  steepoptimizeold(int);
int  steepoptimize(int);
void predtorder();
void setrange();
void setpredt(int,int);
void setpredt(Res *,int,int);
void setpredt(int);
void setallbackbone();
void setdssp();
void detectbadbound(FILE *fp);
void detectbadtriangle();
void writeatmbound(FILE *fp,int nn);
void setmodel(Pdb *);
int setalldist();
int setfmtcloseresidue();
int setcloseresidue();
int setclosedistance();
int findpair(HashTable **,int,int);
int isnextres(Res *,Res *);
void setstartsize(int n);
void defaultaction();
void boundsmooth(int);
void adjust();
void searchstructure(int);
void initialstructure();
void printoutdistance();
void printoutimproper();
void settag();
void setbuffertag();
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
void clearredundancy(); 
void addbounds(Atm *,Atm *,float,float);
void addbounds(int,int,float,float,float,float);
void setsize(int);
int readdist(char *filename);
int  optimize(int);
int  optimize();
void setoptimize();
void createnewbound();

//not important
int selectoptimize(int,int);
void setbackboneaction();
void setbuddy();
void getsequential(Atm **,Res *);
void getring(Atm **,Res *);
int setcloseresidueold();
int setsecondarybounds();
int setnearnextbounds();
int  setinternaldist(); //1-2,1-3,ring
int quasiringdist(Res *r);
int setphipsidist();
int set14dist();
int set14dist(Res *r);
int setn14dist();
int setn14dist(Res *r0,Res *r);
int setsidechainrestriction();
int setsidechainrestriction(Res *r);
int setquasiringdist();
int sethbondbounds();
int setexternalbounds();
//data member
//atom list
Atm       **atmat;
Atm       **atoms;     //all atom list
int        *predt;     //flag for all atoms whose position to be predicted, default is 1, to be predicted
int        *tag;       //the atm list start and stop position, the first atoms
int         natom;     //number of atoms
Atm       **improper;  //improper dihedral list
int         nimp;
float      *weight;
float      *dist;      //constrained distance
float      *range;
float      *curve;     //a(r-b)^2+c, only store a.b is the mid point.c is the well of lowest energy
int        *atmpair;   //atom pair id
int        *bond;
int         ndist;     //number of constraint;
int         arraysize; //totol array size;
Pdb        *model;
int         seed;
int         flag;
int         chiropt;
int        *dssp;    
int	    maxnit;
float       npower;     
ModTop     *modtop;
DistSecond *distsecond;
Bound      *next;  
};
#endif
