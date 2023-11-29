#ifndef _Tres
#define _Tres

class EngCoeff;
class SimAtm;
class Tatm;
class Rotamer;
class DistPopularBin;
class ModConstant;
class EngTable;
class Tres
{
public:
Tres();
~Tres();
void   setenghetatm(char *);
Tres  *addhetm(char *);
Tatm  *addhetmatm(char *);
Tatm  *updatehetmatm(char *);
Tres  *updatehetm(char *);
Tatm **findtwofarringatm(Tatm *);
Tatm *findfarringatm(Tatm *a);
float *getdipolevector(Atm *,int);
void  switchcharge(int);
void  printchnparserhelp();
void  printxyzseqhelp();
void  printezchixyzhelp();
void  printchixyzhelp();
void  printseqxyzhelp();
void  printnalgnhelp();
void  setdipolecharge();
void  setsimplecharge();
void  calcbackrotamertorsion(char *);
void printalgnhelp();
void printhbondhelp();
void printchihelp();
void  printaddhhelp();
void  printsurfacehelp();
char *getcappedsequence(char);
char *getcappedsequence(char*);
char *getsequence();
ConTable *findcontable(char *);
EngTable *findengtable(char *);
void prepareengtable();
void setcontable();
void createengtable(int);
void writesimatm(FILE *);
void setsimplesimatm(char *);
void  setsimatm();
int   getaanum();
char *delnonstandard(char *);
int  getresid(char );
void set3backrotamer(char *);
Pdb *findrotamer(char *s);
Pdb *findrotamername(char *s);
void setdefaultsimatm();
void readdistpopular(char *,char *);
void readmore(char *);
void readmore(char *back,char *side);
void read(char *);
void setrotm();
int  getorder();
int  isaa(char);
Tatm *findanyatmwithname(char *); 
char *swap(char); // get the three character residue name
char  swap(char *); // get the one character residue name
Tres *operator[](int); // get the nth residue
Tres *operator[](char); // get the residue from its name
Tres *operator[](char *); // get the residue from its name
Tatm *isatm(int n);       // nth atom
Tatm *istatm(char *s);       // nth atom
Tatm *isatm(char *s);     // atom of name s
Tatm *isatm(char *s,int m,int n,int t);//t:the nth atom
void  writechargeradius(FILE *);
void  writedelphicharge(char *);
void  writedelphisize(char *);
void  assigncharge(char *f);
int   ispolaraa();
int   isnonpolaraa();
void  ispolar();
int   getuse(char *s,int m,int n);     //m:ignore;n:extend
int   rotable(Tatm *);
void  vectorm(float *,float *,float *);
float vectord(float *,float *);
void  vectorx(float *,float *,float *);
float distance(Atm *,Atm *);
float distance(float *,float *);
float distance(float *);
float distance(float,float,float,float,float,float);
float distsqr(float *);  
float distsqr(float *,float *);
float distsqr(float,float,float,float,float,float);
float angle(float *x,float *y,float *z); 
float dihedral(float *,float *,float *,float *);//four coodinates
float space(float *,float *,float *);
float maxlen(int,int);
void  maxlen(int);
void  surface(int,float); 
void  surface(int,float,int);
void  copy(float *from, float *to,int ); 
void  transfer(float *,int);
void  smoothf(char *);
void  setdonaldcharge();
void  setchargezero();
void  setnewengcoeff(float,float);
void  setengcoeff(float,float);
void  setengcoeff(int);
//void createvdwtable();
//void createpairdeftable();
//void createengtable();
int  id;
int  flag; //1**, charm,2** amber,**0,all atom,**1,extend,**2, no-hydro
char name;
char name3[4];
int  number;//total number
int  nhydr; //number of H
int  numrot; //total number of rotable
Tres *hetm;
Tres *next;
Tatm *tatm;
float charge;//charges of this residue
float area;  //maximum surface of extended residue
float mlen;  //maximum length of the residue
float *head;
float *smooth;
int    nsmt;
int    smt;
int    logg;
int    smoothclash;
int    smoothmin;
int    selfvdw;
FILE  *logf;
Pdb   *rotamer;       //coordinate rotamer library
ConTable *contable;
SimAtm *simatm;       //18 atom types
EngTable *engtable;   //energy table 
EngCoeff *engcoeff;   //energy coefficient for vdw
DistPopularBin *popbin;
ModConstant *constant;
};
#endif
