#ifndef _Chn
#define _Chn

class Pdb;
class Res;
class Atm;
class Tres;
class Lattice;
class Chn{
public:

Chn();
~Chn();
Chn(Chn *);
float calcareadiff(int);
void writerescard(FILE *);
Res *findsmallres(int);
int isdatabaselinked(Res*,Res*);
int checkpdb(int);
void headerpdbfix(Res *s,int n);
void buildcatracehbond();
void increaseid(Res *,int,int);
void writeoldmore(FILE *s);
void writessbondlist(FILE *fp);
void clearhbond();
//member function
void setdsspstr(int n);
void dihedraloption(FILE *fp,int m);
void dihedral(Res *,int);
int  calcnumber(Res *,int);
int  isalllinked(int,int);
int  isalllinked(Res *,int);
Res *findmaxclash(Res *,int);
float *gettransfer(Res *,int);
void addmissingatms(Res *,int);
void setresenergy(Res *,int);
void setresenergy(Res *,int,int);
float getresenergy(Res *,int,int);
int islinked(Res *,Res *);
void printenergy(FILE *fp,Res *,int);
void setbackbonetorsion(Res *,int);
void setfreenextonmore();
void buildcontact(Lattice *lat,Res *r0,int n,float dcut);
void deletecontact();
void setlastresorder();
void replace(char *,int,int);
void insertresidues(char *,int at,int direct);
void writeoutrotamer(FILE *fpp,char c);
void writesecondary(FILE *fp);
int isssbond(Res *,Res *);
void setallnear(Res *,int);
void allnearbond(int);
void setreversezero();
void setreversethread();
Res *lastres();
Res *isres0(int);
Res *isresoid(int);
Atm *getatmbyoid(int);
float getarea();
int  isnearnext(Res *,Res *);
void setthreestatesec();
void setoid(int);
void dihedral();
void dihedral(char *);
void setdsspstr();
void findturnstructure();
void findbetasheets();
void testladder(Res *,Res *,Res *);
void findalphahelix(int);
void setsecstr();
void addhatoms(int f);
void addhatoms();
void deletehbondlist();
void buildhbond(int,int);
void buildhbond();
void writehbondlist(FILE *);
void writehbondparm();
void buildssbond(int,int);
void buildssbond();
void setseqcardnogap();
void setseqcard();
void transform(int,int);
void addatmid0(int);
void addsidechain();
void addresid0(int);
void addresid(int);
void addresoid(int);
char *getseqn(Res *,int);
char *getseqn();
void initial();
Chn *clonechain();
void addresiduesonly(int from,int to,char *seqn);
float *center(Res *,int,int,int);
void configure();
void setflg(int);
void setflgr(int);
void setflgr(char c,int f);
void write(char *,Res *,int,int);
void write(char *,Res *,int);
void write(char *);
void write(FILE *fp,Res *,int);
void write(FILE *fp,char);
void writeold(char *);
void writeold(FILE *);
void write(FILE *fp);
void dihedral(FILE *fp);
void dihedral(FILE *fp,int);
void dihedral(FILE *fp, char);
void dihedral(FILE *fp, char,int);
void setoff(Res *,int,float);
void setoff(float *,float);
void transfer(Res *,Res *,int);
void transfer(float *,Res *,int,int);
float *gettransfertemp(Res *,int);
void transfertemp(float *,Res *,int,int);
void transfer(float *,int);
void transfer(Res *,float *,int);
void transfer(Chn *s);
void surfaceold(int,float);
void surface(int,float);
void surface(char *,int,float);
void onlybc();
void secstr(FILE *fp); //create the secondary structure
void create(char *);
void create(Tres *);   //create the chain of this twenty aa
void create(Chn *s);   //create the chain of sequence of chain s
void create(Res *,int); //create the chain of sequence of chain s segment 
void transfer(float);  //transfer coordianates
void copytemp(float *,int);
void transform(Res *,int,int);
void transform(int );//transform the chain into different format;
void header(int,int);
void header();
void header(Res *,int);
char *getseq();
char *get2d();
int getseq(char *);    //get the sequence of this chain
int ssbond();          //calculate the s-s bond
int manyatm();
int manyres();
int lastid(int);
void seqout(FILE *fp);
int  setflg(int *,Res *,int,int,int);
void setflg(Res *,int,int);
void setflg(int *,Res *,int,int);
void setflgr(Res *,int,int);
void setflgr(int *,Res *,int,int);
int  good(int begin,int end,int m); 
Res *operator[](int resn); 
Res *isres(int);       //the sequential residue 
Res *isres(char, int);
Atm *isatm(int);
void bondleng();
void idequal(Res *r1,Res *r2,int n);
void buryangle(char *,int,int);
void buryangle(FILE *,int,int);
void setchnonly();
//data

Res *res;
int  ishet;
char id;
int number; //number of residue
Chn *next;
Pdb *pdb;
float area;
int start;
int   nemp;
float energy;
float *temp;
char *seqcard;
Chn *more;
};
#endif
