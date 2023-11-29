#ifndef _Pdb
#define _Pdb

// this class(plus chn,res,atm) is to store protein information

class Aln;
class Chn;
class Atom;
class Pdb{
public:

Pdb();
Pdb(Pdb*);
~Pdb();

//member function
void writepdbused(char *);
void writepdbused(FILE *fp);
void cleanemptychn();
void setrescard(int);
int  checkrescard();
void adaptindex(Res *r,int n);
int isconnected(Res *a,Res *b);
void indexreswithcard();
void indexreswithcard(Chn *c);
//read pdb file 
int ishetatmchainexist(Chn* c);
int  isresidchar(char *);
void checkhetatmchainid();
void relinkhetatmchain(Chn *,Chn *);
void removehetatmchain(Chn *);

int checkpdb(int);
void parsehetatmline(char *);
void addmorechain(Chn *,Chn *);
void relinkchain(Chn *,Chn *);
void setallatmid0();
void checkpdbchainid();
void addxyz(float *);
void getmaxxyz(float *);
void headerpdbfix();
void setseqcard();
void writehbondlist(FILE *);
void setseqcardnogap();
void hetatmconfigure();
void configuremoreatmid0();
void writeseq(FILE *);
void writeseq(char *);
void writeoldmore(char *);
void setatmoid();
char *getname(char *);
char *getname();
Atom **getAtomAll();
void writesurface(char *);
void writesurface(FILE *);
void setresid0back();
void setatmid0back();
void writehbondlist(char *s);
void surface(int,float);
void setallnear();
void setthreestatesec();
void deletecontact();
Atm *getatmbyid0(int);
Atm *getatmbyoid(int);
float getarea();
void writeoutrotamer(FILE *fpp,char c);
void transform(int);
void transform(int,int);
void removechain(Chn *);
void header();
int  manyres();
int  maxresid0();
int  maxatmid0();
int  ischainexist(Chn *);
Res *isres(int);
void readseqcard();
void readseqcard(char);
void giveresmoreid();
void transfer(float *,int);
void setlastresorder();
void setoid(int);
Pdb *clonepdb();
int  manyatm();
int  getresnum();
int  manychain();
void center();
void allresid0();
void addhatoms(int);
void addsidechain();
void setflg(int f); 
void dihedraloption(FILE *fp,int);
void dihedraloption(char *,int);
void dihedral(FILE *fp);
void initial();
void configureatmid0();
void configure();
void allresid();
void setoff(float *xyz,float cut);
void setoff(Res *r0,int n,float cut);
int getidnumber(Chn *);
int read(char *filnam,char ch); //ch: '0' all; '1' the first one; 
int readmore(char *filnam);
void convrt(Chn *,FILE *); //convert the pdb file to charmm format
void writedelphi(char *,Res *r,int,int,int);
//write pdb file
void writerescard(char *);
void writeold(FILE *);
void writeold(char *);
void write(char *);
void write(FILE *fp);
int  dssp0(char *filnam);
void dssp(char *filnam);
void setflgr(char,int);
void setflgr(int);
//search for chain
Chn *lastchain();
Chn *ischain(int ch);
Chn *ischain(char ch);
Chn *operator[](int ch); 
Chn *operator[](char ch);
//member data
char *rotamername;
char **info;
char **endinfo;
char *name;
char *token;
int number; //number of chain
Pdb *next; // x-ray file
Chn *chn;
Chn *hetatm;
float area;
int isinfo;
int delchn;
int ishet;
};

#endif
