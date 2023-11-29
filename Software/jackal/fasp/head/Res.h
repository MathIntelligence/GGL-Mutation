#ifndef _Res
#define _Res

class Chn;
class Atm;
class Tres;
class HBondList;
class Res{

public:

Res();
Res(Res *);
Res(Tres *);
~Res();

//member function
float directrmsdanyway(Res *,int,int);
int  sethnatomxyz(Atm *);
int  issidechaincomplete();
int  isbackbonecomplete();
void addhnatoms(int);
void addmissingatms();
void writeoldmore(FILE *fp);
void dihedraloption(FILE *fp,int m);
void deletecontact();
void setbackbonetorsion();
void hooksidechain();
void transferbackbone(Res *);
void setfreenextonmore();
void addambgt(char *);
void writessbondlist(FILE *);
void writehbondlist(FILE *);
void addhatoms(int);
void setoid(int);
int ishbondedsource(Res *);
int ishbondedwithatm(Atm *,Atm *);
int ishbondedwithres(Res *);
int ishbondedwithres(int);
int setsidechainhatomxyz(Atm*);
int setbackbonehatomxyz(Atm*);
void addbackbonehatoms(int);
void addsidechainhatoms(int);
float colonycoeff();
float directrmsd(Res *,int,int);
void mutateResidue(char);
void addsidechain();
void configure();
void configureaddhatom();
Res *cloneres();
void initial();
void write(char *);
void writerescard(FILE *fp);
void writeold(FILE *fp);
void write(FILE *fp,int a,int b); //write atoms between a and b
void write(FILE *fp);     //write out this residue
void write(FILE *fp,int); //write out residue in different format
void dihedral();
void dihedral(FILE *fp);  //write out dihedral information
void dihedral(FILE *fp,int);
int  ssbond(Res *r1,Res *r2,float s_s,float s_cb);
float clashnoring(int,int);
float clash(int f);       // clash inside the residue
float clash(Atm *,int);   //residue clash with this atom
float *center(int ,int);
void setflg(int f);
void transfer(Res *,int);
void transfer(Res *);
void transfer(Atm *,float *,int f);
void transfer(float *,int f);//f=1 to res;
void transfertemp(Res *);
void transfertemp(float*,int);
void transfer(float); 
void copytemp(float *,int);
void copytemp(Res *); // copy temp
void copytemp(float *);
void giveid();  //give id to more residues
void bondleng();
//void givedad();//assign dad to more residues
int  many();//number of more residues;
Res *findbestmore(float,float);
Res *ismore(int); //get the more residue of id n
Atm *operator[](char *);
Atm *operator[](int);// get the atom from nth atom
Atm *isatm(char *);
Atm *isatmid(int);
Atm *isatm(int);
float bury(int,int);
Res *isres(int);
void attachAtm(Atm *);
void takeoffatm(Atm *);
int  getuse(char *,int,int);
int  getuseambgt(char *,int,int);
int  getid();
int  rotable(Atm **,int,int);
int  rotable(int,int);
int  good(int begin,int end,int m);
void addhbondlist(Atm *,Atm *,float,int);
void buryangle(FILE *fp,int,int);
float rmsd(Res *,int,int);
float ambgtrmsd(int);
float getarea();
//member data

Chn *chn;
Res *last;
Res *next;
Res *more; // the other conformation of this residue
Atm *atm;
Tres *tres;
HBondList *hbond;
char name;   
int id;      // reside id with sequence gap
int id0;     // sequential number of reside id 
int oid;
int number;  // number of atoms in the Res
int nhydr;   // number of h atoms
int nummore;
int regular; // regular residue
float  area;//solvent accessible area
char   sec; //secondary structure
Atm   *ambgt; //ambiguity choice for residues;
float *temp; //temparory buffer
int    nemp; //number of this buffer
float  energy;
float  flag; 
float  rotchance;
float  rotenergy;
int    rotdegree;
int    addhanyway;
char   rescard[6];
};
#endif

