#ifndef _Scap
#define _Scap

class Chn;
class Res;
class Atm;
class LatSurfv;
class Scap
{
public:

Scap();
~Scap();
void printeznewhelp();
void setcbetaflg();
void printnewhelp();
Atm *findfaratom(Res *);
Atm *createnewatm(Res *,char *);
Atm *createnewatmold(Res *,char *);
Res *findminimaltorsionmore(Res *);
Res *findminimalbordermore(Res *);
Res *findminimalmore(Res *);
void  calcbeta();
void  calcbeta(Res *);
float calcdelphi(Res *r);
float calcularea(Chn *c);
float calcularea(Chn *c,float);
float colonycoeff(Res *);
float colengcoeff(Res *);
float colengcoeff(Res *,float *);
float averageeng(Res *);
float engrmsd(Res *);
float engrmsd(Res *,float *);
void doincludeitself(Res *);
void initialside(int);
void adddisulfide();
void addbackboneh();
void hookside();
void setlist(char,int,char,char);
void readlist(char *);
float minimize(Res *,int *,int *,int);
float minimize();
float iterate(int);
float scpred(int);
float scpred(char *,char *,int); //side chain prediction
float mutate(); //side chain prediction
void remove(Res*,int,float);
void removefast(Res*,int,float);
void removelayermode(Res*,int,float);
void sidelib(char *filnam);
Res  *crtmore(Res *,int nth);  //insert nth more
int  crtmore(Res *,int,float);
void  bbdep(char *filnam);//read bbdep data
void  bbind(char *filnam);//read indep data 
void  writedelphiparm();
float clash(Lattice *,Res *,int,int); 
float clash(Lattice *,Atm *); //clash with lattice
float clashfast(Lattice *,Atm *); //clash with lattice

void insert(Res *);
void averagermsd();
void write();
void printhelp();
void rmsd(Chn *,float,FILE *);
void rmsd(Chn *,float,FILE *,int);
void writeburylist(char *,float);
int *combination(int);
void setsolventorder();
void setsolventorderanyway();
//void setsolventburyanyway();
//void setsolventbury();

//member data

Pdb     *pdb;
Pdb     *pdborg;
int      includeself; //include itself as one rotamer and set intial as original
int      initside;
float ***rotamer;
int     *rotnum;
int      flag;
float    cut;
float   *ener;
float   *armsd;
int      seqn;
int      torsionring;
int      ring;
char     force[10];//1 torsional,2 charge,3 hydrogen,4 all hydrogen;
float    tormax;
float    bmax;
float    mine;
float    cutbb;
int      arbt;
int      nummore;
int      rott;  //rot=1, angle rotamer; rot=0, coordinate rotamer;
int      cbeta;
int      nout;
int      colony;     // colony=1, only colony used in sidechain minimization, colony=2, reassembled in order
int      bordertorsion;
int      singletorsion;
float    nncolony; //colony for backbone only
float    ncolony;  //cofficient of colony energy upon rmsd
int      colonyline; //mutiple initial conformations minimized, the energy compared is colony energy
int      resultout;  //all results out
int      disulfide;  //calculate disulfide bond or not
int      coleng;     //new colony energy formula which takes account energy distribution 
float    coleffect;  //decides the importance of colony energy
int      solventeff; //randomize the iteration order according to the solvent burial percent
Res    **solventorder;   //put all the residues in the array, iterate them sequentially
float   *solventbury;     //solvent burial
int      solventcolony; //the colony effects depends on the burial
float    iterateweight; //the weight of iteraction based on burial
int      iweight;
int      outinitial;
//float    *solventdirectbury;
float    dielectric;
float    dsqrt;
int      outenforce; //increase vdw energy for the surface
int      outforcemag;
int      outcharge;
int      solventorderrandom;
int      printouteng;
int      hydrophobic;
int      adaptvdw;
float    probe;
int	 emilcharge;
int      emilnewcharge;
int 	 emilhydro;
int      delphi;
char    *delfile;
int      delphiself;
float    delphiscale;
int      delphilinit;
int      delphigsize;
float    indi;
float    exdi;
float    hbondeng;
float    hydroeng;
int      hydroflg;
int      mutationonly;
int      allcharges;
int      niterate;
int     *cbetaflg;
int      proring;
int      fastmode;
int      fastmodeH;
int      layermode;
float    prbrad;
float    perfil;
//float    ermsd;
int delphionself;
LatSurfv *latsurfv;
Lattice  *lattice0;
Lattice  *lattice;
Lattice  *lattice1;
protected:
};

#endif
