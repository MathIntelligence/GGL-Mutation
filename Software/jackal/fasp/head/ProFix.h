#ifndef _ProFix
#define _ProFix

//loop prediction program
class Dace;
class Mutate;
class Bound;
class Pdb;
class Chn;
class Lattice;
class SegBed;
class Disc;
class ProFix
{
public:
           ProFix();
           ProFix(ProFix *);
          ~ProFix();
void setpdbfixbound(int flg);
float** zipperfix0(float dis,int flg);
float **predtpdbfix();
float *	   fixpdbsegment(Res **);
int        iscompositesegment();
float *    scapfix(float **);
float **   cutlargermsd(float **xyz_all,float rcut);
float **   cutlargermsd(float **xyz_all,float rcut,int);      
void       setboundtag(int);
void       clear();
void       updatexyzsave(int*);
int       nolarge(float **,float);
float *	calcenergycoeff(float *ent,int num);
float**       create7d(int);
float**       create6d(int);
float**       create5d(int);
float      getenergy(Res *s,int to);
float     *avgminrmsd(float **,float *);
float      clashold(Res*,char *);
float	   clashold(Atm*,char *);
float     **create4d(int);
float     **create3d(int);
float     **create2d(int);
float     **allavgminimize(float **);
float     **indavgminimize(float **);
void	   clearxyz(float **,int);
Pdb *      mycreate(Res *,int);
float     **avgminimize(float **);
float     **addup(float **,float *);
float     **addup(float **,float **);
float     **newfixsegment(int *);
float     *myfixsegment(int *);
int 	   resort(SegBed *);
void       readyminfix(Res *,int);
float    **getnewxyz(float **,Res *,int);
float    **minfix(float **);
float      minimizefix(float *,int);
float*     boundminimizefix(float **,int);
float*     minimizefix(float **,int);
float    **zipperit(float);
float    **zipperit();
float      predtmin();
float      calcavgrmsd(float **,int );
void       setbound(int,int);
void       setbound(int);
int        ifrotatmexist(Atm *);
float    **gofix();
float     *minimizesidechain(float **,int);
float    **zipperallpossible();
void       setallnear();
void       hooksidechain();
void 	   readynewfix();
void 	   readyfix();
float      boundmin(int);
float **   predtfix();
Res 	  *getmidresidue(Res *,Res *);
int        findmidresidue(int,int);
int 	   findmidresidue(Res *,Res *);
int        randmzChianglefix(int,int);
void       randmzfix(ProFix *,ProFix *);
void	   randmzfix();
void	   segmentfix();
float **   gofix(int);
float *	   getaveragepost(float **);
float	   clashfix(Atm *);
float **   zipperboth(float);
float **   zipperfix(float,int);
float     *ensortfix(float **xyz); //energy sort
float     *ensortold(float **xyz,char *);
int        fixsegment(int *);
int        linkviaboth(int,int,int);
int	   getdirection();
void       initial();
void       ready();  //preparing for loop
int	   randmzChiangle(int);
float      sidechainmin(int);
void       checkomega(FILE *);       
int	   checkomega();
float      closeg(ProFix *,ProFix *); // tweak two segment to consistency
float      closeg(); // tweak segment to constraint
float      *trim(float **,int haa,int);
float    **around(float **);
float    **zipper(float,int);  //random loop generator
float     *expand(float **,int flg); //expanding nearby
float      expand(float *,int flg);  //expanding nearby
float    **expand(float *,int,int,int); //loop searching around neighbor and return all found
float     *peel(float **,int got);
float    **fill(float **,float);
float    **trim(float **);
float     *trim(float **,float *,float d); //trim those large than d
float     *trim(float **,float *,int haa);//trim those no more than part
float     *trimold(float **,int haa); //trim those no more than part
float     *trim(float **,int haa); //trim those no more than part
float    **generateold(float **,int,float);
float    **generate(float **,int,float);
float    **predt0save();
float    **predt0old();
float    **predt0();
float    **predt(int);
float    **predtold(int);
float    **predtsave(int);
float    **copy(float **);
float    **copy(float **,int);
float    **refill(float **,int need);
void       testclosure();
float      segmin(int); //loop minimizer
float      segmin0(float,int); //loop minimizer
int        noclose(float **xyz,float d);  //noclose than d
float      noclose(float **xyz,float *xyz0); //find closest
float      minimize(float *xyz_all,int flg);
float     *minimize(float **xyz_all,int flg); //minimization
float    **cross(float **xyz,int cot,float cut);   //loop cross
float      ensort(float *xyz); 
float     *ensort(float **xyz); //energy sort
void       setoff(float **,float);
int        setend(int); //satisfying end meet
//int        readmore(char *,char *);  
int        create();
void       create(Res *,int);        
void       create(char *seqn);
void       randmz(int rand);
void       randmz();
void       writeseg(float **,char *,int);
void       write(float **,char *,int);
void       printhelp();
float      scpred();
void       crtmore();
float     *scap(float **);
float      scap();
float      rmsdmutate(int);
float      rmsd(int);
float      clash(Res *,int from,int to);
float      clash(Res *,int to);
float      clash(Res *);
float      clash(Atm *);
void      randmzHelix(int);
void 	  randmzSheet(int);
int      **allconf(int);
float    **zipper();

//data member
Pdb       *pdborg;
Pdb       *pdb;
int        start;
int        end;
Pdb       *segment;
Lattice   *lat;
float      step;
float      cutoff;
ProFix     *last;
ProFix     *next;
int        direct; //direction
int        arbt;  //arbitrary number of first attemp
int        id;
float    **near;
float    **init;
float     *xyzsave;
int        part;
int        revs;
int        pdbout;
float      outlet;
int        flat;  //control minimization
int        cycle;
char       secd;
char	   cid;
int	   part0;
float      discut;
int        nomin;
int	   databaseonly;
int        numclose;
Atm      **rotatm;
Bound     *bound;
Mutate    *mutate;
int        ranseed;
int        fapr; 
int        flex; //determine when the linking is not impossible
int        flexrot; //flexible rotation       
int	   onlybound;
int        onlysidechain;
int	   smoothclash;
int        onlyenergy;
//int        charge;
int        rmsdinsert;
int        randcoil;
int	   logg;
int        onlyxyz;
char       boundtag[100];//'H', hydrogen 'P': position from segment
Chiangle  *chiangle;
Disc      *disc;
Dace      *dace;
char       dforce[100];
int        addoriginal;
int        dconst;
int        rbackbone;
};
#endif
