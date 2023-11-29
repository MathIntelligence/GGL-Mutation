#ifndef _Loopy
#define _Loopy

//loop prediction program
class Pdb;
class Chn;
class Lattice;

class Loopy
{
public:
           Loopy();
           Loopy(Loopy *);
          ~Loopy();
void       initial();
void       ready();  //preparing for loop
int	       randmzChiangle(int);
void       checkomega(FILE *);       
int	   checkomega();
float      closeg(Loopy *,Loopy *); // tweak two segment to consistency
float      closeg(); // tweak segment to constraint
float      *trim(float **,char *,int haa,int);
float    **around(float **);
float    **zipper(int,float,int);  //random loop generator
float     *expand(float **,int flg); //expanding nearby
float      expand(float *,int flg);  //expanding nearby
float    **expand(float *,int,int,int); //loop searching around neighbor and return all found
float     *peel(float **,char *flg,int got);
float    **fill(float **,float);
float    **trim(float **);
float     *trim(float **,float *,float d); //trim those large than d
float     *trim(float **,float *,int haa);//trim those no more than part
float     *trim(float **,char *,int haa); //trim those no more than part
float    **generate(float **,int,float);
float    **predt0();
float    **predt(int);
float    **copy(float **);
float    **copy(float **,int);
float    **refill(float **,char *flg,int need);
void       testclosure();
float      segmin(int); //loop minimizer
float      segmin0(float,int); //loop minimizer
int        noclose(float **xyz,float d);  //noclose than d
float      noclose(float **xyz,float *xyz0); //find closest
float      minimize(float *xyz_all,int flg);
float     *minimize(float **xyz_all,int flg); //minimization
float    **cross(float **xyz,int cot,float cut);   //loop cross
float      ensort(float *xyz,char *flg); 
float     *ensort(float **xyz,char* flg); //energy sort
void       setoff(float **,float);
int        setend(int); //satisfying end meet
int        readmore(char *,char *);  
int        create();
void       create(Res *,int);        
void       create(char *seqn);
void       randmz(int seed,int rand);
void       randmz(int);
void       writeseg(float **,char *,int);
void       write(float **,char *,int);
void       printhelp();
float      scpred(char *);
void       crtmore(char *);
float     *scap(float **,char *flg);
float      scap(char *flg);
float      rmsd(int);
float      clash(Res *,int from,int to,char *);
float      clash(Res *,int to,char *);
float      clash(Res *,char *);
float      clash(Atm *,char *);
void      randmzHelix(int,int);
void 	  randmzSheet(int,int);
int      **allconf(int);
float    **zipper();

//data member
Pdb       *pdb;
int        start;
int        end;
Pdb       *segment;
Lattice   *lat;
float      step;
float      cutoff;
Loopy     *last;
Loopy     *next;
int        direct; //direction
int        arbt;  //arbitrary number of first attemp
int        id;
float    **near;
float    **init;
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
//
int        fapr; 
Chiangle  *chiangle;
};
#endif
