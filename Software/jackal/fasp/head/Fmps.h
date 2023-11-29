#ifndef _Fmps
#define _Fmps

//loop prediction program
class Pdb;
class Chn;
class Lattice;

class Fmps
{
public:
           Fmps();
          ~Fmps();
void       scap();
void       initial();
void       ready();  //preparing for loop
float      segmin(int); //loop minimizer
void       setoff(float **,float);
int        create();
void       create(Res *,int);        
void       minall(int);
void       minimize(int);
int        setend(int);
void       minseg(int);
float      clashall();
float      clashall(Res *,int);
void       printhelp();
float      rmsd(int);
float      clash(Res *,int from,int to);
float      clash(Res *,int to);
float      clash(Res *);
float      clash(Atm *);

//data member
char      *force;
Pdb       *pdb;
int        start;
int        end;
Pdb       *segment;
Lattice   *lat;
float      step;
float      cutoff;
int        direct; //direction
int        pdbout;
float      outlet;
int        flat;  //control minimization
int        cycle;
char	   cid;
int        tep; 
int        flag;
int        side;
};
#endif
