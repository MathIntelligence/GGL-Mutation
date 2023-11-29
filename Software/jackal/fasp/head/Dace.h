#ifndef _Dace
#define _Dace

class Dace
{
public:
Dace();
~Dace();
void create(float *,Res *,int);
void update1(float *);
void update0(float *);
void update(float *);
float getrmsd(Dace *,Dace *);
void calcallenergy();
void create(float **,Res *,int);
void create(Res *,int);
void ready(Pdb *,int);
void clear();
Chn *segchn;
int flag;
char cid;
float energy;
float energyeff;
float weight;
float rmsd;
float rmsdeff;
Pdb *pdb;
Atm **atoms;
char force[100];
Disc *disc;
Dace *next;
Dace *more;
float cutoff;
float coeff;
};
#endif
