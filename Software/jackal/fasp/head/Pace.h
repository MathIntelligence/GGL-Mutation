#ifndef _Dace
#define _Dace

class Dace
{
public:
Dace();
~Dace();
void update(float *);
float getrmsd(Dace *,Dace *);
void calcallenergy();
void create(float **,Res *,int);
void create(Res *,int);
void Dace::ready(Pdb *,int);
void clear();
Chn *segchn;
int flag;
int cid;
float energy;
float weight;
float rmsd;
Pdb *pdb;
Atm **atoms;
Disc *disc;
Dace *next;
Dace *more;
float cutoff;
};
#endif
