#ifndef _MapNew
#define _MapNew
class Chn;

class MapNew
{
public:
MapNew();
~MapNew();
void ready(Pdb *);
void clear();
void reform(float *,int,int,int);
void restart();
void initial();
void deletenode();
int cover(float *);
int cover(Atm *);
int cover(float,float,float);
void puton();
void puton(Chn *);
void puton(Atm *);
void puton(Res *,int);
void putoff();
void putoff(Atm *);
void putoff(Res *,int);
void putoff(Chn *);
int ***node;
Pdb *pdb;
float offmin[3],offmax[3];
int   edge[3];
float grdsiz;
float probe;
float leg;
};
#endif
