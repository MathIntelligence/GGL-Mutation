#ifndef _Map
#define _Map
class Chn;

class Map
{
public:
Map();
~Map();
void ready(Pdb *);
void clear();
void reform(float *,int,int,int);
int cover(float *);
int cover(Atm *);
void puton();
void puton(Chn *);
void puton(Atm *);
void puton(Res *,int);
int ***node;
Pdb *pdb;
float offmin[3],offmax[3];
int   edge[3];
float grdsiz;
};
#endif
