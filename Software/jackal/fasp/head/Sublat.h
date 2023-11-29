#ifndef _Sublat
#define _Sublat

class  Lattice;
class  Cell;
class Sublat
{
public:

Sublat(Lattice *);
~Sublat();
void ready(float);
void ready();
void ready(Chn *);
int getcell(Atm *,float,int);
int getcell(Chn *);
void getmap(Chn *,int);
void getmap(Res *,int);
void getmap(int);
void addmap(Atm *,float,int);
int clash(Res *,int,int);
int clash(Res *,int);
int clash(Atm *);
//member data

Lattice *lattice;
Atm    *atm;
Atm    **obtain;
Sublat *next;
float grdsiz; //grid step
int   mgrd[3];
float mlen[3]; //maximum radius
int *oncell;//the energy on cell,mapped or not
int total;    //total number of this lattice
int flag;
int nget;
};

#endif
