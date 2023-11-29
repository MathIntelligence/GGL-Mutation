#ifndef _Stralg
#define _Stralg

class Chn;
class Res;
class Stralg
{
public:
Stralg();
~Stralg();
float superimposestill(int flg);
float superimposesimple(Res *,int,Res *,int,int);
float simplermsd(Chn *,Chn *,int);
float superimpose(Pdb*,Pdb *,int);
float superimpose(Chn *,Chn *,int);
float superimpose(int);
void  randmz(Chn *,float,float,float);
float *getcenter(Res **,int);
void setcenter(float *,Res **,int);
float *center(Res **,int);
float *xyza;
float *xyzb;
Res **alga;
Res **algb;
int   flag; //1: all; 2: chain only 3: residue only
int   onlyaligned;
};
#endif
