#ifndef _AtmGeom
#define _AtmGeom

class AtmGeom
{
public:
AtmGeom();
~AtmGeom();

void addbounds(float *x,float d,float off);
void addbounds(float *x,float d,float off,int w);
float evaluatemin();
float evaluate();
float evaluaterange();
float calcscore(float *);
void findmincoo(float);
void findmincoo(int,float);
void findcoo();
void findautocoo();
float dist[1000];  //distance constraint
float xyz[1000];   //coordinate
float weight[300];
float coo[3]; 
int num;
};
#endif
