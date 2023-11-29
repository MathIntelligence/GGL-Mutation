#ifndef _EngCoeff
#define _EngCoeff

class Atm;
class EngCoeff
{
public:
EngCoeff();
~EngCoeff();
void printoutbound(float,float,float,float);
void clear();
void setnewcoeff(float,float);
void setcoeff(float,float);
void setdots(float x,float y,float w);
void setcoeff(float t);
float getvdwvalue(float);
float getvdw1storder(float);
float getvdw2storder(float);
float getcrgvalue(float,float);
float getcrg1storder(float,float);
float getcrg2storder(float,float);
float getdelta(float,float,float,float);
float getcrgvalue(float);
float getcrg1storder(float);
float getcrg2storder(float);

void printnewout();
void printout();
void printoutorder(float delt);
void printoutcharge(float);
float getmatch(); 
float get2match(); 
float getdelta(float,float,float,float,float,float,float,float,float);
int flag;
float *dots;
float *coeff;
int size;
EngCoeff *next;
};

#endif
