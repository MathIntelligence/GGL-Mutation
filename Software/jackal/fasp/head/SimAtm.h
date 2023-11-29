#ifndef _SimAtm
#define _SimAtm

class Atm;
class SimAtm
{
public:
SimAtm();
SimAtm(SimAtm*);
SimAtm(char,int,float,float,float);
~SimAtm();
SimAtm *ifexist(char,float,float,float);
SimAtm *ifexist(char,float,float);
SimAtm *ifexist(char,float);
SimAtm *find(int);
int getnumber();
void write(FILE *);
int id;
int flag;
char name;
float charge;
float radius;
float epslon;
SimAtm *next;
};

#endif
