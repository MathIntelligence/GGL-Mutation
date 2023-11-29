#ifndef _HBondList
#define _HBondList

class HBondList
{
public:
HBondList();
~HBondList();
void add(Atm *a,Atm *b,float,int);
void fill(Atm *a,Atm *b,float,int);
HBondList* get(Atm *,int);
Atm *donor;
Atm *acceptor;
float energy;
int type;  //h bond, salt bridge, ss bond
HBondList *next;
};
#endif
