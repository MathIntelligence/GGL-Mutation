#include"source.h"
Dipole::Dipole()
{
int i;
for(i=0;i<3;i++) xyz[i]=0;
for(i=0;i<3;i++) pot[i]=0;
db=0;
atm=0;
crg[0]=0;
crg[1]=0;
}
Dipole::~Dipole()
{
if(crg[0]) delete crg[0];crg[0]=0;
if(crg[1]) delete crg[1];crg[1]=0;
}
