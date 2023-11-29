#ifndef _Dipole
#define _Dipole

class Charge; 
class Dipole
{
public:
Dipole();
~Dipole();
float xyz[3];
float pot[3];
float db;
Atm *atm;
Charge *crg[2];
};

#endif
