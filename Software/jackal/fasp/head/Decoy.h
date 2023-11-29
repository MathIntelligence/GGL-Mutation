#ifndef _Decoy
#define _Decoy

class Decoy
{
public:
Decoy();
~Decoy();
void createdecoyfly(int,int,Chn *);
void createdecoy(int,int);
void randomskip(int);
void write(char *);
Chn *randmz(Chn*,int);
Chn *native;
Chn *decoy;
int  algn;
};
#endif
