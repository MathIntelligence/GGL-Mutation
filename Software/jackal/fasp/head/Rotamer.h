#ifndef _Rotamer
#define _Rotamer

class Rotamer
{
public:
Rotamer();
~Rotamer();
void calcrotamer(int argc, char **argv);
void calcrotamer0(Res **rrr, int nnn,int step);
void calcrotamer1(Res **rrr, int nnn,int step);
void help();
int selfenergy;

};

#endif
