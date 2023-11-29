#ifndef _ResBound
#define _ResBound

class ResBound
{
public:
ResBound();
~ResBound();
void setresidue(Tres *t);
DistBound *distbound;
Tres *tres;
ResBound *next;
};
#endif
