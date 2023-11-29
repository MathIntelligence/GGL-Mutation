#ifndef _MinBound
#define _MinBound

class Atm;
class MinBound
{
public:
MinBound();
~MinBound();
void optimize();
int flag;
Bound *bound;
MinBound *next;
};

#endif
