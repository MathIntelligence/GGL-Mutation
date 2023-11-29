#ifndef _Cell
#define _Cell

class Atm;
class Cell
{
public:
Cell();
~Cell();
Atm *atm;
Cell *next;
int used;
int id;
int flag;
};

#endif
