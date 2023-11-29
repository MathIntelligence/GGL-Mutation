#include"source.h"
Cell::Cell()
{
atm=0;
next=0;
id=0;
used=-1;
flag=0;
}
Cell::~Cell()
{
 //if(next) delete next;
 next=0;
}
