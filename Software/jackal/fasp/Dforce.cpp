#include"source.h"
Dforce::Dforce()
{
next=0;
flag=0;
}
Dforce::~Dforce()
{
 if(next) delete next;
 next=0;
}
