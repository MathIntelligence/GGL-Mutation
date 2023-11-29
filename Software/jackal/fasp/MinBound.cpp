#include"source.h"
MinBound::MinBound()
{
next=0;
bound=0;
flag=0;
}

MinBound::~MinBound()
{
 if(next) delete next;
 next=0;
}
