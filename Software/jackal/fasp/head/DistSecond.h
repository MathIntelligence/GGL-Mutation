#ifndef _DistSecond
#define _DistSecond

class DistSecond
{
public:
DistSecond();
~DistSecond();
void addbounds(int i1,int i2,float a1,float a2,float a3,float w);
void addbounds(float a1,float a2,float a3,float w);
DistSecond * createDistSecond(int i1,int i2);
DistSecond * findDistSecond(int i1,int i2);
float *dist;
int    size;
int    ndist;
int    first;
int    second;
DistSecond *next;
};
#endif
