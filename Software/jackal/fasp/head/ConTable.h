#ifndef _ConTable
#define _ConTable

class ConTable
{
public:
ConTable();
~ConTable();
void add(ConTable *c);
float getsqrtvalue(float);
ConTable *setsqrt(float,float);
ConTable *setsqrt(float,int);
float *table;
float distance;
float resolution;
char *name;
int  size;
ConTable *next;
};

#endif
