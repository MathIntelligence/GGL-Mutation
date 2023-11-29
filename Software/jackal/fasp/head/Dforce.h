#ifndef _Dforce
#define _Dforce

class Dforce
{
public:
Dforce();
~Dforce();
float *calcenergy(float **);
float *ensort(float **);
Dforce *next;
Pdb *pdb;
int flag;
char force[100];
float cutoff;
};
#endif
