#ifndef _Qsort
#define _Qsort

//including most used subroutines
class Qsort
{
public:
Qsort();
~Qsort();
void sort(double *,int,int *);
void sort(float *,int,int *); 
void back(float *,int,int *);
int  slnpd(float **,float,int,int);
int  slnpd(float **matrix,int *order,int wide,int length);
double matrix(float **,int);
};
#endif
