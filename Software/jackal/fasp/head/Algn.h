#ifndef _Algn
#define _Algn

class Algn{

public:

Algn();
~Algn();

void ready();
void defalgnid();
void defalgnidz();
int  getroutelength();
void setsequence(char *target,char *query);
void setmatrix(float *);
float *getscoretable(char *);
float calctotalscore();
void alignment();
void defalgn();
void defalgnwithidmat();
float **getscore(float *matx); 
float getgapvalue(int **,int,int,int,int);
char **output(FILE *);
int   getorder(char);
//member data

char   *target;
char   *query;
int    *routine;
float **matrix;
float   opncost;
float   gapcost;
char   *aaorder;
int     gap;
};
#endif
