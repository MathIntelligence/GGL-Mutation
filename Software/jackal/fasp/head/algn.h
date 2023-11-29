#ifndef _Algn
#define _Algn

class Algn{

public:

Algn();
~Algn();

void seqgive(char *target,char *query,int flg);

//flag=0,structure guided gapping with normal scoring;
//flag=1,sequnce guessed gapping with normal scoring;
//flag=2,structure guided gapping seeking maximum buried hydrophobic core;
//       without normal scoring;
      
void alignment(short *,int gap);
char **output(FILE *,int flg);

//member data

char *target,*query;
float score;
int *routine,*prfm;
char **result;
int flag,slen;

float opncost, gapcost;
int outfmt;
int tlen,qlen;
void scoring(short *matx, float **dmn); 
float gapping(int **,int,int,int,int);
void destroy();
};


#endif
