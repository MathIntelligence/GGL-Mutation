#ifndef _Seqchn
#define _Seqchn

class Seqchn{

public:
Seqchn();
~Seqchn();
void write(FILE *fp,int flag);//flag,0 profile,1 sequence only;
char *seqn;
char *name;
Seqchn *next;
int flag;
};

#endif
