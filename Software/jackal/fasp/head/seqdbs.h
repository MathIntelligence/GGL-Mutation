#ifndef _Seqdbs
#define _Seqdbs

class Seqchn;
class Seqdbs{

public:
Seqdbs();
~Seqdbs();
void readdbs(char *,int flag);
//flag: 0 profile sequence;flag: 1 sequence
int chn_num;
Seqchn *seqchn;

protected:
void addseq(Seqchn *,char**,int,int );   
};
#endif
