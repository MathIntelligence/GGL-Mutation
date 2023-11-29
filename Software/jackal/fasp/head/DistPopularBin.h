#ifndef _DistPopularBin
#define _DistPopularBin

class DistPopularBin
{
public:
DistPopularBin();
~DistPopularBin();
void setnearnextbond();
DistPopular *getDistPopular(char *ss);
DistPopular *getnearnext();
DistPopular *getfaraway();
DistPopular *gethbond();
char *name;
DistPopular *dist;
DistPopularBin *next;
};
#endif
