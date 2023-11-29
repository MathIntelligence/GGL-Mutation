#ifndef _NearDistCut
#define _NearDistCut

class Bound;
class DistPopular;
class NearDistCut
{
public:
NearDistCut();
~NearDistCut();
NearDistCut *next;
void cut();
void read(char *,char *,char *);
void clear();
Bound *bound;
Pdb   *pdb;
DistPopular *dist;
};

#endif
