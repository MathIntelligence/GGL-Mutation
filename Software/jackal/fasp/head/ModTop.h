#ifndef _ModTop
#define _ModTop

class Pdb;
class StrFmt;
class ModTop
{
public:
ModTop();
~ModTop();

void sethbondconstraint(StrFmt*);
void setmodel();
void setbound();
//not important
void setbound(StrFmt*);
void addhbondconstraint(HBondList *);
StrFmt *strfmt;
Bound *bound;
int flag;
Pdb *model;
};
#endif
