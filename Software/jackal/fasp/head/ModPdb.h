#ifndef _ModPdb
#define _ModPdb

class Pdb;
class ModAlgnFmt;
class ModPdb
{
public:
ModPdb();
~ModPdb();

void sethbondconstraint(ModAlgnFmt*);
void setmodel();
void setdistbound();
void setdistbound(ModAlgnFmt*);
void addhbondconstraint(HBondList *);
ModAlgnFmt *modalgn;
DistMutateBound *distbound;
int flag;
Pdb *model;
};
#endif
