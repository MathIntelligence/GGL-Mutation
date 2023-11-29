#ifndef _PdbFix
#define _PdbFix

class PdbFix
{
public:
PdbFix();
~PdbFix();
void printhelp();
void printhelpfix();
void fixlostresidues(Pdb *,int);
void fixlostresidues(Pdb *);
int checkpdb(Pdb *s);
int  isbackboneatom(Res *first,Res *last);
Res **maxenergystem(Res *r,Res *t);
void sidechainminimal(Pdb *);
void mynewfix0(Res *first, Res *last);
void pushaway(Pdb *);
void tuneheader(Pdb *);
void linkallres(Res *first,Res *last);
int  directlink(Res *,Res *,int );
int  islinked(Res *,Res *);
int  isbreaker(Res *,Res *);
void addsidechains(Res *r1,Res *r2);
float getsidechainrmsd(Res *r0,Res *r);
void addsidechains(Pdb *);
void myfixnow();
void minimize(Res *first,Res *last);
Chn *createchain(Chn *);
float getrmsdmore(Res *,Res *,int);
void addhnatoms(Pdb *);
void ready();
float energy(Res *,Res *);
void setatmat(Pdb*);
void setatmat();
void mynewreversefix(Chn *chn,Res *rr0,int rend);
void superimposesegment(Res *r1,Res*r2,int flg);
float getrmsd3(Res *r,Res *t,int flg);
float getrmsd(Res *r,Res *t,int flg);
float getendrmsd(Res *r,Res *t,int flg);
void mydirectnewfix(Res *first,Res *last);
void mynewfix(Chn *);
void mynewfix(Res *first,Res *last);
void mynewfixmore(Chn *);
void mynewfixmoremore(Chn *);
//
Res **findbackbonelostsegment(Res *);
void myfix0(Chn *);
void myfix();
void minimize(Chn *);
void minimize0(Chn *);
void myfix(Chn *);
void sethbonddssp(Chn *);
Res **findsegment0(Res *,int);
Res **findsegment(Res *,int);
void fixsegment0(Res **stem);
void fixsegment(Res **stem);
void fix();
void fixsave();
void fixold();
void fix(Chn *);
void fixold(Chn*);
void fixtemp();
Chn *createchain0(Chn *);
void setfixedatm(Chn *,Bound *);
void setfixedatmxyz(Chn *,Bound *);
void shiftcenter(Res *,Res *);
Pdb *pdborg;
Pdb *pdb;
PdbFix *next;
Res **atmat;
int exact;
char *seqn;
int onlybackbone;
int ctrip;
int fast;
};
#endif
