#ifndef _ModAlgnFmt
#define _ModAlgnFmt

class Tres;
class Rcs;

class ModAlgnFmt {

public:

ModAlgnFmt();
~ModAlgnFmt();
void setmutatesec();
int   getstructurenumber();
int   getcreditofhbond(Atm *a,Atm *a0);
float *gethbondbounds(Res *r,Res *r0);
float *gethbondbounds(Atm *a,Atm *a0);
float getalgnweight(int);
void  setsitescore();
float getalgnscore();
int   isconservedpair(int n,int n0);
void setmutateresn();
void setsecstruct();
void setconservedres();
void setallconservedres();
void takeofzero();
float *getfardistbounds(Atm *a,Atm *a0);
float *getdistbounds(Atm *a,Atm *a0);
float *getdistbounds(Atm *a,Atm *a0,float low,float high);
int  ishbondexist(Res *r, Res *r0);
ModAlgnFmt *findfirstsequencefmt();
ModAlgnFmt *findsequencefmt();
void printoutnoseqmesg();
void printoutallnoseqmesg();
void clearemptyseq();
void setmutatesidepdb();
ModAlgnFmt *getparentModAlgnFmt();
ModAlgnFmt *getsequenceModAlgnFmt();
ModAlgnFmt *getlongestsequence();
void setdefaultqueryseq();
void setalldefaultseq();
ModAlgnFmt *getrootModAlgnFmt();
void setdefaultstructseq();
void checklength();
void checkalignmentlength();
void setseqnstructure(Pdb *);
void setseqnstructure();
void getallstructure();
void readModelAlgnFmt(char *fileName);
void setresn();
void deleteallchain();
void deleteallgap();
void getstructure();
void createstructure(char *);
void createstructure();
void deletegap();
void removenonstandard();
void removeallnonstandard();
void deletechain();
void setstartend();
void setmutatesidechainpdb();
void prepare();

//data member
char *fileName; //filename
char *code;     //pdb code
int flag;       //indicate different choices;
int start;	//start id
int end;        //end id
char cid;       //chain id;
char *seqn;     //sequence no gap
char *seqngap;  //original sequence with gap
Res  **resn;    //each sequence char corresponding to each residue
int  *match;    //seqngap to seqn
int  *compare;  //seqn to seqngap;
float *sitescore;
char source;    //N,X or M
char *token;    //sequnece or structure
Pdb  *pdb;       //pdb
Pdb  *mutatepdb; //sidechain mutated pdb
ModAlgnFmt *more; //the same file list
ModAlgnFmt *last; //in the same more list
ModAlgnFmt *up;   //go up to the last next;
ModAlgnFmt *next; //in the different file list
};
#endif
