#ifndef _StrFmt
#define _StrFmt

class Tres;
class Rcs;
class Mutate;

class StrFmt {

public:

//essential functions
StrFmt();
~StrFmt();
void setseqnlimit(int);
void reallocatmatch(int);
void treatbreaker0();
void checknonesegmentcomposite();
int getblockids();
int *getresalncomposite(int);
void checkresalncomposite();
int iscodeexist(char *);
void setdefaultseqn();
void initial();
StrFmt(StrFmt*);
StrFmt *getStrFmtWithId(char);
char *getseqn();
void checkseqcomname();
void createComStrFmt();
void setupatmid0();
void beforewritefinal();
void writefinal();
void optimizecom();
void optimize();
void setshift(float,int);
void checkseqname();
void reorderstrftm();
void setotherfmtloop(int,float,int);
void printoutall(FILE *);
void setotherfmt(int,float,int);
void setotherfmt0(int);
void addhatoms();
void setsize(int);
void checkresn();
void reindexcomposite();
void setmutate(char *,int);
void exec();
void writeseqxyz();
void checklegalseqxyz();
void treatcomposite();
void printoutonly(FILE *);
void printout(FILE *);
void printhelpold();
void printminhelp();
void printrefinehelp();
void printhelp();
void checklegal();
void checkcont();
void treatallbreaker();
void combine();
void treatbreaker();
StrFmt *getStrFmt(char *);
int isbreaker();
int findpost(Res *);
int tunealign();
float getsimplealgnbound(int );
void setallseqnresn();
void setseqnresn();
int  tuneterminal();
int  iscomposite();
void tuneallalign();
float calcavgscore();
void smoothinggap();
void setmatch();
void setalloldgap();
void prepareprofix(int);
void prepare();
void preparesuperimpose();
void structuresuperimpose(char *out,int flg);
void checklegalsuperimpose();
float getsimplealgnweight(int);
void readModelAlgnFmt(char *fileName);
void setnest0();
void setnest();
void setalldssp();
char *setstate(char *);
char *setnewdssp();
char *setnewdsspsimple();
void setdssp();
void setgapvalue();
void setmatrixname(char *);
void setmutatedone();
float findmatrixminvalue();
float getnewalgnweight(int);
float getalgnweight(int);
void setnewboundscore();
void setnewsitescore();
void  setsitescore();
float getalgnscore();
void setsecstruct();
void setconservedres();
void setallconservedres();
void takeofzero();
StrFmt *findfirstsequencefmt();
StrFmt *findsequencefmt();
StrFmt *findstructurefmt();
void printoutnoseqmesg();
void printoutallnoseqmesg();
void clearemptyseq();
StrFmt *getparentStrFmt();
StrFmt *getsequenceStrFmt();
StrFmt *getlongestsequence();
void setdefaultqueryseq();
void setalldefaultseq();
StrFmt *getrootStrFmt();
void setdefaultstructseq();
void checklength();
void checkalignmentlength();
void setseqnstructure(Pdb *);
void setseqnstructure();
void getallstructure();
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

//not essential
void setmutatesidepdb();
void setmutateresn();
void setmutatesec();
int   getstructurenumber();
int   getcreditofhbond(Atm *a,Atm *a0);
float *gethbondbounds(Res *r,Res *r0);
float *gethbondbounds(Atm *a,Atm *a0);
void setmutatesidechainpdb();
float *getfardistbounds(Atm *a,Atm *a0);
float *getdistbounds(Atm *a,Atm *a0);
float *getdistbounds(Atm *a,Atm *a0,float low,float high);
int  ishbondexist(Res *r, Res *r0);
int   isconservedpair(int n,int n0);

//data member
char   *fileName; //filename
char   *code;     //pdb code
int     flag;       //indicate different choices;
int     start;	//start id
int     end;        //end id
char    cid;       //chain id;
char   *seqn;     //sequence no gap
char   *seqngap;  //new one
char   *oldgap;   //original sequence with gap
Res   **resn;    //each sequence char corresponding to each residue
int    *match;    //seqngap to seqn
int    *compare;  //seqn to seqngap;
float  *sitescore;
float   avgscore;
char    source;    //N,X or M
char    zipcode;
char   *token;    //sequnece or structure
Pdb    *pdb;       //pdb
float   gap;
char   *matrix;
Mutate *mutate;
int     tune;
int     addh;
int 	logg;
int     iscom;  //if composite 
int    	resaln[100];
int     numaln;
int     minfast;//1, smooth minimize on all atoms;2 smooth minimize only on sidechains
int     onlyrefine;
int     refstart;
int     refend;
int     alnback;
int     profix;
//float   shiftcut;
StrFmt *other; 
StrFmt *last; //in the same more list
StrFmt *up;   //go up to the last next;
StrFmt *more; //the same file list
StrFmt *next; //in the different file list
};
#endif
