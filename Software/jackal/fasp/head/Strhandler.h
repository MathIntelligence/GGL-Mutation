#ifndef _Strhandler
#define _Strhandler

class Strhandler
{
public:
  Strhandler();
 ~Strhandler();
static void clearbadascii(char *);
static void clearendemptyspace(char *);
static char **clearendchar(char **,char *);
static char *lastindexof(char *,char *);
static int gettotnum(float **s);
static void ungets(char *,FILE *);
static int  notuntilstrlast(char *,char *);
static int  notuntilcharlast(char *,char *);
static int  notuntilstr(char *,char *);
static int  notuntilchar(char *,char *);
static int    getasciinum(char *); 
static void   setzero(char *,int);
static char **myallocat(char **,int);
static  char *myallocat(char *,int);
static char * gethttptime(long int );
static char * straddup(char *,char *,char *);
static char * straddup(char *,char *);
static char * assemblefilepath(char **,char *);
static char * assemblefilepath(char **);
static char **attacharray(char **,char **);
static char **addchararray(char **,char **);
static char **addchartoarray(char **,char *);
static char **delchartoarray(char **,char *);
static int   gettotnum(char **);
static int   gettotnum(int **);

static char **shrink(char **,int);
static int   writeout(char **ss,FILE *fp);
static int   ifexist(char **,char *);
static int   filefound(char *);
static char *getfile(char *);
static char *getpath(char *);
static int   attachfilebylineifnoexist(char *,char *);
static int   attachfile(char *,char **);
static int   attachfilebyline(char *,char *);
static int   attachfile(char *,char *);
static int   pushfile(char *,char **);
static int   pushfile(char *,char *);
static int   findid(char *,char *);
static int   findid(char **,char *);
static int   prnfile(char *,char *);
static int   prnfilebyline(char *,char **);
static int   prnfilebyline(FILE *, char **);
static char **getdelfilebyline(char *,char *);
static int   delfilebyline(char *,char *);
static char *opnfile(char *);
static char **opnfilebylinesimple(char *); 
static char **opnfilebyline(char *);
static char **opnfilebyline(char *,int,int);
static char **opnfilebyline(FILE *);
static char *opnfile(FILE *);
static char *clearchar(char *tmp,char *);  //kick out all space
static char *clearfirstchar(char *tmp,char *); //
static char *clearlastchar(char *tmp,char *);
static char *clearendchar(char *tmp,char *);
static char *lowercase(char *tmp);
static char *uppercase(char *tmp);
static int   stlen(char *);
static char **opnarray(int);
static int  **opnarrayint(int);
static char *gettimeymd(int);
static char *gettimeymd_s(int);
static char **opnstr(char **tmp); //create sting array
static int  *opnint(int *,int);
static int  *opnint(int);
static char *opnstr(int);
static char *opnstr(char *tmp);  // create string
static char *changestr(char *tmp,char *,char *); //clear \n
static char *changestr(char *tmp,char *,char *,int);
static int  *intdel(int *);
static int  **intdel(int **);
static float *floatdel(float *);
static float **floatdel(float **);
static char *strdel(char *);
static char **strdel(char **);
static char *strcasestr(char *,char *);
static int  *pairitem(char *,char *);
static char **pairbytokensimple(char *,char *);
static char **pairbytoken(char *,char *);
static char *cutflag(char *flg,char *add);
static char *addflag(char *flg,char *add);
static int   isstdchar(char *);
static int   isstdchar(char);
static int   isgoodchar(char);
static int   islowerchar(char *);
static int   isupperchar(char *);
static int   ischar(char *);
static int   isfloat(char *);
static int   isinteger(char *);
static int   isempty(char *,char *);
static char *getpart(char *,int,int);
static char *getstr(int);
static char *getstr(float);
static int   strncasecmplast(char *,char *,int);
static int   mystrncasecmp(char *,char *,int);
static int   strncmplast(char *,char *,int);
static int   isemailendchar(char );
static char* getbetweenwithend(char *,char *,char *);

static char* getbetween(char *,char *,char *);
static char* getbetweeninclude(char *,char *,char *);
static int   getposition(char *,char *);
static char *insert(char *,int,char *);
static char *insert(char *,char *,char *);
static char *insertbefore(char *,char *,char *);
static char *cuttextbetween(char *,char *,char *);
static char *cuttextinclude(char *,char *);
static char *strstrlast(char *,char *);
static int   viewablechar(char *);
static char *translate(char **,char *);
static int   istimeright(int);
static int   istimeright(int,int);
static int   istimeright(int,int,int);
static char *changestrloop(char *,char *, char *,int ,int );
static char *getsegment(char *,int,int);

private:
static void init();
};
#endif
