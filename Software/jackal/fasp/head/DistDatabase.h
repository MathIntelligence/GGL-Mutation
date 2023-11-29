#ifndef _DistDatabase
#define _DistDatabase

class DistDatabase
{
public:
DistDatabase();
DistDatabase(Tres *,char *, char *,int);
DistDatabase(Tres *);
~DistDatabase();
void write(FILE *fp,int ff);
void addbounds(float);
void addbounds(Tres *,char *, char *,int,float);
void write(FILE *);
float dist[2];  
float *popular;
char *aim,*des; // the atom type;
int resn; //the destination residue sequence id -1,+1, etc
Tres *tres;
float cutoff;
int size;
DistDatabase *more;
DistDatabase *next;
};
#endif
