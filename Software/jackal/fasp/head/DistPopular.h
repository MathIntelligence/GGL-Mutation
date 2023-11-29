#ifndef _DistPopular
#define _DistPopular

class DistPopular
{
public:
DistPopular();
DistPopular(Tres *,char *, char *,int);
DistPopular(Tres *);
~DistPopular();
void setbond();
DistPopular *getDistPopular(Tres *,char *, char *,int);
DistPopular *getDistPopular(Tres *);
DistPopular *getDistPopular(char *, char *,int);
float findmostpopulardistance();
float *findmostpopulardp(float);
float findmax();
void write(FILE *fp);
void writelarge(FILE *fp);
void addbounds(int,float,float,float *,int);
void read(char *);
void addbounds(Tres *,char *, char *,int,int,float,float,float *,int);
int   bond;
float dist[2];  
float *popular;
char  *aim,*des; // the atom type;
int resn; //the destination residue sequence id -1,+1, etc
Tres *tres;
int size;
int num;
char *file;
DistPopular *more;
DistPopular *next;
};
#endif
