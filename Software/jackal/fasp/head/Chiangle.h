#ifndef _Chiangle 
#define _Chiangle

class Chiangle {

public:

Chiangle();
~Chiangle();
void setuprecord(int,float);
void buildallstructure();
void buildstructure();
void readone(char *file);
void read(char *file);
Chiangle *get(int);
float getangle(int,int);
char *cleardoublespace(char *);
float *parseangle(char *);
int parseanglenum(char *);
void addanyway(Chiangle*);
void add(float *,int);
void add(Chiangle*);
int  getsize();
int ifexist(Chiangle *);
void size(int);
void resize(int);
char *name;    //sequence residue
float **angle; //the angle
int  *numa;    //the number of chi angles of each residue
int *oid;
char *chain;
Chiangle *next;
int number;
int omega;
int rotdir;
int endres; //end residue ids
int wavedir;
char *filename;
Pdb *pdb;
};
#endif
