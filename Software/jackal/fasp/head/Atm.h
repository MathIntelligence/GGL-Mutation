#ifndef _Atm
#define _Atm

class Contact;
class Res;
class Tatm;
class Lattice;

class Atm{
public:

Atm();
Atm(Atm *);
Atm(Tatm *);
~Atm();
//member function
Atm *ischildatm(char *);
Atm  *cloneatm();
void deletecontact();
int dsspdefinedhbond(Atm *s);
void writerescard(FILE *);
int isnewbondnear(Atm *s);
int isbondnear(Atm *s);
int ishbondexist(Atm *s);
int ishbondexist(Atm *s,float,float,float);
int getnumberhbonds();
float gettorsionangle();
void initial();
void writeold(FILE *fp,char);
void writeold(FILE *fp);
void write(FILE *fp);  //write the atom to fp
void write(FILE *fp,int); //write the atom in different format
int isbond(Atm *s);       //s is bonded to this atom or not?
int isbond(char *s);      //char s is bonded atom or not
int isnear(Atm *s);
int isnear(Atm *s,int i);    //detect how far away in bond of s and this
int allnear(int dis,int flg); //flg 0: +-0 res;1, none;-1, +/- res
float ishbond(Atm *s,float,float,float); //is hbond?
void writehbondparm(Atm *s);
void transfer(Atm *s);            //set the coordinat to that of s
void transfer(float *s,int flg); //set the coordinate based on flg
void transfer(float s); //set the coordiante to s
float dihedral();
float dihedral(FILE *fp); 
float dihedral(FILE *fp,int m);//m==-1 all,m>=0 one
float clash(Atm *);  //clash energy between a and this
float torsion(int);  //torsional eneregy of this->parent...
float surface(Lattice *lat,int m,float probe);
void  bondleng();
//member data
Res  *res;
Atm  *next;    
Atm  *bond[4]; //bonded atoms
Atm  *hbond[2]; //hydrogen bond
Atm **near;  //the neighboring atoms
int  *nnear; // how far the neighboring
Tatm *tatm;   // standard atom pointer
Contact *contact;
char  name[5]; //atom's name
int   id;  //the id in tres
int   id0; //the sequential number of atom, no gap in whole protein
int   oid;
float bondlen[4]; //bond length
float xyz[3]; //coordinates
float occup;  //occupancy
float bfact;  //B-factor
float chi;   //atoms chi angle
float energy;
float area; //atoms solvent accessible area
int   good; //1 useful, 0 ignore
int   mask;
int   flag;
int   nemp;
float *temp;
};
#endif
