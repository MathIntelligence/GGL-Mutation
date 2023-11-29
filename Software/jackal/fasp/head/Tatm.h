#ifndef _Tatm
#define _Tatm

class Tres;
class Eng;
class SimAtm;
class Tatm
{
public:
Tatm();
~Tatm();
int isnear(Tatm*,int);
void bondleng();
int getnumberhbonds();
int id;
Tres *tres;
Tatm *next;
Tatm *bond[4];
float bondlen[4]; //bond length;
int isbackbone();
int rotate;   // rotatable
int balance;  //1 balance,0 no
int hbond;    // 1 donor,2 acceptor, 3 both, 0 none of these 
int nbond;    // number of bond
int polarkeep; //polar keep
int keep;     // extended or not   
int nonp;     // belong to nonpolar residue
char name[5]; // atom name
char type[5]; //atom type name
float xyz[3]; 
float area;   //maximum surface area of extended
int ntat;     //rotable sequential number
int rotm;
int ishn;     // is HN or CD in pro
int ispolar;  // 1 nonpoalr;0: polar atom
int ring;
int sidecenter;
SimAtm *sim; //pointer to simatm in TRES
Eng *eng;     //energy parameter
};

#endif
