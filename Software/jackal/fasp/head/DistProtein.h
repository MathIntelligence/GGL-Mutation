#ifndef _DistProtein
#define _DistProtein

class DistProtein
{
public:
DistProtein();
~DistProtein();
void buildneighboratomdatabase(char *file);
void testhbonddatabase(char *);
void indONbuildallhbonddatabase(char *file);
void buildallfarawaytruedatabase(char *file);
void buildnewloopdatabase(char *file);
void buildallnewloopdatabase(char *file);
void buildthenewloopdatabase(char *file);
void buildloopdatabase(char *file);
void buildhbonddatabase(char *file);
void buildallhbonddatabase(char *file);
void reordernearnext();
void buildnearnextdatabase(char *);
void buildnearnextfulldatabase(char *);
void buildssyesdatabase(char *);
void buildssnodatabase(char *);
void buildfarawaydatabase(char *);
void buildfarawaytruedatabase(char *);
void buildnearnextfulldatabaseh(char *file);
void indbuildallhbonddatabase(char *file);
void writetresnum(FILE *);
Tres *tres;
DistDatabase *distdb;
int tresnum[200];
};
#endif
