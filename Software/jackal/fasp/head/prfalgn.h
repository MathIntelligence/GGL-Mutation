#ifndef _Prfalgn
#define _Prfalgn

class Protein;
class Prfalgn{

public:

Prfalgn(char *s);
Prfalgn(Protein *protein);
~Prfalgn();
void prfmatrice(char *);
void matrice();
void align(char *sequence);
char **output(FILE *);

//member data
Protein *protein;
int *important;
int score;

//member data
char *target;
char *query;
float **matrix;
int *routine;
char **result;
};

#endif
