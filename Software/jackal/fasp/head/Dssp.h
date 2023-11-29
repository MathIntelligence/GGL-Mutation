#ifndef _Dssp
#define _Dssp

class Dssp
{
public:
Dssp();
~Dssp();
void add(int);
int key;
Dssp *next;
typedef struct backbone {
  char aaident;
  char sheetlabel, aa;
  char threelettercode;
  char ss[1];
  long partner[1];
  long access;
  double alpha, kappa;
  long atompointer, nsideatoms;
} backbone;

};
#endif
