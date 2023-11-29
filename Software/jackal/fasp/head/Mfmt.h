#ifndef _Mfmt
#define _Mfmt

class Mfmt
{
public:
Mfmt();
~Mfmt();
void setmfmt(StrFmt *);
void setfmt(StrFmt *);
void printhelp();
char **aln;
char **name;
int sec;
void writeout(int);
};
#endif
