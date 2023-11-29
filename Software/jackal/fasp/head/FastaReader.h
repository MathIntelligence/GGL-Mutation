#ifndef _FastaReader
#define _FastaReader

class Tres;
class Rcs;
class FastaReader {

public:

FastaReader();
~FastaReader();

void readFastaAlign(char *fileName);
int  num;
char *fileName;
char **name;
char **seqn;
FastaReader *next;
};

#endif
