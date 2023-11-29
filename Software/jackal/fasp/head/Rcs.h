#ifndef _Rcs
#define _Rcs

class Rcs{
public:
  Rcs(char*); 
  ~Rcs();
  char *operator[](char *);
  int number;
  char **name; 
  char **path;
};
#endif
