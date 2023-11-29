#ifndef _ATOM
#define _ATOM


#include <iostream.h>


class Atom {

public:

  Atom();
  Atom(float,float,float,float);
  Atom(float,float,float);
  ~Atom(){}; 

  Vector pos;
  float radius;
  int id;
} ;
  

#endif


