#ifndef _ATOM
#define _ATOM


#include <iostream.h>

#include "vector3d.h"

namespace Troll { 

class Residue;


class AtomName {

public:
  AtomName(char*,Residue*);

  int operator==(char*);
  int operator==(AtomName&);
  int operator=(char*);
  operator char*() { return name; }

  char name[5];
  Residue *residue;

} ;



class Atom {

public:


// Construction/destruction.

  Atom(Residue*,char*,float b=0.0);


// Atom identity.

  AtomName name;
  Residue *residue;
  Vector pos;
  int Is(char*);

// Throwaway variable

  int param;


// Atom properties

  float Area(unsigned long flag=A_Current);
  float radius,charge,bfactor,area,occupancy;

  unsigned long style,color;

} ;
  

::ostream& operator<<(::ostream&,Atom&);
::ostream& operator<<(::ostream&,AtomName&);

}



#endif


