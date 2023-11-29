#include "source.h"



Atom::Atom(float x,float y,float z) {
	pos.x=x;
	pos.y=y;
	pos.z=z;	
	radius=0;
	id=0;
}

Atom::Atom() {

	pos.x=0;
        pos.y=0;
        pos.z=0;
        radius=0;
	id=0;
}

Atom::Atom(float x,float y,float z,float t) {
        pos.x=x;
        pos.y=y;
        pos.z=z;
        radius=t;
	id=0;
}


