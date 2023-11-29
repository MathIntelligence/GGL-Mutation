#include"source.h"

Dssp::Dssp()
{
	key=-1;
	next=0;
}

Dssp::~Dssp()
{
	if(next) delete next;
}
//
/*
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
*/
void Dssp::add(int n) {

	Dssp *x;

	backbone c;
	c.access=1;
	cerr<<c.access<<endl;
	if(key==-1) {
		key=n;
		return;
	}

	for(x=this;x->next;x=x->next);

	if(x==0) {
		cerr<<"could not add to hashtable, error"<<endl;
		return;
	}

	x->next=new Dssp();
	x->next->key=n;
}
