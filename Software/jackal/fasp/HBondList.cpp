#include"source.h"

HBondList::HBondList()
{
	donor=0;
	acceptor=0;
	energy=0;
	type=0;
	next=0;
}

HBondList::~HBondList()
{
	if(next) delete next;
}

void HBondList::add(Atm *a,Atm *b,float e,int t) {

	if(donor==0&&acceptor==0) {
		fill(a,b,e,t);
		return;
	}

	HBondList *c;

	for(c=this;c->next;c=c->next);

	c->next=new HBondList();

	c->next->fill(a,b,e,t);
}

void HBondList::fill(Atm *a,Atm *b,float e,int t) {

	donor=a;
	acceptor=b;
	energy=e;
	type=t;
}

HBondList *HBondList::get(Atm *a,int n) {

	HBondList *c;

	int i=0;
        for(c=this;c;c=c->next) {

		if(c->donor!=a&&c->acceptor!=a) continue;

		if(n==i) return c;
	
		i++;
	}

	return 0;
}
