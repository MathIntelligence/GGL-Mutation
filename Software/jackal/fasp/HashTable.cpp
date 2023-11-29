#include"source.h"

HashTable::HashTable()
{
	key=-1;
	next=0;
}

HashTable::~HashTable()
{
	if(next) delete next;
}


void HashTable::add(int n) {

	HashTable *x;

	if(key==-1) {
		key=n;
		return;
	}

	for(x=this;x->next;x=x->next);

	if(x==0) {
		cerr<<"could not add to hashtable, error"<<endl;
		return;
	}

	x->next=new HashTable();
	x->next->key=n;
}
