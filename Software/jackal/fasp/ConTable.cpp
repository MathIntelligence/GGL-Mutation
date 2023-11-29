#include"source.h"
ConTable::ConTable()
{
	table=0;
	distance=0;
	resolution=0;
	name=0;
	next=0;
	size=0;
}

ConTable::~ConTable()
{
 if(next) delete next; next=0;
 if(table) delete [] table;table=0;
 if(name) delete [] name;name=0;
}

void ConTable::add(ConTable *c) {
	ConTable *t;
	for(t=this;t->next;t=t->next);
	t->next=c;
}

ConTable *ConTable::setsqrt(float d,float r) {

	ConTable *t=new ConTable();
	

	t->distance=d;

	t->resolution=1/r;

	int n=(int)(d/r+10);

	t->name=strdup("sqrt");

	t->table=new float[n];

	for(int i=0;i<n;i++) {
		float e=(i+0.5)*r;
		t->table[i]=sqrt(e);
	}

	return t;
}

ConTable *ConTable::setsqrt(float d,int rr) {

        ConTable *t=new ConTable();

	float r=1./rr;
        t->distance=d;

        t->resolution=rr;

        int n=(int)(d*rr+10);

        t->name=strdup("sqrt");

        t->table=new float[n];

        for(int i=0;i<n;i++) {
                float e=(i+0.5)*r;
                t->table[i]=sqrt(e);
        }

	t->size=n;

        return t;
}


float ConTable::getsqrtvalue(float d){
	
	if(d<2||d>distance) return sqrt(d);
	int n=(int)(d*resolution+0.5);	
	return table[n];
}
