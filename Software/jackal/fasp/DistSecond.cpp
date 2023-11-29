#include"source.h"

DistSecond::DistSecond()
{
//weight=0;
dist=0;
size=0;
ndist=0;
first=0;
second=0;
next=0;
}

DistSecond::~DistSecond()
{
	//if(weight) delete [] weight;
	if(dist) delete [] dist;
	if(next) delete next;
}

 

void DistSecond::addbounds(float a1,float a2,float a3,float w) {
        
	if(ndist+1>size) {
		int n=ndist+100;
        	dist=(float *)realloc(dist,4*sizeof(float)*n);
        	size=n;
	}
	dist[4*ndist]=a1;       
        dist[4*ndist+1]=a2;
        dist[4*ndist+2]=a3;
        dist[4*ndist+3]=w;
        ndist++;

}
 
void DistSecond::addbounds(int i1,int i2,float a1,float a2,float a3,float w) {

	
	if(i1>i2) {
		int i3;
		i3=i1;
		i1=i2;
		i2=i3;
	}

	if(first==0&&second==0) {
		first=i1;second=i2;
		addbounds(a1,a2,a3,w);
	}
	else if(first==i1&&second==i2) {
		addbounds(a1,a2,a3,w);
	}
	else {
		DistSecond *t=findDistSecond(i1,i2);
		if(t==0) {
			t=createDistSecond(i1,i2);	
		}
		t->addbounds(i1,i2,a1,a2,a3,w);
	}
}

DistSecond *DistSecond::createDistSecond(int i1,int i2) {
	
	DistSecond *t;

	for(t=this;t->next;t=t->next);

	t->next=new DistSecond();

	t->first=i1;
	t->second=i2;
	return t;
} 

DistSecond *DistSecond::findDistSecond(int i1,int i2) {

	DistSecond *t;
	for(t=this;t;t=t->next) {
		if(t->first==i1&&t->second==i2) return t;
	}
	return 0;
}
