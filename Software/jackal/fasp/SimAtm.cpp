#include"source.h"
SimAtm::SimAtm()
{
next=0;
id=0;
flag=0;
charge=0;
radius=0;
epslon=0;
name=' ';
}
SimAtm::SimAtm(char nc,int n,float c,float r,float e){
next=0;
id=n;
flag=0;
charge=c;
radius=r;
epslon=e;
name=nc;
}

SimAtm::SimAtm(SimAtm *s) {
next=0; 
id=s->id;
flag=s->flag;
charge=s->charge;
radius=s->radius;
epslon=s->epslon;
name=s->name;
}

SimAtm::~SimAtm()
{
 if(next) delete next;
 next=0;
}

SimAtm *SimAtm::ifexist(char n,float c,float r,float e){

	SimAtm *s;
	for(s=this;s;s=s->next) {
		if(s->name==n&&s->charge==c&&s->radius==r&&s->epslon==e) return s;
	}
	return 0;
}

SimAtm *SimAtm::ifexist(char n,float c,float r){

        SimAtm *s;
        for(s=this;s;s=s->next) {
                if(s->name==n&&s->charge==c&&s->radius==r) return s;
        }
        return 0;
}

SimAtm *SimAtm::ifexist(char n,float r){

        SimAtm *s;
        for(s=this;s;s=s->next) {
                if(s->name==n&&s->radius==r) return s;
        }
        return 0;
}

SimAtm *SimAtm::find(int d){

        SimAtm *s;
        for(s=this;s;s=s->next) {
                if(s->id==d) return s;
        }
        return 0;
}
int SimAtm::getnumber(){

        SimAtm *s;int n=0;
        for(s=this;s;s=s->next) {
                n++;
        }
        return n;
}

void SimAtm::write(FILE *fp) {

	fprintf(fp,"%c %5i %8.3f %8.3f %8.3f\n",name,id,charge,radius,epslon);
}
