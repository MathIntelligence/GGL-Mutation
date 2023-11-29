#include"source.h"

DistPopularBin::DistPopularBin()
{
dist=0;
next=0;
}

DistPopularBin::~DistPopularBin()
{
	if(dist) delete dist;
	if(next) delete next;
}

DistPopular *DistPopularBin::getnearnext(){
	
	DistPopularBin *t;

	for(t=this;t;t=t->next) {
		if(strcmp(t->name,"internal")==0) return t->dist;
	}
	return 0;
}

DistPopular *DistPopularBin::getDistPopular(char *ss){

        DistPopularBin *t;
	if(ss==0) return 0;
        for(t=this;t;t=t->next) {
                if(strcmp(t->name,ss)==0) return t->dist;
        }
        return 0;
}
void DistPopularBin::setnearnextbond(){

        DistPopularBin *t;

        for(t=this;t;t=t->next) {
                if(strcmp(t->name,"internal")==0) {
			t->dist->setbond();
		}
        }
}

DistPopular *DistPopularBin::getfaraway() {

	DistPopularBin *t;

        for(t=this;t;t=t->next) {
                if(strcmp(t->name,"external")==0) return t->dist;
        }

	return 0;
}

DistPopular *DistPopularBin::gethbond() {

        DistPopularBin *t;

        for(t=this;t;t=t->next) {
                if(strcmp(t->name,"hbond")==0) return t->dist;
        }

        return 0;
}

