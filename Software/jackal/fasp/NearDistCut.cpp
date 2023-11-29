#include"source.h"
NearDistCut::NearDistCut()
{
bound=0;
dist=0;
pdb=0;
}

NearDistCut::~NearDistCut()
{
if(bound)delete bound;
if(dist) delete dist;
if(pdb) delete pdb;
}

void NearDistCut::read(char *b,char *d,char *p) {

	dist=new DistPopular();
	dist->read(b);

	bound=new Bound();
	


	pdb=new Pdb();
	pdb->read(p,'1');

	bound->model=pdb;
	bound->ready();
	bound->readdist(d);
}

void NearDistCut::clear() {

	DistPopular *c,*d;

	for(c=dist;c;c=c->next)
	for(d=c;d;d=d->more) {
		d->size=0;
	}

	for(int i=0;i<bound->ndist;i++) {

		int n1=bound->atmpair[i*2];
		int n2=bound->atmpair[i*2+1];
		Atm *a1=pdb->getatmbyid0(n1);
		Atm *a2=pdb->getatmbyid0(n2);
		if(a1==0||a2==0) continue;
		if(fabs(a1->res->id-a2->res->id)>=2) continue;

		Atm *b1,*b2;

		if(a1->res->id>a2->res->id) {
			b1=a2;
			b2=a1;
			if(b1->tatm->id>4) continue;
		}
		else if(a1->res->id<a2->res->id) {
			b1=a1;b2=a2;
			if(b1->tatm->id>4) continue;
		}
		else {
			if(a1->tatm->id>a2->tatm->id){
				b1=a2;	
				b2=a1;
			}
			else {
				b1=a1;
				b2=a2;
			}
		}
			
		
		DistPopular *t=dist->getDistPopular(b2->res->tres);	
			
		DistPopular *d;
		for(d=t;d;d=d->more) {
			if(d->resn!=b1->res->id-b2->res->id) continue;
			if(b2->res==b1->res) {
				if(strcmp(d->aim,b1->name)) continue;
				if(strcmp(d->des,b2->name)) continue;
			}
			else {
				if(strcmp(d->aim,b2->name)) continue;
				if(strcmp(d->des,b1->name)) continue;
			}
			if(bound->dist[3*i]-d->dist[1]>1||bound->dist[3*i]-d->dist[0]<-1) {
				cerr<<"error.."<<n1<<" "<<n2<<" "<<endl;
			}
			d->size++;
			break;
		}				

	}

}

void NearDistCut::cut() {

	
	DistPopular *c,*d;

        for(c=dist;c;c=c->next)
        for(d=c;d;d=d->more) {
                d->size=0;
        }


	for(c=dist;c;c=c->next)
        for(d=c;d;d=d->more) {
		if(d->resn==0) d->size=1;
		else if(d->resn==-1) {
			Tatm *b=TRES.findanyatmwithname(d->des);
			if(b==0) continue;
			if(b->id<=4) d->size=1;
		}
        }

}
