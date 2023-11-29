#include"source.h"
Dace::Dace()
{
next=0;
flag=0;
pdb=0;
segchn=0;
energy=0;
weight=0;
rmsd=0;
atoms=0;
more=0;
cid=0;
cutoff=20;
}

Dace::~Dace()
{
 if(disc) delete disc;disc=0;
 if(next) delete next;next=0;
 if(atoms) delete [] atoms;atoms=0;
 if(segchn) delete segchn;segchn=0;
 if(next) delete next;next=0;
 if(more) delete more;more=0;
 rmsd=0;
 energy=0;
 pdb=0;
 flag=0;
 weight=0;
}

void Dace::clear() {
 if(atoms) delete [] atoms;atoms=0;
 if(segchn) delete segchn;segchn=0;
 if(disc) delete disc;disc=0;
 rmsd=0;
 energy=0;
 pdb=0;
 flag=0;
 weight=0;
}

void Dace::ready(Pdb *s,int c0) {
   	clear();
	pdb=s;
	if(pdb==0) return;
	int n=pdb->maxatmid0()+100;
	int i;
	atoms=new Atm*[n];
	for(i=0;i<n;i++)atoms[i]=0;
	cid=c0;
}

void Dace::create(Res *rr,int nn) {

	if(segchn) delete segchn;segchn=0;
	segchn=new Chn;
	segchn->create(rr,nn);
	Res *r;
	Atm *a;
	for(r=segchn->res;r;r=r->next) 
	for(a=r->atm;a;a=a->next) {
		atoms[a->id0]=a;
	}	
} 

void Dace::create(float **xyz_all,Res *rr,int nn) {
	
	if(more) delete more;more=0;

	int kep=0;while(xyz_all[kep])kep++;

	int i;

	Dace *ss=0;

	for(i=0;i<kep;i++) {
		if(ss==0) {
			more=new Dace;
			ss=more;
		}		
		else {
			ss->more=new Dace;
			ss=ss->more;
		}
		ss->ready(rr->chn->pdb,rr->chn->id);		
		ss->create(rr,nn);
		ss->segchn->transfer(xyz_all[i],1);
		Res *r;
		Atm *a;
		for(r=ss->segchn->res;r;r=r->next) 
		for(a=r->atm;a;a=a->next) {
			ss->atoms[a->id0]=a;
		}	
	}

	ss->create(rr,nn);	
}

void Dace::calcallenergy() {
	Dace *s;
	if(disc==0) disc=new Disc;
	strcpy(disc->force,"uDd");
	disc->setupall(pdb,8,1);
	Chn *chn=(*pdb)[cid];
	if(chn==0) return;
	Res *rr=chn->isres0(segchn->res->id0);
	if(rr==0) return;
	int nn=segchn->lastres()->id0;
	float *xyzorg=chn->gettransfer(rr,nn);
	for(s=this;s;s=s->more) {
		chn->transfer(rr,s->segchn->res,nn);
		s->energy=disc->clash(rr,nn);		
	}
	chn->transfer(xyzorg,rr,nn,1);
	if(xyzorg) delete [] xyzorg;xyzorg=0;
}

void Dace::update(float *xyz) {

	Chn *chn=(*pdb)[cid];	
	if(chn==0||segchn==0) return;
	segchn->transfer(xyz,1);
	Res *rr=chn->isres0(segchn->res->id0);
	if(rr==0) return;
	int nn=segchn->lastres()->id0;
	float *xyzorg=chn->gettransfer(rr,nn);
	chn->transfer(rr,segchn->res,nn);
	energy=disc->clash(rr,nn);	
	chn->transfer(xyzorg,rr,nn,1);
	if(xyzorg) delete [] xyzorg;xyzorg=0;

	Dace *s;
	for(s=this;s;s=s->more) {
		s->rmsd=getrmsd(this,s);		
	}
	
	float avgrmsd=0;
	int k=0;
	for(s=this->more;s;s=s->more) {
		avgrmsd+=s->rmsd;
		k++;		
	}
	if(k==0) k=1;
	avgrmsd=avgrmsd/k;
	if(avgrmsd==0) avgrmsd=1;
	avgrmsd=1/avgrmsd/avgrmsd;
	float xx=-log(0.5);
	float ee=0.3*8.31;
 

	
	for(s=this;s;s=s->more) {
		s->rmsd=s->rmsd*s->rmsd*avgrmsd*xx;
	}
	
	float emin=energy;
	for(s=this;s;s=s->more) {	
		if(s->energy>emin) {
			emin=s->energy;
		}
	}

	float  eavg=0;
	k=0;
	for(s=this;s;s=s->more) {	
		if(s->energy-emin>20) continue;
		eavg+=s->energy;
		k++;
	}
	if(k==0) k=1;
	eavg=eavg/k;
	float eff=1;
	if(eavg-emin>cutoff) {
		eff=cutoff/(eavg-emin);
	}	

	for(s=this;s;s=s->more) {	
		float t=eff*s->energy/ee-s->rmsd*s->rmsd*avgrmsd*xx;
		s->weight=exp(-t);
	}
}


float Dace::getrmsd(Dace *aa,Dace *bb) {

	if(aa==bb) return 0;
	if(aa==0||bb==0) return 0;
	Res *r;
	Atm *a;

	float d=0;
	int k=0;
	for(r=aa->segchn->res;r;r=r->next)
	for(a=r->atm;a;a=a->next) {
		if(a->id0>4) continue;
		Atm *b=bb->atoms[a->id0];
		if(b==0) continue;
		if(a->res->name!=b->res->name) continue;
		d+=TRES.distsqr(a->xyz,b->xyz);
		k++;
	}
	if(k==0) return 0;
	return sqrt(d/k);
}
