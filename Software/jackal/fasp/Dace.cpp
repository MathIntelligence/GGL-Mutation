#include"source.h"
Dace::Dace()
{
energyeff=0;
next=0;
flag=0;
pdb=0;
segchn=0;
energy=0;
weight=0;
rmsd=0;
rmsdeff=0;
atoms=0;
more=0;
cid=0;
cutoff=30;
disc=0;
strcpy(force,"udD");
coeff=1;
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
 energyeff=0;
 flag=0;
 weight=0;
 rmsdeff=0;
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
 cid=0; 
 energyeff=0;
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
void Dace::create(float *xyz_all,Res *rr,int nn) {

	float *xyz[2];
	xyz[0]=xyz_all;
	xyz[1]=0;
	create(xyz,rr,nn);
}

void Dace::create(float **xyz_all,Res *rr,int nn) {
	
	if(more) delete more;more=0;

	int kep=0;while(xyz_all[kep])kep++;

	int i;

	Dace *ss=0;
	float *xyz=rr->chn->gettransfer(rr,nn);
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
	rr->chn->transfer(xyz,rr,nn,1);
	ready(rr->chn->pdb,rr->chn->id);	
	create(rr,nn);	
	if(xyz) delete [] xyz;xyz=0;
}

void Dace::calcallenergy() {
	Dace *s;
	if(disc==0) disc=new Disc;
	strcpy(disc->force,force);
	disc->setupall(pdb,8,1);
	Chn *chn=(*pdb)[cid];
	if(chn==0) return;
	Res *rr=chn->isres0(segchn->res->id0);
	if(rr==0) return;
	int nn=segchn->lastres()->id0;
	float *xyzorg=chn->gettransfer(rr,nn);
	disc->setrange();	
	for(s=this;s;s=s->more) {
		chn->transfer(rr,s->segchn->res,nn);
		disc->setup(rr,nn,1); 
		disc->setdipoleascharge(rr,nn); 
		s->energy=disc->clash(rr,nn);		
	}
	chn->transfer(xyzorg,rr,nn,1);
	if(xyzorg) delete [] xyzorg;xyzorg=0;
}

void Dace::update0(float *xyz) {

	//find the chain
	Chn *chn=(*pdb)[cid];	

	if(chn==0||segchn==0) return;

	//copy into the segment
	segchn->transfer(xyz,1);

	//find start and end
	Res *rr=chn->isres0(segchn->res->id0);
	if(rr==0) return;
	int nn=segchn->lastres()->id0;

	//save orginal 
	float *xyzorg=chn->gettransfer(rr,nn);

	//transfer
	chn->transfer(xyz,rr,nn,1);

	//set up disc
	disc->setrange();
	disc->setup(rr,nn,1);
	disc->setdipoleascharge(rr,nn);
 
	//calc energy
	energy=disc->clash(rr,nn);	

	////set back
	//chn->transfer(xyzorg,rr,nn,1);
	//if(xyzorg) delete [] xyzorg;xyzorg=0;
	
	//calc rmsd
	Dace *s;
	for(s=this;s;s=s->more) {
		s->rmsd=getrmsd(this,s);		
	}
	
	//get average rmsd	
	float avgrmsd=0;
	int k=0;
	for(s=this->more;s;s=s->more) {
		avgrmsd+=s->rmsd;
		k++;		
	}
	if(k==0) k=1;
	avgrmsd=avgrmsd/k;
	if(avgrmsd==0) avgrmsd=1;

	//set up rmsd eff
	avgrmsd=1/avgrmsd/avgrmsd;
	float xx=-log(0.5);
	
 	
	for(s=this;s;s=s->more) {
		s->rmsdeff=s->rmsd*s->rmsd*avgrmsd*xx;
	}

	
	
	float ee=0.3*8.31;

	//find minimal energy
	float emin=energy;
	for(s=this;s;s=s->more) {	
		if(s->energy<emin) {
			emin=s->energy;
		}
	}
			
	if(emin>-cutoff/2) {
		emin=emin+cutoff/2;	
		for(s=this;s;s=s->more) {			 
			s->weight=s->energy-emin;		 
		}	
 	}
	else {
		for(s=this;s;s=s->more) {			 
			s->weight=s->energy;		 
		}
	}

	emin=weight;
	for(s=this;s;s=s->more) {	
		if(s->energy<emin) {
			emin=s->weight;
		}
	}
		
	//get average energy
	float  eavg=0;
 
	k=0;
	for(s=this;s;s=s->more) {	
		if(s->weight-emin>20) continue;
		eavg+=s->weight;
		k++;
	}
	if(k==0) k=1;
	eavg=eavg/k;
			 
	//calc coefficient of the energy
	float eff=1;
	if(eavg-emin>cutoff) {
		eff=cutoff/(eavg-emin);
	}	
		
	//calc weight
	for(s=this;s;s=s->more) {	
		float t=-eff*s->weight/ee-s->rmsdeff;
		//if(t>0) t=0;
		s->weight=exp(t)/ee;
	}

	//set back
	chn->transfer(xyzorg,rr,nn,1);
	if(xyzorg) delete [] xyzorg;xyzorg=0;
}


void Dace::update1(float *xyz) {

	//find the chain
	Chn *chn=(*pdb)[cid];	

	if(chn==0||segchn==0) return;

	//copy into the segment
	segchn->transfer(xyz,1);

	//find start and end
	Res *rr=chn->isres0(segchn->res->id0);
	if(rr==0) return;
	int nn=segchn->lastres()->id0;

	//save orginal 
	float *xyzorg=chn->gettransfer(rr,nn);

	//transfer
	chn->transfer(xyz,rr,nn,1);

	//set up disc
	disc->setrange();
	disc->setup(rr,nn,1);
	disc->setdipoleascharge(rr,nn);
 	more->energy=-10; //test....
	//calc energy
	energy=disc->clash(rr,nn);	
	cerr<<"energy:"<<energy<<endl;
	//calc rmsd
	Dace *s;
	for(s=this;s;s=s->more) {
		s->rmsd=getrmsd(this,s);
		cerr<<"rmsd:"<<s->rmsd<<endl;		
	}
	
	//get average rmsd	
	float avgrmsd=0;

	int k=0;
	for(s=this;s;s=s->more) {
		avgrmsd+=s->rmsd;
		k++;		
	}
	if(k==0) k=1;
	avgrmsd=avgrmsd/k;
	if(avgrmsd==0) avgrmsd=1;

	//set up rmsd eff
	avgrmsd=1/avgrmsd/avgrmsd;
	float xx=-log(0.5);
			
	for(s=this;s;s=s->more) {
		s->rmsdeff=s->rmsd*s->rmsd*avgrmsd*xx;
	}	
	
	float ee=0.3*8.31;

	for(s=this;s;s=s->more) {
		float t=-s->energy/ee-s->rmsdeff+energy/ee;	 
		s->energyeff=t;
	}

	float emax=energyeff;
	for(s=this;s;s=s->more) {
		if(s->energyeff>emax) 	emax=s->energyeff;
	}
 
	for(s=this;s;s=s->more) {
		s->energyeff=s->energyeff-emax;
	}

	float d=0;
	for(s=this;s;s=s->more) {
		d+=exp(s->energyeff);		
	}
	energy=energy-emax*ee-log(d)*ee;
	
	for(s=this;s;s=s->more) {
		s->weight=ee*exp(s->energyeff)/d; 
	}
 	
	//set back
	chn->transfer(xyzorg,rr,nn,1);
	if(xyzorg) delete [] xyzorg;xyzorg=0;
}

void Dace::update(float *xyz) {

	//find the chain
	Chn *chn=(*pdb)[cid];	

	if(chn==0||segchn==0) return;

	//copy into the segment
	segchn->transfer(xyz,1);

	//find start and end
	Res *rr=chn->isres0(segchn->res->id0);
	if(rr==0) return;
	int nn=segchn->lastres()->id0;

	//save orginal 
	float *xyzorg=chn->gettransfer(rr,nn);

	//transfer
	chn->transfer(xyz,rr,nn,1);

	//set up disc
	disc->setrange();
	disc->setup(rr,nn,1);
	disc->setdipoleascharge(rr,nn);
  
	//calc energy
	energy=disc->clash(rr,nn);	
	cerr<<"energy:"<<energy<<endl;
	//calc rmsd
	Dace *s;
	for(s=this;s;s=s->more) {
		s->rmsd=getrmsd(this,s);
		cerr<<"rmsd:"<<s->rmsd<<endl;		
	}
	
	//get average rmsd	
	float avgrmsd=0;

	int k=0;
	for(s=this;s;s=s->more) {
		avgrmsd+=s->rmsd;
		k++;		
	}
	if(k==0) k=1;
	avgrmsd=avgrmsd/k;
	if(avgrmsd==0) avgrmsd=1;

	//set up rmsd eff
	avgrmsd=1/avgrmsd/avgrmsd;
	float xx=-log(0.5);
			
	for(s=this;s;s=s->more) {
		s->rmsdeff=s->rmsd*s->rmsd*avgrmsd*xx;
		s->coeff=avgrmsd*xx;
	}	
	
	float eav=0;

	k=0;
	for(s=this->more;s;s=s->more) {
		eav+=fabs(s->energy-energy);
		k++;
	}

	if(k==0) k=1;
	eav=eav/k;

	float ee=0.3*8.31;

	float d=0;
	if(eav>10) d=10/eav; 	
	else d=1;
	for(s=this;s;s=s->more) {
		s->energyeff=-d*s->energy/ee-s->rmsdeff;
	}
 	coeff=d;	
 	
	d=0;
	for(s=this;s;s=s->more) {
		d+=exp(s->energyeff);		
	}

	for(s=this;s;s=s->more) {
		if(s==this) s->weight=s->coeff*exp(s->energyeff)/d; 
		else        s->weight=ee*s->coeff*exp(s->energyeff)/d; 
	}
 	
	/*
	//calculate energy
	for(s=this;s;s=s->more) {
		s->energyeff=-s->energy/ee*coeff-s->rmsdeff;
	}
	
	float emax=energyeff;
	for(s=this;s;s=s->more) {
		if(s->energyeff>emax) 	emax=s->energyeff;
	}
 
	for(s=this;s;s=s->more) {
		s->energyeff=s->energyeff-emax;
	}	

	d=0;
	for(s=this;s;s=s->more) {
		d+=exp(s->energyeff);		
	}

	energy=-ee*emax-ee*log(d);
	*/

	//set back
	chn->transfer(xyzorg,rr,nn,1);
	if(xyzorg) delete [] xyzorg;xyzorg=0;
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
		if(a->tatm->id>4) continue;
		Atm *b=bb->atoms[a->id0];
		if(b==0) continue;
		if(b->tatm->id>4) continue;
		if(a->tatm->id!=b->tatm->id) continue;
		if(a->res->name!=b->res->name) continue;
		d+=TRES.distsqr(a->xyz,b->xyz);
		k++;
	}
	if(k==0) return 0;
	return sqrt(d/k);
}
