#include "source.h"

ResBound::ResBound() {
	tres=0;
	next=0;
	distbound=0;
}

ResBound::~ResBound() {
	if(next) delete next;
	if(distbound) delete distbound;	
}


void ResBound::setresidue(Tres *t) {

	if(t==0) return;

	tres=t; 

	distbound=new DistBound();

	Tatm *a;
	
	int nax=-1000;
	for(a=t->tatm;a;a=a->next) {
		nax=max(nax,a->id);
	}
	
	
	distbound->setstartsize(nax*nax+100);
	
	setresdistancebound()
}

void ResBound::setresdistancebound() {

	if(TRES.rotamer==0||TRES.rotamer->next==0) {
		setsimpledistancebound();	
		return;
	}

	Res *r=TRES.rotamer->next->chn->findres(tres->name,0);

	if(r==0) {
		setsimpledistancebound();
	}
	else 
		setcomplexdistancebound(r);
	}
}


void ResBound::setsimpledistancebound() {

	//int nax=-100;
        //for(Tatm *a=tres->tatm;a;a=a->next) nax=max(nax,a->id);

        for(Tatm *a=tres->tatm;a;a=a->next)
        for(Tatm *b=a->next;b;b=b->next) {
		if(a->isnear(b,3)!=-1) continue;
		float e=a->eng->radius+b->eng->radius;
		distbound->addbounds(a->id,b->id,e,e*0.8,1000);	
	}
	return;
}


void ResBound::setcomplexdistancebound(Res *r) {
	
	int nax=-100;
	for(Tatm *a=tres->tatm;a;a=a->next) nax=max(nax,a->id);

	for(Tatm *a=tres->tatm;a;a=a->next) 
	for(Tatm *b=a->next;b;b=b->next) { 
			
		float fmin=1000;
		float fmax=-1000;
		for(Res *r=r->more;r;r=r->more) {
			Atm *aa=r->isatmid(a->id);
			Atm *bb=r->isatmid(b->id);
			float d=TRES.distance(aa->xyz,bb->xyz);
			fmin=min(d,fmin);
			fmax=max(d,fmax); 
		}
		if(fmin==1000||fmax==-1000) continue;
		fmin=fmin-0.05;
		fmax=fmax+0.05;
		distbounds->addbounds(a->id,b->id,(fmin+fmax)/2,fmin,fmax);
	}

	//ca
	for(Tatm *a=tres->tatm;a;a=a->next) {
		float fmin=1000;
        	float fmax=-1000;
        	for(Res *r=r->more;r;r=r->more) {
              		Atm *aa=r->isatmid(a->id);
              		float d=TRES.distance(aa->xyz,r->temp);
              		fmin=min(d,fmin);
              		fmax=max(d,fmax);
        	}
		if(fmin==1000||fmax==-1000) continue;
            	fmin=fmin-0.05;
              	fmax=fmax+0.05;
              	distbounds->addbounds(a->id,nax+1,(fmin+fmax)/2,fmin,fmax);
	}

	//c
	for(Tatm *a=tres->tatm;a;a=a->next) {
                float fmin=1000;
                float fmax=-1000;
                for(Res *r=r->more;r;r=r->more) {
                        Atm *aa=r->isatmid(a->id);
                        float d=TRES.distance(aa->xyz,r->temp+3);
                        fmin=min(d,fmin);
                        fmax=max(d,fmax);
                }
                if(fmin==1000||fmax==-1000) continue;
                fmin=fmin-0.05;
                fmax=fmax+0.05;
                distbounds->addbounds(a->id,nax+2,(fmin+fmax)/2,fmin,fmax);
        }

	//o     
        for(Tatm *a=tres->tatm;a;a=a->next) {
                float fmin=1000;
                float fmax=-1000;
                for(Res *r=r->more;r;r=r->more) {
                        Atm *aa=r->isatmid(a->id);
                        float d=TRES.distance(aa->xyz,r->temp+3);
                        fmin=min(d,fmin);
                        fmax=max(d,fmax);
                }
                if(fmin==1000||fmax==-1000) continue;
                fmin=fmin-0.05;
                fmax=fmax+0.05;
                distbounds->addbounds(a->id,nax+3,(fmin+fmax)/2,fmin,fmax);
        }

	//n
	/*
	for(Tatm *a=tres->tatm;a;a=a->next) {

		float fmin=1000;
                float fmax=-1000;
                for(Res *r=r->more;r;r=r->more) {
		 
                        Atm *aa=r->isatmid(a->id);
                        float d=TRES.distance(aa->xyz,r->temp+3);
                        fmin=min(d,fmin);
                        fmax=max(d,fmax);
                }
                if(fmin==1000||fmax==-1000) continue;
                fmin=fmin-0.05;
                fmax=fmax+0.05;
                distbounds->addbounds(a->id,nax+3,(fmin+fmax)/2,fmin,fmax);
	}
	*/
}


void ResBound::setresidue() {

	Tres *t;

	ResBound *rb=0;
	for(t=&TRES;t;t=t->next) {
		if(rb==0) {
			rb=this;
		}
		else {
			rb->next=new ResBound();
			rb=rb->next;
		}
		rb->setresidue(t);
	}
}
