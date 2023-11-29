#include "source.h"

ModTop::ModTop() {
  bound=0;
  strfmt=0;
  model=0;
  flag=0;
}

ModTop::~ModTop() {
   if(model&&(flag&TRES.constant->pdbundeletable)==0) delete model;
   if(bound) delete bound;
   if(strfmt) delete strfmt;
}

void ModTop::setmodel() {

	StrFmt *c;//,*c0;

	if(model) {delete model;model=0;}
	if(strfmt==0) {
		cerr<<"the model algnment file is zero"<<endl;
		exit(0);	
	}
	//set the model
	c=strfmt->findfirstsequencefmt();
	
	if(c==0) {
		cerr<<"the model algnment file is zero"<<endl;
                exit(0);
	}
	model=c->pdb;
	flag|=TRES.constant->pdbundeletable;

}

void ModTop::setbound() {
	
	if(bound) {delete bound;bound=0;}
	if(bound==0) bound=new Bound();

	bound->model=model;
	bound->modtop=this;
	bound->flag|=TRES.constant->pdbundeletable;
	bound->ready();	
	bound->setdssp();
	bound->setalldist();
	
	bound->writeatmbound("ie");
} 

void ModTop::addhbondconstraint(HBondList *h) {

}


void ModTop::sethbondconstraint(StrFmt *fmt) {

	int i;
	StrFmt *c;
	HBondList *h;

	StrFmt *sq=fmt->findsequencefmt();

	if(sq==0||sq->seqngap==0) return;

	int n=strlen(sq->seqngap);

	for(i=0;i<n;i++) {
		if(sq->seqngap[i]=='-') continue;
		for(c=fmt;c;c=c->more) {
			if(c->seqngap[i]=='-') continue;
			
			if(c->resn[i]->hbond) {

				for(h=c->resn[i]->hbond;h;h=h->next) {
					addhbondconstraint(h);	
				}
			}
		}
	}
	
}


void ModTop::setbound(StrFmt *fmt) {

	StrFmt *c;

	StrFmt *sq=fmt->findsequencefmt();
	if(sq==0||sq->seqngap==0) return;

	int n=strlen(sq->seqngap);
	
	int i,j;
	int n1,n2;
	int m1,m2;
	Atm *a1,*a2;
	Atm *b1,*b2;
	for(i=0;i<n;i++) {
		for(j=i+1;j<n;j++) {
			if(sq->seqngap[i]=='-'||sq->seqngap[j]=='-') continue;
			m1=sq->match[i];
			m2=sq->match[j];
			//b1=sq->resn[m1]->atm->next;
			//b2=sq->resn[m2]->atm->next;
	     		for(c=fmt;c;c=c->more) {
				if(c->seqngap[i]=='-'||c->seqngap[j]=='-') continue;
				if(c->token&&strcmp(c->token,"sequence")==0) continue;
				n1=c->match[i];
				n2=c->match[j];
				for(int ii=0;ii<=4;ii++) 
				for(int jj=0;jj<=4;jj++) {
					a1=c->resn[n1]->isatmid(ii); 
					a2=c->resn[n2]->isatmid(jj);
					b1=sq->resn[m1]->isatmid(ii);
					b2=sq->resn[m2]->isatmid(jj);
					if(a1==0||a2==0) continue;
					if(abs(b1->res->id-b2->res->id)<=1) continue;
					float d=TRES.distance(a1,a2);
					//cerr<<a1->name<<" "<<a2->name<<" "<<c->code<<d<<endl;
					int ie=b1->tatm->hbond*b2->tatm->hbond;
					float x=(b1->tatm->eng->radius+b2->tatm->eng->radius)*0.9;
					if(ie==2||ie==3||ie==6) x=3.0;
					if(strcmp(b1->tatm->name," SG ")==0&&strcmp(b2->tatm->name," SG ")==0) x=2.0;
					if(b1->id0<b2->id0) bound->addbounds(b1->id0,b2->id0,d,x,d+10,2.0);
					else                bound->addbounds(b2->id0,b1->id0,d,x,d+10,2.0);
				}
	     		}
		}
	}
}
