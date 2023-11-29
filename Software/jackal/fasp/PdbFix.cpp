#include"source.h"
PdbFix::PdbFix()
{
next=0;
pdborg=0;
pdb=0; 
exact=1; //fixed atoms exact as original
atmat=0;
onlybackbone=0;
fast=3;
ctrip=0;
}

PdbFix::~PdbFix()
{
 if(next) delete next;next=0;
 if(pdborg) delete pdborg;pdborg=0;
 if(pdb)  delete pdb;pdb=0; 
 if(atmat) delete [] atmat;atmat=0;
}

void PdbFix::setatmat(Pdb *s) {

	if(atmat) delete [] atmat;atmat=0;
	int n=s->maxresid0()+100;
	atmat=new Res*[n];
	int i;
	for(i=0;i<n;i++) atmat[i]=0;

	Chn *c,*c0;
	for(c=s->chn;c;c=c->next){
		c0=pdborg->ischain(c->id);
		if(c0==0) continue;
		Res *r,*r0;
		for(r=c->res;r;r=r->next) {
			r0=c0->isresoid(r->oid);
			if(r0==0) continue;
			atmat[r->id0]=r0;
		}
	}	
}

void PdbFix::addhnatoms(Pdb *s){

	Res *r;
	Chn *c;
	for(c=s->chn;c;c=c->next)
	for(r=c->res;r;r=r->next) {
		if(r->temp==0) continue;
		r->addhnatoms(1);
	}
	s->configure();
}

void PdbFix::addsidechains(Pdb *s){
		
	Res *r;
	Chn *c;	
	//
	for(c=s->chn;c;c=c->next)
	for(r=c->res;r;r=r->next) {
		if(r->temp==0) continue; //backbone not complete
		Tatm *a;
		int all=1;
		for(a=r->tres->tatm;a;a=a->next) {
			if(a->name[1]=='H') continue;
			Atm *b=r->isatmid(a->id);
			if(b==0) {all=0;break;}
		}
		if(all==1) {
			r->addhatoms(1);
			continue; //all heavy atom there!
		}
		if(r->name!='G') {
			r->addsidechain();
			pdb->configure();			
			pdb->setflgr(-99999);
			r->flag=10000;
			Scap scprd;
        		scprd.includeself=1;
        		scprd.singletorsion=1;
        		scprd.colonyline=1;
        		scprd.colony=2;
        		scprd.ncolony=1;
        		scprd.nncolony=1;
        		scprd.nummore=1000;
        		scprd.bmax=2;
        		scprd.tormax=2;
        		scprd.ring=1;
        		scprd.pdb=s;
        		scprd.mutate();
			//scprd.pdb=0;					
			Res *rm;
			Res *r0=atmat[r->id0];
			float rmax=10000;
			float txyz[200];
			for(rm=r;rm;rm=rm->more) {
				float d=getsidechainrmsd(r0,rm);
				if(d<rmax) {
					rmax=d;
					rm->transfer(txyz,0);
				}		
			}
			r->transfer(txyz,1);
			r->setfreenextonmore();			
		}
		r->addhatoms(1);	
	}
	s->configure();
}

void PdbFix::addsidechains(Res *r1,Res *r2){
		
	Res *r;
 
	//	 
	for(r=r1;r;r=r->next) {
		if(r->id0>r2->id0) break;
		 
		if(r->name!='G') {
			r->addsidechain();
			pdb->configure();			
			pdb->setflgr(-99999);
			r->flag=10000;
			Scap scprd;
        		scprd.includeself=1;
        		scprd.singletorsion=1;
        		scprd.colonyline=1;
        		scprd.colony=2;
        		scprd.ncolony=1;
        		scprd.nncolony=1;
        		scprd.nummore=1000;
        		scprd.bmax=2;
        		scprd.tormax=2;
        		scprd.ring=1;
        		scprd.pdb=r1->chn->pdb;
        		scprd.mutate();
			//scprd.pdb=0;					
			Res *rm;
			Res *r0=atmat[r->id0];
			float rmax=10000;
			float txyz[200];
			for(rm=r;rm;rm=rm->more) {
				float d=getsidechainrmsd(r0,rm);
				if(d<rmax) {
					rmax=d;
					rm->transfer(txyz,0);
				}		
			}
			r->transfer(txyz,1);
			r->setfreenextonmore();			
		}
		r->addhatoms(1);	
	}	
}

void PdbFix::sidechainminimal(Pdb *s){
	Scap scprd;
        scprd.singletorsion=1;
        scprd.colonyline=1;
        scprd.colony=2;
        scprd.ncolony=1;
        scprd.nncolony=1;
        scprd.nummore=100;
        scprd.bmax=2;
        scprd.tormax=2;
        scprd.ring=1;
        scprd.pdb=s;
        scprd.includeself=1;
        scprd.resultout=0;
	
	Res *r;
 	Chn *c;
	//
	int ntot=0;
	pdb->setflgr(-99999);		
	for(c=s->chn;c;c=c->next)
	for(r=c->res;r;r=r->next) {
		if(strchr("AGP",r->name)) continue;		
		Res *r0=atmat[r->id0];
		if(r0==0) {
			r->flag=10000;
			ntot++;
			continue;
		}
			
		Atm *a;
		int jj=0;
		for(a=r0->atm;a;a=a->next) {
			if(a->tatm->name[1]=='H') continue;
			if(a->tatm->id>4) {
				jj=1;
				break;
			}
		}
		if(jj==0) {
			r->flag=10000;
			ntot++;
			continue;
		}
	}	
	if(ntot) {
		scprd.scpred(1);
	}
	scprd.pdb=0;
	for(c=s->chn;c;c=c->next) c->setfreenextonmore();
	if(TRES.logg>3) s->write("side");
}
float PdbFix::getsidechainrmsd(Res *r0,Res *r){

	Atm *aa,*bb;
	float d=0;
	int n=0;
	for(aa=r0->atm;aa;aa=aa->next) {
		if(aa->tatm->name[1]!='H'&&aa->tatm->id<=4) continue;
		if(aa->tatm->name[1]=='H'&&aa->tatm->bond[0]->id<=4) continue;
		bb=r->isatmid(aa->tatm->id);
		if(bb==0) continue;
		d+=TRES.distsqr(aa->xyz,bb->xyz);
		n++;
	}		
	if(n==0) n=1;
	d=sqrt(d/n);
	 	
	return d;
}

void PdbFix::tuneheader(Pdb *s) {
	
	Chn *c;
	//check header file;
	for(c=pdb->chn;c;c=c->next) { 
		Res *r0=c->lastres();
		Res *r;
		for(r=r0;r;r=r->last) {
			if(r->temp==0) continue;
			if(r->last&&r->id-r->last->id==1) {
				if(r->last->temp==0) {
					delete [] r->temp;r->temp=0;
					r->nemp=0;
				}
			}
			else if(r->last&&r->id-r->last->id>1) {
				delete [] r->temp;r->temp=0;
				r->nemp=0;
			}	
		}
	}

}



int PdbFix::checkpdb(Pdb *s){

	Chn *c;
	Res *r;
	int n=0;
	cerr<<"checking missing atoms..."<<endl;
	for(c=s->chn;c;c=c->next)
	for(r=c->res;r;r=r->next) {
		Tatm *a;
		for(a=r->tres->tatm;a;a=a->next){
			Atm *b=r->isatmid(a->id);
			if(b==0) {
				n++;
				cerr<<"Atom missing: "<<r->name<<r->oid<<"--"<<a->name<<"!"<<endl;
			}	
		}
	}
	cerr<<"the total number of missing atoms: "<<n<<endl;
	if(n==0) {
		cerr<<"there is no missing atoms. "<<endl;
		return n;
	}
	else {
		cerr<<"add back missing atoms..."<<endl;
		return n;
	}
}
void PdbFix::fixlostresidues(Pdb *s) {
	fixlostresidues(s,10000);
}
void PdbFix::fixlostresidues(Pdb *s,int lent) {
	
	StrFmt *nest=new StrFmt();
	StrFmt *sc=0;
	Chn *c;
	int ntot=0;
	for(c=s->chn;c;c=c->next) {
		s->readseqcard(c->id);
		char *seq=c->seqcard;
		c->seqcard=0;
		if(seq==0) {
			cerr<<endl;
			cerr<<"no sequence card found for chain";
			if(c->id!=' ') cerr<<": "<<c->id<<endl;
			else	cerr<<endl;
			cerr<<"no residues added"<<endl;							 
		}
				
		char *seqn=c->getseqn(c->res,c->lastres()->id0+1000);
		if(seqn==0||strlen(seqn)==0) {
			cerr<<endl;
			cerr<<"no atom coordinates exist in chain";
			if(c->id!=' ') cerr<<": "<<c->id<<endl;
			else	cerr<<endl;							 
		}

		if(seq==0) seq=c->getseqn(c->res,c->lastres()->id0+1000);
		if(seqn) delete [] seqn; seqn=0;
		if(seq) {
			seqn=strdup(seq);
		}	
		if(sc==0) {
			sc=nest;
		}
		else {
			sc->next=new StrFmt();
			sc=sc->next;
		}		
		sc->token=strdup("structure");
		sc->cid=c->id;
		sc->seqngap=seqn;		
		sc->pdb=new Pdb(pdb);
		sc->pdb->configure();
		sc->code=strdup("unknown");
		sc->more=new StrFmt();
		sc->more->token=strdup("sequence");
		sc->more->seqngap=seq;
		char test[100];
		sprintf(test,"test%i",ntot);
		sc->more->code=strdup(test);
		ntot++;
	}
	for(sc=nest;sc;sc=sc->next) {
		StrFmt *sc0;
		sc->profix=1;
		for(sc0=sc->more;sc0;sc0=sc0->more) {
			sc0->profix=1;
		}
        }
	if(ntot) {
		int out=0;
		int tune=0;
		//int logg=1;
		int fapr=3;
		int seed=18120;
		int sharp=0;
		int rmsd=(int)(2.0);
		int asloop=0;
		int test=0;
		int seglen=2;
		//int nopt=0;
		int restraint=1;
		nest->tune=tune;
		nest->removeallnonstandard();
		nest->clearemptyseq();
		nest->alnback=0;
		nest->prepare();	
		nest->setmutate("fapr",fapr);
		nest->setmutate("seed",seed);
		nest->setmutate("test",test);
		nest->setmutate("sharp",sharp);
		nest->setmutate("out",out);
		nest->setmutate("rmsd",rmsd*100);
		nest->setmutate("asloop",asloop);
		nest->setmutate("seglen",seglen);
		nest->setmutate("restraint",restraint);
		nest->exec();
	}	
	
	delete pdb->chn;
	pdb->chn=0;
	c=0;
	

	for(sc=nest;sc;sc=sc->next) {
		if(c==0) {
			pdb->chn=sc->mutate->mpdb->chn;
			sc->mutate->mpdb->chn=0;
			pdb->chn->pdb=pdb;
			c=pdb->chn;
		}
		else {
			c->next=sc->mutate->mpdb->chn;
			sc->mutate->mpdb->chn=0;
			c->next->pdb=pdb;
			c=c->next;
		}
	}
	delete nest;
}

void PdbFix::ready() {
	char **info;
        char **endinfo;
	Chn *c,*c1;
	//
	info=pdb->info;
        endinfo=pdb->endinfo;
	pdb->info=0;
	pdb->endinfo=0;
	pdborg=new Pdb(pdb);
	pdb->configure();
	pdborg->configure();
	setatmat(pdb); 
	pdb->headerpdbfix();
	//tuneheader(pdb);		
        Res *re;	
        for(c=pdb->chn;c;c=c->next) 
        for(re=c->res;re;re=re->next) re->addhanyway=1;
        //
	if(onlybackbone==0) addsidechains(pdb);
	if(TRES.logg>3) pdb->write("temp1.pdb");
	//
	delete pdborg;pdborg=0;
	pdborg=new Pdb(pdb);
	pdborg->configure();
	pdborg->headerpdbfix();
	pdborg->info=info;
	pdborg->endinfo=endinfo;
	//tuneheader(pdborg);
	if(TRES.logg>3) pdborg->write("temp.pdb");
	pdb->configure();
	if(TRES.logg>3) pdb->write("temp.pdb"); 
	//
	c1=0;
	for(c=pdb->chn;c;c=c->next) {  
		//
		Chn *c0=createchain(c);	
			
		c0->setallnear(c0->res,10000);
		c0->pdb=pdb;
		//
		c0->next=c->next;
		c->next=0;		
		delete c;c=0;
		if(c1==0) pdb->chn=c0;
		else 	  c1->next=c0;
		 	
		c=c0;
		//	 		 
		c1=c; 
	}
        if(TRES.logg>3) pdb->write("temp.pdb");
	pdb->configure();
	if(TRES.logg>3) pdb->write("temp.pdb");
	for(c=pdb->chn;c;c=c->next) {
		c1=pdborg->ischain(c->id);
		c->start=c1->start;
	}
	setatmat();
}
 
Chn *PdbFix::createchain(Chn *c) {

	Rotate rot;
	//create a new chain with the same sequecen as the old

	char *seqn=c->getseqn();
	Chn *c0=new Chn();
	c0->create(seqn);
	
	if(TRES.logg>3) c0->write("s");

	//configure them

	c0->configure();
	c->configure();
	c0->id=c->id;
	c0->start=c->start;
		
	//assign old id to new chains
	Res *r,*r0;	 
	for(r=c->res;r;r=r->next) {
		r0=c0->isres0(r->id0);
		if(r0==0||r0->name!=r->name) {
			cerr<<"program detected errors when trying to fix missing atoms..."<<endl;
			exit(0);
		}		 
		r0->oid=r->oid; //old id
		r0->id=r->id;   //new id
		strcpy(r0->rescard,r->rescard);
	}
	return c0;
}

void PdbFix::myfix() {

        Chn *c;
 	setatmat();
	pdb->transform(8);//backbone including HN
        pdb->configure();
        pdborg->setlastresorder();

        for(c=pdb->chn;c;c=c->next) { 	
                c->clearhbond();
                mynewfix(c);              	
                sethbonddssp(c);
                mynewfix(c);
                minimize(c);	
        }

        if(TRES.logg>3) pdb->write("ss1");
}

void PdbFix::pushaway(Pdb *s) {
	float coo[3];
	pdborg->getmaxxyz(coo);
	int i;
	for(i=0;i<3;i++) coo[i]=100*coo[i]+10000;
	s->addxyz(coo);
}

int PdbFix::isbackboneatom(Res *first,Res *last){

	Res *r;
	
	for(r=first;r;r=r->next) {
		if(r->id0>last->id0) break;
		Atm *aa;
		int jj=0;
		for(aa=r->atm;aa;aa=aa->next){
			if(aa->tatm->name[1]=='H') continue;
                        if(aa->tatm->id>4) continue;                        	
                        jj++;
                }
		if(jj==0) return 0;
	}
	return 1;
}

void PdbFix::myfixnow() {
	Rotate rot;
        Chn *c;
        pdborg->setlastresorder();
	pdb->setlastresorder();
	//pushaway(pdb);
	//start
	//copy every reliable segments
	for(c=pdb->chn;c;c=c->next) {
		Res *r;
		c->clearhbond();
		for(r=c->lastres();r;r=r->last) {
			Res *r0=atmat[r->id0];
			if(r0->temp) {
				r->transfer(r0);
				r->transfertemp(r0);	
				continue;
			}
		}			
	}
	//end
		 
        for(c=pdb->chn;c;c=c->next) { 	
		Res *r;
		c->clearhbond();
		for(r=c->lastres();r;r=r->last) {
			Res *r0=atmat[r->id0];
			if(r0->temp) continue;
			
			Res *rend;
			for(rend=r;rend;rend=rend->last) {
				r0=atmat[rend->id0];				
				if(r0->temp) {
					rend=rend->next;
					break;
				}
				if(rend->last&&rend->id-rend->last->id!=1) break;
			} 		
			if(rend==0) rend=c->res;  				
			Res *t=r;
			//
			Res *fr0,*fr1,*fr2;
			//fr0
			fr2=c->res;
			if(rend->last&&rend->last->id-rend->id==-1&&atmat[rend->last->id0]->temp) {
				fr0=rend->last->last;
				rend->last->last=0;
				c->res=rend->last;	
			}
			else {
				fr0=rend->last;
				rend->last=0;
				c->res=rend;
			}
			//fr1
			if(t->next&&t->next->id-t->id==1&&atmat[t->next->id0]->temp) {
				fr1=t->next->next;
				t->next->next=0;
			}
			else {
				fr1=t->next;
				t->next=0;
			}
			int bc=isbackboneatom(rend,t);	
			int start=c->res->id;int start0=c->res->id0;		
			if(t->id0-rend->id0+1<5&&(rend->last||t->next)&&bc==0) {
				if(rend->last)  rot.link(rend->last,rend,rend->id0,1);
				if(t->next) rot.link(t,t->next,t->id0,0);
				
				int n=directlink(rend,t,1); 
		
				if(n==0&&t->next) mynewfix(rend,t);
				else if(n==0&&rend->last)  mydirectnewfix(rend,t);
				else if(n==0) mynewfix(rend,t);
			}
			else { 
				if(rend->last)  rot.link(rend->last,rend,rend->id0,1);
				if(t->next) rot.link(t,t->next,t->id0,0);
				if(t->next) mynewfix(rend,t);
				else if(rend->last)  mydirectnewfix(rend,t);
				else mynewfix(rend,t);                			
			}
			Res *re;
			for(re=c->res;re;re=re->next) {
				re->id=start;re->id=start0;
				start++;start0++;
			}
			if(rend->last&&rend->last->id-rend->id==-1&&atmat[rend->last->id0]->temp) {
				rend->last->last=fr0;
				c->res=fr2;	
			}
			else {
				rend->last=fr0;
				c->res=fr2;
			}
			//fr1
			if(t->next&&t->next->id-t->id==1&&atmat[t->next->id0]->temp) {					
				t->next->next=fr1;
			}
			else {
				t->next=fr1;
			}		
			//rend->last=fr0;
			//t->next=fr1;
			pdb->configure();
			if(TRES.logg>3) pdb->write("sa");
			//r=rend->last;  
			r=rend; 
			if(r==0) break;  
		}
		if(TRES.logg>3) pdb->writeold("ss0");
		sethbonddssp(c);
		for(r=c->lastres();r;r=r->last) {
			Res *r0=atmat[r->id0];
			if(r0->temp) continue;
			
			Res *rend;
			for(rend=r;rend;rend=rend->last) {
				r0=atmat[rend->id0];				
				if(r0->temp) {
					rend=rend->next;
					break;
				}
				if(rend->last&&rend->id-rend->last->id!=1) break;
			} 		
			if(rend==0) rend=c->res;  				
			Res *t=r;
			//
			Res *fr0,*fr1,*fr2;
			//fr0
			fr2=c->res;
			if(rend->last&&rend->last->id-rend->id==-1&&atmat[rend->last->id0]->temp) {
				fr0=rend->last->last;
				rend->last->last=0;
				c->res=rend->last;	
			}
			else {
				fr0=rend->last;
				rend->last=0;
				c->res=rend;
			}
			//fr1
			if(t->next&&t->next->id-t->id==1&&atmat[t->next->id0]->temp) {
				fr1=t->next->next;
				t->next->next=0;
			}
			else {
				fr1=t->next;
				t->next=0;
			}
			int bc=isbackboneatom(rend,t);
			int start=c->res->id;
			int start0=c->res->id0;			
			if(t->id0-rend->id0+1<5&&(rend->last||t->next)&&bc==0) {
				if(rend->last)  rot.link(rend->last,rend,rend->id0,1);
				if(t->next) rot.link(t,t->next,t->id0,0);
				
				int n=directlink(rend,t,0); 
		
				if(n==0&&t->next) mynewfix(rend,t);
				else if(n==0&&rend->last)  mydirectnewfix(rend,t);
				else if(n==0) mynewfix(rend,t);
			}
			else { 
				if(rend->last)  rot.link(rend->last,rend,rend->id0,1);
				if(t->next) rot.link(t,t->next,t->id0,0);
				if(t->next) mynewfix(rend,t);
				else if(rend->last)  mydirectnewfix(rend,t);
				else mynewfix(rend,t);                			
			}
			if(TRES.logg>3) pdb->write("s1");
			linkallres(c->res,c->lastres());
			if(TRES.logg>3) pdb->write("s2");
			minimize(rend,t);
			if(0) {
				int ei=0;
				int nn;
				if(ctrip==1) nn=(5-fast)*2;
				else         nn=(5-fast);
				while(ei<nn) {
					shiftcenter(rend,t);					 
					float d=energy(rend,t);					
					if(d<0.05) {						 
						break;
					}
					minimize(rend,t);
					ei++;
				} 
			}

			if(0) {
				int ei=0;
				int nn;
				if(ctrip==1) nn=(5-fast)*2;
				else         nn=(5-fast);
				while(ei<nn) {
					shiftcenter(rend,t);
					Res **stem=maxenergystem(rend,t);
					float d=energy(stem[0],stem[1]);					
					if(d<0.05) {
						if(stem) delete [] stem;stem=0;
						break;
					}
					minimize(stem[0],stem[1]);
					if(stem) delete [] stem;stem=0;
					ei++;
				} 
			}
			if(TRES.logg>3) pdb->write("s3");	
			Res *re;
			for(re=c->res;re;re=re->next) {
				re->id=start;
				re->id0=start0;
				start++;start0++;
			}
			if(rend->last&&rend->last->id-rend->id==-1&&atmat[rend->last->id0]->temp) {
				rend->last->last=fr0;
				c->res=fr2;	
			}
			else {
				rend->last=fr0;
				c->res=fr2;
			}
			//fr1
			if(t->next&&t->next->id-t->id==1&&atmat[t->next->id0]->temp) {					
				t->next->next=fr1;
			}
			else {
				t->next=fr1;
			}		
			//rend->last=fr0;
			//t->next=fr1;
			pdb->configure();
			if(onlybackbone==0) addsidechains(rend,t);
			if(TRES.logg>3) pdb->write("sa");
			//r=rend->last;  
			r=rend; 
			if(r==0) break;  
		}
        }
	 
	if(onlybackbone==0) sidechainminimal(pdb); 
	else pdb->transform(8);
        if(TRES.logg>3) pdb->writeold("ss1");
	pdb->configure();
}

void PdbFix::shiftcenter(Res *first,Res *end) {

	Res *r;

	for(r=first;r;r=r->next) {

		if(r->id0>end->id0) break;

		Res *r0=atmat[r->id0];
		if(r0==0) continue;
	
		Atm *a0,*a;
		float coo[3];
		float coo0[3];
		int i;
		for(i=0;i<3;i++) coo[i]=0;
		for(i=0;i<3;i++) coo0[i]=0;
		int ii=0;
		for(a0=r0->atm;a0;a0=a0->next) {		
			if(a0->tatm->name[1]=='H'&&a0->tatm->bond[0]->id>3) continue;
			if(a0->tatm->name[1]!='H'&&a0->tatm->id>4) continue;
			a=r->isatmid(a0->tatm->id);
			if(a==0) continue;
			for(i=0;i<3;i++) coo[i]+=a->xyz[i];
			for(i=0;i<3;i++) coo0[i]+=a0->xyz[i];
			ii++;	
		}
		if(ii==0) continue;
		for(i=0;i<3;i++) coo[i]=coo[i]/ii;
		for(i=0;i<3;i++) coo0[i]=coo0[i]/ii;
		 
		float d=TRES.distance(coo,coo0);
		if(TRES.logg>3) cerr<<r->name<<r->oid<<" :" <<d<<endl; 
		 
		if(d<0.1) {
			float dx=(coo0[0]-coo[0]);
			float dy=(coo0[1]-coo[1]);
			float dz=(coo0[2]-coo[2]);
			Atm *aa;
			for(aa=r->atm;aa;aa=aa->next) {
				aa->xyz[0]+=dx;
                                aa->xyz[1]+=dy;
                                aa->xyz[2]+=dz;
			}
			int i;
			for(i=0;i<3;i++) {
				r->temp[i*3+0]+=dx;
				r->temp[i*3+1]+=dy;
				r->temp[i*3+2]+=dz;
			}		 
		}	  
	}
}

void PdbFix::sethbonddssp(Chn *chn) {
        chn->clearhbond();
        chn->header();
        chn->buildhbond();
        chn->setdsspstr();
        chn->buildssbond();
        chn->setthreestatesec();
        if(TRES.logg>3) chn->writesecondary(stderr);     	
}


void PdbFix::setatmat() {

	if(atmat) delete [] atmat;atmat=0;
	int n=pdb->maxresid0()+100;
	atmat=new Res*[n];
	int i;
	for(i=0;i<n;i++) atmat[i]=0;

	Chn *c,*c0;
	for(c=pdb->chn;c;c=c->next){
		c0=pdborg->ischain(c->id);
		if(c0==0) continue;
		Res *r,*r0;
		for(r=c->res;r;r=r->next) {
			r0=c0->isresoid(r->oid);
			if(r0==0) continue;
			atmat[r->id0]=r0;
		}
	}	

}

void PdbFix::mynewfix(Chn *chn) {

	Rotate rot;

	chn->setlastresorder();
	
	Res *r0=chn->lastres();  
	Res *r;
	
	Chn *rotamer=TRES.findrotamer("backbone")->chn;
	if(rotamer==0) {
		cerr<<"backbone rotamer not found..."<<endl;
		exit(0);
	}

	Res *t=0;	
	for(r=r0;r;r=r->last) {			
		float end=0;
		if(t==0||t->next==0||t->next->id-t->id!=1) end=1; 
		Res *r1=rotamer->isres(r->name,0);
		Res *r2;
		
		float dmin=100000;
		
		int cis=0;
		if(t&&atmat[t->id0]&&atmat[r->id0]) {
			Atm *ca0=atmat[t->id0]->isatmid(1);
			Atm *ca1=atmat[r->id0]->isatmid(1);
			if(ca0&&ca1) {
				float e=TRES.distance(ca0,ca1);
				if(e<3) cis=1;
			}
		}	
		 
		Res *rs=0;
		float txyz[200];
		float thead[200];

		//save t conformation
		if(t) {
			t->transfer(txyz,0);
			t->transfertemp(thead,0);
		}
		 
		for(r2=r1->more;r2;r2=r2->more) {
			if(r->sec=='h'&&t&&r->last&&t->sec=='h'&&r->last->sec=='h') {
				if(!(fabs(r2->atm->next->chi+60)<40&&fabs(r2->atm->next->next->chi+45)<30)) continue;				 
			}
			else if(r->sec=='e'&&t&&r->last&&t->sec=='e'&&r->last->sec=='e') {
				if(!(fabs(r2->atm->next->chi+110)<70&&fabs(r2->atm->next->next->chi-130)<90)) continue;	
			}
			r->transfer(r2);
			r->copytemp(r2);
			if(t) {
				t->transfer(txyz,1);
				t->transfertemp(thead,1);
				rot.link(r,t,r->id0,0);
				if(cis) {					
					float rt=t->atm->gettorsionangle();						
					rot.rotateomega(t->atm,rt,r->id0,0);						
				}
																		
			}
			float d;
			if(end) d=getendrmsd(r,t,1);  
			else    d=getrmsd(r,t,1);
			if(d<dmin) {
				dmin=d;
				rs=r2;	 	
			}
		}
		 
		 
		r->transfer(rs);
		r->copytemp(rs);	
		
		if(TRES.logg>3) cerr<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
		
		if(t) {
			t->transfer(txyz,1);
			t->transfertemp(thead,1);
			rot.link(r,t,r->id0,0);	
		}
		 		
		//rotate continously to find the local minima
 		if(!t) {
			r->transfer(txyz,0);
			r->transfertemp(thead,0);
		}
		else {
			chn->transfer(txyz,r,t->id0,0);
			chn->transfertemp(thead,r,t->id0,0);
		}
		if(end) dmin=getendrmsd(r,t,1);
		else	dmin=getrmsd(r,t,1);

		Atm *g1,*g2,*g3;
		 		 
		re200:
		float dmin0=dmin;	
		int aa=5;
		while(t&&t->id-r->id==1) {
			if(dmin<0.01) break;
			int a,b,c;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			for(c=-1;c<2;c++) {
				float e,f,g;
				e=a*aa;f=b*aa;g=c*aa;				
				chn->transfer(txyz,r,t->id0,1);
				chn->transfertemp(thead,r,t->id0,1);
				g1=t->isatmid(2);
				g2=t->isatmid(1);
				g3=t->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;
								
				rot.rotate(g1,e,r->id0,0);
				rot.rotate(g2,f,r->id0,0);
				rot.rotateomega(g3,g,r->id0,0);						
				
				float d;
				if(end) d=getendrmsd(r,t,1); //1, all fixed atm, 0 this atom	
				else    d=getrmsd(r,t,1);
				if(d<dmin) {
					float rt=g3->gettorsionangle();
					if(fabs(rt)>160&&cis==0) {	
						if(fabs(dmin-d)>0.01) nn++;
						dmin=d;
						chn->transfer(txyz,r,t->id0,0);
						chn->transfertemp(thead,r,t->id0,0);
					}
					else if(fabs(rt)<20&&cis==1) {	
						if(fabs(dmin-d)>0.01) nn++;
						dmin=d;
						chn->transfer(txyz,r,t->id0,0);
						chn->transfertemp(thead,r,t->id0,0);
					}
				}
			}
 
			if(nn==0) aa=aa/2;
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);
			if(aa==0) break;
		}

		if(TRES.logg>3) cerr<<"rotate t:"<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;		
		aa=5;
		while(t&&t->id-r->id==1) {
			if(dmin<0.01) break;
			int a,b,c;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			for(c=-1;c<0;c++) {
				float e,f;
				e=a*aa;f=b*aa;
				chn->transfer(txyz,r,t->id0,1);
				chn->transfertemp(thead,r,t->id0,1);
				g1=r->isatmid(2);
				g2=r->isatmid(1);
				g3=r->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;				
				rot.rotate(g1,e,r->id0,0);
				rot.rotate(g2,f,r->id0,0);
				
				float d;
				if(end) d=getendrmsd(r,t,1); //1, all fixed atm, 0 this atom	
				else    d=getrmsd(r,t,1);
				if(d<dmin) {
					dmin=d;
					chn->transfer(txyz,r,t->id0,0);
					chn->transfertemp(thead,r,t->id0,0);
					nn++;
				}
			}
			if(nn==0) aa=aa/2;
						 
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);
			if(aa==0) break;
		}
		if(!t) {
			r->transfer(txyz,1);
			r->transfertemp(thead,1);
		}
		else {
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);	
		}
		if(TRES.logg>3) cerr<<"rotate r: "<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
	 
		//rotate the bond so that it is more closer to the exact
		//start exacting the bond to the original
		
		//end
		if(fabs(dmin-dmin0)>0.05) {
			dmin0=dmin;
			goto re200;
		}
		Res *r0=atmat[r->id0];
		Atm *a0,*a;
		if(exact&&r0&&r0->next&&r0->next->id-r0->id==1&&t&&t->id-r->id==1) {
			 
			a0=r0->isatmid(1);		
			
			a=r->isatmid(1);
			
			 
			if(a0&&a) {
				float d=TRES.distance(a0,a);
				if(d<0.1) {
					float dx=(a0->xyz[0]-a->xyz[0]);
					float dy=(a0->xyz[1]-a->xyz[1]);
					float dz=(a0->xyz[2]-a->xyz[2]);
					Atm *aa;
					for(aa=r->atm;aa;aa=aa->next) {
						aa->xyz[0]+=dx;
                                		aa->xyz[1]+=dy;
                                		aa->xyz[2]+=dz;
					}
					int i;
					for(i=0;i<3;i++) {
						r->temp[i*3+0]+=dx;
						r->temp[i*3+1]+=dy;
						r->temp[i*3+2]+=dz;
					}
					//cerr<<d<<" "<<d0<<" "<<TRES.distance(a,b)<<endl;
				}
			}
		}
		
		if(r->last) rot.link(r->last,r,r->last->id0,0);
		if(t&&end==1&&r->last)  superimposesegment(r->last,t,1);
		else if(t&&end==1) superimposesegment(r,t,1);
		if(TRES.logg>3) pdb->write("ss");
		
		t=r;			
	}	
}

void PdbFix::mynewfix(Res *first, Res *last) {

	Rotate rot;
	Chn *chn=first->chn;
	chn->setlastresorder();
	
	Res *rr0=last;  
	Res *r;
	
	Chn *rotamer=TRES.findrotamer("backbone")->chn;
	if(rotamer==0) {
		cerr<<"backbone rotamer not found..."<<endl;
		exit(0);
	}

	Res *t=chn->isres(rr0->id+1);	
	for(r=rr0;r;r=r->last) {
		if(r->id0<first->id0) break;			
		float end=0;
		if(t==0||t->next==0||t->next->id-t->id!=1) end=1; 
		Res *r1=rotamer->isres(r->name,0);
		Res *r2;
		
		float dmin=100000;
		
		int cis=0;
		if(t&&atmat[t->id0]&&atmat[r->id0]) {
			Atm *ca0=atmat[t->id0]->isatmid(1);
			Atm *ca1=atmat[r->id0]->isatmid(1);
			if(ca0&&ca1&&t&&t->name=='P') {
                                float e=TRES.distance(ca0,ca1);
                                if(e<3.3) cis=1;
                        }
		}	
		 
		Res *rs=0;
		float txyz[200];
		float thead[200];

		//save t conformation
		if(t) {
			t->transfer(txyz,0);
			t->transfertemp(thead,0);
		}
		 
		for(r2=r1->more;r2;r2=r2->more) {
			if(r->sec=='h'&&t&&r->last&&t->sec=='h'&&r->last->sec=='h') {
				if(!(fabs(r2->atm->next->chi+60)<40&&fabs(r2->atm->next->next->chi+45)<30)) continue;				 
			}
			else if(r->sec=='e'&&t&&r->last&&t->sec=='e'&&r->last->sec=='e') {
				if(!(fabs(r2->atm->next->chi+110)<70&&fabs(r2->atm->next->next->chi-130)<90)) continue;	
			}
			r->transfer(r2);
			r->copytemp(r2);
			if(t) {
				t->transfer(txyz,1);
				t->transfertemp(thead,1);
				rot.link(r,t,r->id0,0);
				if(cis) {					
					float rt=t->atm->gettorsionangle();						
					rot.rotateomega(t->atm,rt,r->id0,0);						
				}
																		
			}
			float d;
			if(end) d=getendrmsd(r,t,1);  
			else    d=getrmsd(r,t,1);
			if(d<dmin) {
				dmin=d;
				rs=r2;	 	
			}
		}
		 
		 
		r->transfer(rs);
		r->copytemp(rs);	
		
		if(TRES.logg>3) cerr<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
		
		if(t) {
			t->transfer(txyz,1);
			t->transfertemp(thead,1);
			rot.link(r,t,r->id0,0);	
		}
		 		
		//rotate continously to find the local minima
 		if(!t) {
			r->transfer(txyz,0);
			r->transfertemp(thead,0);
		}
		else {
			chn->transfer(txyz,r,t->id0,0);
			chn->transfertemp(thead,r,t->id0,0);
		}
		if(end) dmin=getendrmsd(r,t,1);
		else	dmin=getrmsd(r,t,1);

		Atm *g1,*g2,*g3;
		 		 
		re200:
		float dmin0=dmin;	
		int aa=5;
		while(t&&t->id-r->id==1) {
			if(dmin<0.01) break;
			int a,b,c;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			for(c=-1;c<2;c++) {
				float e,f,g;
				e=a*aa;f=b*aa;g=c*aa;				
				chn->transfer(txyz,r,t->id0,1);
				chn->transfertemp(thead,r,t->id0,1);
				g1=t->isatmid(2);
				g2=t->isatmid(1);
				g3=t->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;
								
				rot.rotate(g1,e,r->id0,0);
				rot.rotate(g2,f,r->id0,0);
				rot.rotateomega(g3,g,r->id0,0);						
				
				float d;
				if(end) d=getendrmsd(r,t,1); //1, all fixed atm, 0 this atom	
				else    d=getrmsd(r,t,1);
				if(d<dmin) {
					float rt=g3->gettorsionangle();
					if(fabs(rt)>160&&cis==0) {	
						if(fabs(dmin-d)>0.01) nn++;
						dmin=d;
						chn->transfer(txyz,r,t->id0,0);
						chn->transfertemp(thead,r,t->id0,0);
					}
					else if(fabs(rt)<20&&cis==1) {	
						if(fabs(dmin-d)>0.01) nn++;
						dmin=d;
						chn->transfer(txyz,r,t->id0,0);
						chn->transfertemp(thead,r,t->id0,0);
					}
				}
			}
 
			if(nn==0) aa=aa/2;
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);
			if(aa==0) break;
		}

		if(TRES.logg>3) cerr<<"rotate t:"<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;		
		aa=5;
		while(t&&t->id-r->id==1) {
			if(dmin<0.01) break;
			int a,b,c;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			for(c=-1;c<0;c++) {
				float e,f;
				e=a*aa;f=b*aa;
				chn->transfer(txyz,r,t->id0,1);
				chn->transfertemp(thead,r,t->id0,1);
				g1=r->isatmid(2);
				g2=r->isatmid(1);
				g3=r->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;				
				rot.rotate(g1,e,r->id0,0);
				rot.rotate(g2,f,r->id0,0);
				
				float d;
				if(end) d=getendrmsd(r,t,1); //1, all fixed atm, 0 this atom	
				else    d=getrmsd(r,t,1);
				if(d<dmin) {
					dmin=d;
					chn->transfer(txyz,r,t->id0,0);
					chn->transfertemp(thead,r,t->id0,0);
					nn++;
				}
			}
			if(nn==0) aa=aa/2;
						 
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);
			if(aa==0) break;
		}
		if(!t) {
			r->transfer(txyz,1);
			r->transfertemp(thead,1);
		}
		else {
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);	
		}
		if(TRES.logg>3) cerr<<"rotate r: "<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
	 
		//rotate the bond so that it is more closer to the exact
		//start exacting the bond to the original
		
		//end
		if(fabs(dmin-dmin0)>0.05) {
			dmin0=dmin;
			goto re200;
		}
		Res *r0=atmat[r->id0];
		Atm *a0,*a;
		if(exact&&r0&&r0->next&&r0->next->id-r0->id==1&&t&&t->id-r->id==1) {
			 
			a0=r0->isatmid(1);		
			
			a=r->isatmid(1);
			
			 
			if(a0&&a) {
				float d=TRES.distance(a0,a);
				if(d<0.1) {
					float dx=(a0->xyz[0]-a->xyz[0]);
					float dy=(a0->xyz[1]-a->xyz[1]);
					float dz=(a0->xyz[2]-a->xyz[2]);
					Atm *aa;
					for(aa=r->atm;aa;aa=aa->next) {
						aa->xyz[0]+=dx;
                                		aa->xyz[1]+=dy;
                                		aa->xyz[2]+=dz;
					}
					int i;
					for(i=0;i<3;i++) {
						r->temp[i*3+0]+=dx;
						r->temp[i*3+1]+=dy;
						r->temp[i*3+2]+=dz;
					}
					//cerr<<d<<" "<<d0<<" "<<TRES.distance(a,b)<<endl;
				}
			}
		}
		
		
		if(r->last&&r->last->id0>=first->id0) rot.link(r->last,r,r->last->id0,0);		
		//if(r->last) rot.link(r->last,r,r->last->id0,0);
		if(rr0->next==0) {	
			if(t&&end==1&&r->last)  superimposesegment(r->last,t,1);
			else if(t&&end==1) superimposesegment(r,t,1);
		}	
		if(TRES.logg>3) pdb->write("ss");
		
		t=r;			
	}	
}



void PdbFix::mynewfix0(Res *first,Res *last) {

	Rotate rot;
	Chn *chn=first->chn;
	chn->setlastresorder();
	
	Res *rr0=last;//chn->lastres();  
	Res *r;
	
	Chn *rotamer=TRES.findrotamer("backbone")->chn;
	if(rotamer==0) {
		cerr<<"backbone rotamer not found..."<<endl;
		exit(0);
	}

	Res *t=first->chn->isres(last->id+1);
	if(t&&t->id-rr0->id!=1) t=0;	
	for(r=rr0;r;r=r->last) {	
		if(r->id0<first->id0) break;		
		float end=0;
		if(t==0||t->next==0||t->next->id-t->id!=1) end=1; 
		Res *r1=rotamer->isres(r->name,0);
		Res *r2;
		if(t&&t->id-r->id!=1) t=0;		

		float dmin=100000;
		
		int cis=0;
		if(t&&atmat[t->id0]&&atmat[r->id0]&&t->id-r->id==1) {
			Atm *ca0=atmat[t->id0]->isatmid(1);
			Atm *ca1=atmat[r->id0]->isatmid(1);
			if(ca0&&ca1&&t&&t->name=='P') {
				float e=TRES.distance(ca0,ca1);
				if(e<3.3) cis=1;
			}
		}	
		 
		Res *rs=0;
		float txyz[200];
		float thead[200];

		//save t conformation
		if(t) {
			t->transfer(txyz,0);
			t->transfertemp(thead,0);
		}
		 
		for(r2=r1->more;r2;r2=r2->more) {
			if(r->sec=='h'&&t&&r->last&&t->sec=='h'&&r->last->sec=='h') {
				if(!(fabs(r2->atm->next->chi+60)<40&&fabs(r2->atm->next->next->chi+45)<30)) continue;				 
			}
			else if(r->sec=='e'&&t&&r->last&&t->sec=='e'&&r->last->sec=='e') {
				if(!(fabs(r2->atm->next->chi+110)<70&&fabs(r2->atm->next->next->chi-130)<90)) continue;	
			}
			r->transfer(r2);
			r->copytemp(r2);
			if(t) {
				t->transfer(txyz,0);
				t->transfertemp(thead,0);
			}
			if(t&&t->id-r->id==1) {
				t->transfer(txyz,1);
				t->transfertemp(thead,1);
				rot.link(r,t,r->id0,0);
				if(cis) {					
					float rt=t->atm->gettorsionangle();						
					rot.rotateomega(t->atm,rt,r->id0,0);
					if(TRES.logg>3) cerr<<t->atm->gettorsionangle()<<endl;						
				}
																		
			}
			float d;
			if(t&&t->id-r->id==1) {
				if(end==1) d=getendrmsd(r,t,1);  
				else  if(end==0) d=getrmsd(r,t,1);
			}
			else {
				if(end==1) d=getendrmsd(r,0,1);  
				else  if(end==0) d=getrmsd(r,0,1);
			}
			if(d<dmin) {
				dmin=d;
				rs=r2;	 	
			}
		}
		 
		 
		r->transfer(rs);
		r->copytemp(rs);	
		
		if(TRES.logg>3) cerr<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
		if(t) {
			t->transfer(txyz,1);
			t->transfertemp(thead,1);
		}
		if(t&&t->id-r->id==1) {
			t->transfer(txyz,1);
			t->transfertemp(thead,1);
			rot.link(r,t,r->id0,0);	
			if(cis) {					
				float rt=t->atm->gettorsionangle();						
				rot.rotateomega(t->atm,rt,r->id0,0);
				if(TRES.logg>3) cerr<<t->atm->gettorsionangle()<<endl;						
			}
		}
		 		
		//rotate continously to find the local minima
		if(t&&t->id-r->id==1) {
			chn->transfer(txyz,r,t->id0,0);
			chn->transfertemp(thead,r,t->id0,0);
		}
 		else if(t) {
			r->transfer(txyz,0);
			r->transfertemp(thead,0);
		}
		
		if(t&&t->id-r->id==1) {
			if(end==1) dmin=getendrmsd(r,t,1);  
			else  if(end==0) dmin=getrmsd(r,t,1);
		}
		else {
			if(end==1) dmin=getendrmsd(r,0,1);  
			else  if(end==0) dmin=getrmsd(r,0,1);
		}

		Atm *g1,*g2,*g3;
		 		 
		re200:
		float dmin0=dmin;	
		int aa=5;
		while(t&&t->id-r->id==1) {
			if(dmin<0.01) break;
			int a,b,c;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			for(c=-1;c<2;c++) {
				float e,f,g;
				e=a*aa;f=b*aa;g=c*aa;				
				chn->transfer(txyz,r,t->id0,1);
				chn->transfertemp(thead,r,t->id0,1);
				g1=t->isatmid(2);
				g2=t->isatmid(1);
				g3=t->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;
								
				rot.rotate(g1,e,r->id0,0);
				rot.rotate(g2,f,r->id0,0);
				rot.rotateomega(g3,g,r->id0,0);						
				
				float d;
				if(end) d=getendrmsd(r,t,1); //1, all fixed atm, 0 this atom	
				else    d=getrmsd(r,t,1);
				if(d<dmin) {
					float rt=g3->gettorsionangle();
					if(fabs(rt)>170&&cis==0) {	
						if(fabs(dmin-d)>0.01) nn++;
						dmin=d;
						chn->transfer(txyz,r,t->id0,0);
						chn->transfertemp(thead,r,t->id0,0);
					}
					else if(fabs(rt)<10&&cis==1) {	
						if(fabs(dmin-d)>0.01) nn++;
						dmin=d;
						chn->transfer(txyz,r,t->id0,0);
						chn->transfertemp(thead,r,t->id0,0);
					}
				}
			}
 
			if(nn==0) aa=aa/2;
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);
			if(aa==0) break;
		}

		if(TRES.logg>3) cerr<<"rotate t:"<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;		
		aa=5;
		while(t&&t->id-r->id==1) {
			if(dmin<0.01) break;
			int a,b,c;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			for(c=-1;c<0;c++) {
				float e,f;
				e=a*aa;f=b*aa;
				chn->transfer(txyz,r,t->id0,1);
				chn->transfertemp(thead,r,t->id0,1);
				g1=r->isatmid(2);
				g2=r->isatmid(1);
				g3=r->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;				
				rot.rotate(g1,e,r->id0,0);
				rot.rotate(g2,f,r->id0,0);
				
				float d;
				if(end) d=getendrmsd(r,t,1); //1, all fixed atm, 0 this atom	
				else    d=getrmsd(r,t,1);
				if(d<dmin) {
					dmin=d;
					chn->transfer(txyz,r,t->id0,0);
					chn->transfertemp(thead,r,t->id0,0);
					nn++;
				}
			}
			if(nn==0) aa=aa/2;
						 
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);
			if(aa==0) break;
		}
		if(t&&t->id-r->id==1) {
			chn->transfer(txyz,r,t->id0,0);
			chn->transfertemp(thead,r,t->id0,0);
		}
 		else if(t) {
			r->transfer(txyz,0);
			r->transfertemp(thead,0);
		}
		if(TRES.logg>3) cerr<<"rotate r: "<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
	 
		//rotate the bond so that it is more closer to the exact
		//start exacting the bond to the original
		
		//end
		if(fabs(dmin-dmin0)>0.05) {
			dmin0=dmin;
			goto re200;
		}
		Res *r0=atmat[r->id0];
		Atm *a0,*a;
		if(exact&&r0&&r0->next&&r0->next->id-r0->id==1&&t&&t->id-r->id==1) {
			 
			a0=r0->isatmid(1);		
			
			a=r->isatmid(1);
			
			 
			if(a0&&a) {
				float d=TRES.distance(a0,a);
				if(d<0.1) {
					float dx=(a0->xyz[0]-a->xyz[0]);
					float dy=(a0->xyz[1]-a->xyz[1]);
					float dz=(a0->xyz[2]-a->xyz[2]);
					Atm *aa;
					for(aa=r->atm;aa;aa=aa->next) {
						aa->xyz[0]+=dx;
                                		aa->xyz[1]+=dy;
                                		aa->xyz[2]+=dz;
					}
					int i;
					for(i=0;i<3;i++) {
						r->temp[i*3+0]+=dx;
						r->temp[i*3+1]+=dy;
						r->temp[i*3+2]+=dz;
					}
					//cerr<<d<<" "<<d0<<" "<<TRES.distance(a,b)<<endl;
				}
			}
		}
		if(t&&t->id-r->id==1) {
			if(r->last&&r->id-r->last->id==1) {
				rot.link(r->last,r,r->last->id0,0);
				if(end==1) superimposesegment(r->last,t,1); 
			}
			else {
				if(end==1) superimposesegment(r,t,1); 
			}	
		}
		else {
			if(r->last&&r->id-r->last->id==1) {
				rot.link(r->last,r,r->last->id0,0);
				if(end==1) superimposesegment(r->last,r,1); 
			}
			else {
				if(end==1) superimposesegment(r,r,1); 
			}
		}
		if(TRES.logg>3) pdb->write("ss");
		
		t=r;			
	}	
}
void PdbFix::mynewfixmoremore(Chn *chn) {

	Rotate rot;

	chn->setlastresorder();
	
	Res *r0=chn->lastres();	
	Res *r;
	
	Chn *rotamer=TRES.findrotamer("backbone")->chn;
	if(rotamer==0) {
		cerr<<"backbone rotamer not found..."<<endl;
		exit(0);
	}
	Res *t=0;
	for(r=r0;r;r=r->last) {

		float end=0;
		if(t==0||t->next==0||t->next->id-t->id!=1) {
			t=r;
			continue;
		}
		if(r->id-t->id!=-1) {
			t=r;
			continue;
		}

		Res *r1=rotamer->isres(r->name,0);
		Res *r2;
		
		float dmin=100000;
		
		int cis=0;
		if(t&&atmat[t->id0]&&atmat[r->id0]) {
			Atm *ca0=atmat[t->id0]->isatmid(1);
			Atm *ca1=atmat[r->id0]->isatmid(1);
			if(ca0&&ca1) {
				float e=TRES.distance(ca0,ca1);
				if(e<3) cis=1;
			}
		}	
		 
 
		float txyz[200];
		float thead[200];

		//save t conformation
		if(t) {
			chn->transfer(txyz,r,t->id0,0);
			chn->transfertemp(thead,r,t->id0,0);
		}
		dmin=getrmsd(r,t,1);
		for(r2=r1->more;r2;r2=r2->more) {
			if(r->sec=='h'&&t&&r->last&&t->sec=='h'&&r->last->sec=='h') {
				if(!(fabs(r2->atm->next->chi+60)<40&&fabs(r2->atm->next->next->chi+45)<30)) continue;				 
			}
			else if(r->sec=='e'&&t&&r->last&&t->sec=='e'&&r->last->sec=='e') {
				if(!(fabs(r2->atm->next->chi+110)<70&&fabs(r2->atm->next->next->chi-130)<90)) continue;	
			}
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);
			r->transfer(r2);
			r->copytemp(r2);
			if(t) {
				//t->transfer(txyz,1);
				//t->transfertemp(thead,1);
				rot.link(r,t,r->id0,0);
				if(cis) {					
					float rt=t->atm->gettorsionangle();						
					rot.rotateomega(t->atm,rt,r->id0,0);						
				}
																		
			}
			float d;
			if(end) d=getendrmsd(r,t,1);  
			else    d=getrmsd(r,t,1);
			if(d<dmin) {
				dmin=d;
				chn->transfer(txyz,r,t->id0,0);
				chn->transfertemp(thead,r,t->id0,0);
				//rs=r2;	 	
			}
		}
		chn->transfer(txyz,r,t->id0,1);
		chn->transfertemp(thead,r,t->id0,1);
		//r->transfer(rs);
		//r->copytemp(rs);	
		
		if(TRES.logg>3) cerr<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
		
		if(t) {
			//t->transfer(txyz,1);
			//t->transfertemp(thead,1);
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);
			rot.link(r,t,r->id0,0);	
		}
		 		
		//rotate continously to find the local minima
 		if(!t) {
			r->transfer(txyz,0);
			r->transfertemp(thead,0);
		}
		else {
			chn->transfer(txyz,r,t->id0,0);
			chn->transfertemp(thead,r,t->id0,0);
		}
		if(end) dmin=getendrmsd(r,t,1);
		else	dmin=getrmsd(r,t,1);

		Atm *g1,*g2,*g3;
		 		 
		re200:
		float dmin0=dmin;	
		int aa=5;
		while(1&&t&&t->id-r->id==1) {
			if(dmin<0.01) break;
			int a,b,c;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			for(c=-1;c<2;c++) {
				float e,f,g;
				e=a*aa;f=b*aa;g=c*aa;				
				chn->transfer(txyz,r,t->id0,1);
				chn->transfertemp(thead,r,t->id0,1);
				g1=t->isatmid(2);
				g2=t->isatmid(1);
				g3=t->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;
								
				rot.rotate(g1,e,r->id0,0);
				rot.rotate(g2,f,r->id0,0);
				rot.rotateomega(t->atm,g,r->id0,0);						
				//float rt=t->atm->gettorsionangle();						
				//rot.rotateomega(t->atm,rt+180+g,r->id0,0);
				//if(fabs(rt-g)>170) rot.rotateomega(g3,g,r->id0,0);	
				//cerr<<rt-g<<" "<<t->atm->gettorsionangle()<<" "<<dmin<<endl;		 
				float d;
				if(end) d=getendrmsd(r,t,1); //1, all fixed atm, 0 this atom	
				else    d=getrmsd(r,t,1);
				if(d<dmin) {
					float rt=t->atm->gettorsionangle();
					if(fabs(rt)>160&&cis==0) {	
						if(fabs(dmin-d)>0.01) nn++;
						dmin=d;
						chn->transfer(txyz,r,t->id0,0);
						chn->transfertemp(thead,r,t->id0,0);
					}
					else if(fabs(rt)<20&&cis==1) {	
						if(fabs(dmin-d)>0.01) nn++;
						dmin=d;
						chn->transfer(txyz,r,t->id0,0);
						chn->transfertemp(thead,r,t->id0,0);
					}
				}
			}
 
			if(nn==0) aa=aa/2;
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);
			if(aa==0) break;
		}

		if(TRES.logg>3) cerr<<"rotate t:"<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;		
		aa=5;
		while(1&&t&&t->id-r->id==1) {
			if(dmin<0.01) break;
			int a,b,c;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			for(c=-1;c<0;c++) {
				float e,f;
				e=a*aa;f=b*aa;
				chn->transfer(txyz,r,t->id0,1);
				chn->transfertemp(thead,r,t->id0,1);
				g1=r->isatmid(2);
				g2=r->isatmid(1);
				g3=r->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;				
				rot.rotate(g1,e,r->id0,0);
				rot.rotate(g2,f,r->id0,0);
				//float rt=r->atm->gettorsionangle();						
				//rot.rotateomega(t->atm,rt+180+g,r->id0,0);
				//if(fabs(rt-g)>170) rot.rotateomega(g3,g,r->id0,0);
				float d;
				if(end) d=getendrmsd(r,t,1); //1, all fixed atm, 0 this atom	
				else    d=getrmsd(r,t,1);
				if(d<dmin) {
					dmin=d;
					chn->transfer(txyz,r,t->id0,0);
					chn->transfertemp(thead,r,t->id0,0);
					nn++;
				}
			}
			if(nn==0) aa=aa/2;
						 
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);
			if(aa==0) break;
		}
		if(!t) {
			r->transfer(txyz,1);
			r->transfertemp(thead,1);
		}
		else {
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);	
		}
		if(TRES.logg>3) cerr<<"rotate r: "<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
		
		//start
		
		if(1||end||r->last==0||r->id-r->last->id!=1) goto re400;

		float ttxyz[200],tthead[200];
			
		t->transfer(ttxyz,0);
		t->transfertemp(tthead,0);	
		rot.link(r->last,r,r->last->id0,0); 
		float ang1;ang1=r->isatmid(1)->gettorsionangle();
		float ang2;ang2=r->isatmid(2)->gettorsionangle();
		float ang3;ang3=t->atm->gettorsionangle();
		
		for(r2=r1->more;r2;r2=r2->more) {
			r->transfer(r2);
			r->copytemp(r2);
				
			t->transfer(ttxyz,1);
			t->transfertemp(tthead,1);
			rot.link(r,t,r->id0,0);			
			if(cis) {					
				float rt=t->atm->gettorsionangle();						
				rot.rotateomega(t->atm,rt,r->id0,0);						
			}	
			float rt;													
			rot.link(r->last,r,r->last->id0,0); 
			rt=t->atm->gettorsionangle();
			rot.rotateomega(t->atm,rt-ang3,r->last->id0,0);
			//a
			Atm *a=r->isatmid(2);
			rt=a->gettorsionangle();
			rot.rotate(a,rt-ang2,r->last->id0,0);
			//a
			a=r->isatmid(1);
			rt=a->gettorsionangle();
			rot.rotate(a,rt-ang1,r->last->id0,0);			 
			// 
			float d;
			if(end) d=getendrmsd(r,t,1);  
			else    d=getrmsd(r,t,1);
			if(d<dmin) {
					dmin=d;					
					chn->transfer(txyz,r,t->id0,0);
					chn->transfertemp(thead,r,t->id0,0);							
			}
		}
		//end
		 
		re400:
		
		if(fabs(dmin-dmin0)>0.05) {
			dmin0=dmin;
			goto re200;
		}
		//end
		if(!t) {
			r->transfer(txyz,1);
			r->transfertemp(thead,1);
		}
		else {
			chn->transfer(txyz,r,t->id0,1);
			chn->transfertemp(thead,r,t->id0,1);	
		}
		//rotate the bond so that it is more closer to the exact
		//start exacting the bond to the original
		Res *r0=atmat[r->id0];
		Atm *a0,*b0,*a,*b;
		if(exact&&r0&&r0->next&&r0->next->id-r0->id==1&&t&&t->id-r->id==1) {
			 
			a0=r0->isatmid(1);		
			b0=r0->next->isatmid(1);
			a=r->isatmid(1);
			b=r->next->isatmid(1);
			/*
			if(a0&&b0&&a&&b) {
				float d=TRES.distance(a0,b0);
				if(d>3.9) d=3.9;
				else if(d<3.6&&d>3.0) d=3.6;
				float d0=TRES.distance(a,b);
				float x,y,z;
                        	float dx,dy,dz;
				dx=(a->xyz[0]-b->xyz[0])/d0;
                        	dy=(a->xyz[1]-b->xyz[1])/d0;
                        	dz=(a->xyz[2]-b->xyz[2])/d0;
                        	x=(d-d0)*dx;
                        	y=(d-d0)*dy;
                        	z=(d-d0)*dz;	
				Atm *aa;
				for(aa=r->atm;aa;aa=aa->next) {
					aa->xyz[0]+=x;
                                	aa->xyz[1]+=y;
                                	aa->xyz[2]+=z;
				}
				int i;
				for(i=0;i<3;i++) {
					r->temp[i*3+0]+=x;
					r->temp[i*3+1]+=y;
					r->temp[i*3+2]+=z;
				}
				cerr<<d<<" "<<d0<<" "<<TRES.distance(a,b)<<endl;
			}
			*/ 
			if(a0&&b0&&a&&b) {
				float d=TRES.distance(a0,a);
				if(d<0.1) {
					float dx=(a0->xyz[0]-a->xyz[0]);
					float dy=(a0->xyz[1]-a->xyz[1]);
					float dz=(a0->xyz[2]-a->xyz[2]);
					Atm *aa;
					for(aa=r->atm;aa;aa=aa->next) {
						aa->xyz[0]+=dx;
                                		aa->xyz[1]+=dy;
                                		aa->xyz[2]+=dz;
					}
					int i;
					for(i=0;i<3;i++) {
						r->temp[i*3+0]+=dx;
						r->temp[i*3+1]+=dy;
						r->temp[i*3+2]+=dz;
					}
					//cerr<<d<<" "<<d0<<" "<<TRES.distance(a,b)<<endl;
				}
			}
		}
		//end


		if(r->last) rot.link(r->last,r,r->last->id0,0);
		//if(t&&end==1&&r->last)  superimposesegment(r->last,t,1);
		//else if(t&&end==1) superimposesegment(r,t,1);
		if(TRES.logg>3) pdb->write("ss");
		
		
		t=r;			
	}	
}

void PdbFix::mynewfixmore(Chn *chn) {

	Rotate rot;

	chn->setlastresorder();
	
		
	Res *r;
	
	Chn *rotamer=TRES.findrotamer("backbone")->chn;
	if(rotamer==0) {
		cerr<<"backbone rotamer not found..."<<endl;
		exit(0);
	}
	Res *t=0;
	
	for(r=chn->res;r;r=r->next) {

		int isnext=0; //if next exist.0, no;1, yes.

		//if(t==0||t->last==0||t->id-t->last->id!=1) {
		if(t==0){
			t=r;
			continue;
		}
		else if(r->id-t->id!=1) {
			t=r;
			continue;
		}
		 
		if(r->next&&r->next->id-r->id==1) isnext=1;

		Res *r1=rotamer->isres(r->name,0);
		Res *r2;
		
		float dmin=100000;
		
		int cis=0;
		if(t&&atmat[t->id0]&&atmat[r->id0]) {
			Atm *ca0=atmat[t->id0]->isatmid(1);
			Atm *ca1=atmat[r->id0]->isatmid(1);
			if(ca0&&ca1) {
				float e=TRES.distance(ca0,ca1);
				if(e<3) cis=1;
			}
		}	
		 
	
		float txyz[200];
		float thead[200];
		
		//save t conformation
		if(t) {
			rot.link(t,r,r->id0,1);
			if(isnext) rot.link(r,r->next,r->next->id0,1);	
			dmin=getrmsdmore(t,r,1);
			if(TRES.logg>3) cerr<<"original: "<<dmin<<endl;
			chn->transfer(txyz,t,r->id0+isnext,0);
			chn->transfertemp(thead,t,r->id0+isnext,0);
			float rt=r->atm->gettorsionangle();
			if(cis==1&&fabs(rt)>20) dmin=10000;
			else if(cis==0&&fabs(rt)<160) dmin=10000;
		}
		
		if(dmin<0.05) {t=r;continue;}

		for(r2=r1->more;r2;r2=r2->more) {
			if(r->sec=='h'&&t&&r->next&&t->sec=='h'&&r->next->sec=='h') {
				if(!(fabs(r2->atm->next->chi+60)<40&&fabs(r2->atm->next->next->chi+45)<30)) continue;				 
			}
			else if(r->sec=='e'&&t&&r->next&&t->sec=='e'&&r->next->sec=='e') {
				if(!(fabs(r2->atm->next->chi+110)<70&&fabs(r2->atm->next->next->chi-130)<90)) continue;	
			}
			chn->transfer(txyz,t,r->id0+isnext,1);
			chn->transfertemp(thead,t,r->id0+isnext,1);	
			r->transfer(r2);
			r->copytemp(r2);
			if(t) {
				//t->transfer(txyz,1);
				//t->transfertemp(thead,1);
				rot.link(t,r,r->id0,1);
				if(isnext) rot.link(r,r->next,r->next->id0,1);	
					
				if(cis) {					
					float rt=r->atm->gettorsionangle();						
					rot.rotateomega(r->atm,-rt,r->id0,1);						
				}														
			}
			float d;
			 
			d=getrmsdmore(t,r,1);
			if(d<dmin) {
				dmin=d;
				chn->transfer(txyz,t,r->id0+isnext,0);
				chn->transfertemp(thead,t,r->id0+isnext,0);				
			}
			if(TRES.logg>3) pdb->write("s0");
		}
		

		//r->transfer(rs);
		//r->copytemp(rs);	
		
		if(TRES.logg>3) cerr<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
		
		if(t) {
			chn->transfer(txyz,t,r->id0+isnext,1);
			chn->transfertemp(thead,t,r->id0+isnext,1);	
			//r->transfer(rs);
			//r->copytemp(rs);		
			//t->transfer(txyz,1);
			//t->transfertemp(thead,1);
			//rot.link(t,r,r->id0,1);	
			//if(isnext) rot.link(r,r->next,r->next->id0,1);
		}
		 		
		//rotate continously to find the local minima
 		
		chn->transfer(txyz,t,r->id0+isnext,0);
		chn->transfertemp(thead,t,r->id0+isnext,0);
		 
		 
		dmin=getrmsdmore(t,r,1);

		Atm *g1,*g2,*g3;
		 		 
		re200:
		float dmin0=dmin;	
		int aa=5;
		while(1&&t&&r->id-t->id==1) {
			if(dmin<0.01) break;
			int a,b,c;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			for(c=-1;c<2;c++) {
				float e,f,g;
				e=a*aa;f=b*aa;g=c*aa;				
				chn->transfer(txyz,t,r->id0+isnext,1);
				chn->transfertemp(thead,t,r->id0+isnext,1);
				g1=t->isatmid(1);
				g2=t->isatmid(2);
				g3=r->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;				
				rot.rotate(g1,e,r->id0+isnext,1);
				rot.rotate(g2,f,r->id0+isnext,1);	
				rot.rotateomega(g3,g,r->id0+isnext,1);				 
				float d;
				d=getrmsdmore(t,r,1);
				if(d<dmin) {
					float rt=r->atm->gettorsionangle();
                                        if(fabs(rt)>160&&cis==0) {
                                                if(fabs(dmin-d)>0.01) nn++;
                                                dmin=d;
                                                chn->transfer(txyz,t,r->id0+isnext,0);
						chn->transfertemp(thead,t,r->id0+isnext,0);
                                        }
                                        else if(fabs(rt)<20&&cis==1) {
                                                if(fabs(dmin-d)>0.01) nn++;
                                                dmin=d;
                                                chn->transfer(txyz,t,r->id0+isnext,0);
						chn->transfertemp(thead,t,r->id0+isnext,0);
                                        }					
				}
			}
 
			if(nn==0) aa=aa/2;
			chn->transfer(txyz,t,r->id0+isnext,1);
			chn->transfertemp(thead,t,r->id0+isnext,1);
			if(aa==0) break;
		}

		if(TRES.logg>3) cerr<<"rotate t:"<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;		
		aa=5;
		while(1&&t&&r->id-t->id==1) {
			if(dmin<0.01) break;
			int a,b;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			/*for(c=-1;c<0;c++)*/ {
				float e,f;
				e=a*aa;f=b*aa;
				chn->transfer(txyz,t,r->id0+isnext,1);
				chn->transfertemp(thead,t,r->id0+isnext,1);
				g1=r->isatmid(1);
				g2=r->isatmid(2);
				g3=r->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;				
				rot.rotate(g1,e,r->id0+isnext,1);
				rot.rotate(g2,f,r->id0+isnext,1);
				float d;
				d=getrmsdmore(t,r,1);
				if(d<dmin) {
					if(fabs(dmin-d)>0.01) nn++;
					dmin=d;
					chn->transfer(txyz,t,r->id0+isnext,0);
					chn->transfertemp(thead,t,r->id0+isnext,0);		
				}
			}
			if(nn==0) aa=aa/2;
						 
			chn->transfer(txyz,t,r->id0+isnext,1);
			chn->transfertemp(thead,t,r->id0+isnext,1);
			if(aa==0) break;
		}
		chn->transfer(txyz,t,r->id0+isnext,1);
		chn->transfertemp(thead,t,r->id0+isnext,1);	
	 
		if(TRES.logg>3) cerr<<"rotate r: "<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
 
		if(fabs(dmin-dmin0)>0.05) {
			dmin0=dmin;
			goto re200;
		}

		chn->transfer(txyz,t,r->id0+isnext,1);
		chn->transfertemp(thead,t,r->id0+isnext,1);	
			
		//rotate the bond so that it is more closer to the exact
		//start exacting the bond to the original
		Res *r0=atmat[r->id0];
		Atm *a0,*a;
		if(exact&&r0&&r0->next&&r0->next->id-r0->id==1&&t&&r->id-t->id==1) {
			 
			a0=r0->isatmid(1);					 
			a=r->isatmid(1);
			 
			
			if(a0&&a) {
				float d=TRES.distance(a0,a);
				if(d<0.1) {
					float dx=(a0->xyz[0]-a->xyz[0]);
					float dy=(a0->xyz[1]-a->xyz[1]);
					float dz=(a0->xyz[2]-a->xyz[2]);
					Atm *aa;
					for(aa=r->atm;aa;aa=aa->next) {
						aa->xyz[0]+=dx;
                                		aa->xyz[1]+=dy;
                                		aa->xyz[2]+=dz;
					}
					if(r->next) {	
						for(aa=r->next->atm;aa;aa=aa->next) {
        						aa->xyz[0]+=dx;
        						aa->xyz[1]+=dy;
        						aa->xyz[2]+=dz;
						}	
						int i;
						for(i=0;i<3;i++) {
							r->next->temp[i*3+0]+=dx;
							r->next->temp[i*3+1]+=dy;
							r->next->temp[i*3+2]+=dz;
						}
					}							 
				}
			}
		}
		//end
 
		if(TRES.logg>3) pdb->write("ss");
				
		t=r;			
	}	
}

void PdbFix::mydirectnewfix(Res *first, Res *last) {

	Rotate rot;
	Chn *chn=first->chn;
	chn->setlastresorder();
	
	Res *r,*r0;

	r0=first;

	Chn *rotamer=TRES.findrotamer("backbone")->chn;
	if(rotamer==0) {
		cerr<<"backbone rotamer not found..."<<endl;
		exit(0);
	}

	Res *t=chn->isres(first->id-1);
	for(r=r0;r;r=r->next) {
		if(r->id0>last->id0) break;
		int isnext=0; //if next exist.0, no;1, yes.

		//if(t==0||t->last==0||t->id-t->last->id!=1) {
		if(t==0){
			t=r;
			continue;
		}
		else if(r->id-t->id!=1) {
			t=r;
			continue;
		}
		 
		if(r->next&&r->next->id-r->id==1) isnext=1;

		Res *r1=rotamer->isres(r->name,0);
		Res *r2;
		
		float dmin=100000;
		
		int cis=0;
		if(t&&atmat[t->id0]&&atmat[r->id0]) {
			Atm *ca0=atmat[t->id0]->isatmid(1);
			Atm *ca1=atmat[r->id0]->isatmid(1);
			if(ca0&&ca1) {
				float e=TRES.distance(ca0,ca1);
				if(e<3) cis=1;
			}
		}	
		 	
		float txyz[200];
		float thead[200];
		
		//save t conformation
		if(t) {
			rot.link(t,r,r->id0,1);
			if(isnext) rot.link(r,r->next,r->next->id0,1);	
			dmin=getrmsdmore(t,r,1);
			if(TRES.logg>3) cerr<<"original: "<<dmin<<endl;
			chn->transfer(txyz,t,r->id0+isnext,0);
			chn->transfertemp(thead,t,r->id0+isnext,0);
			float rt=r->atm->gettorsionangle();
			if(cis==1&&fabs(rt)>20) dmin=10000;
			else if(cis==0&&fabs(rt)<160) dmin=10000;
		}
		
		if(dmin<0.05) {t=r;continue;}

		for(r2=r1->more;r2;r2=r2->more) {
			if(r->sec=='h'&&t&&r->next&&t->sec=='h'&&r->next->sec=='h') {
				if(!(fabs(r2->atm->next->chi+60)<40&&fabs(r2->atm->next->next->chi+45)<30)) continue;				 
			}
			else if(r->sec=='e'&&t&&r->next&&t->sec=='e'&&r->next->sec=='e') {
				if(!(fabs(r2->atm->next->chi+110)<70&&fabs(r2->atm->next->next->chi-130)<90)) continue;	
			}
			chn->transfer(txyz,t,r->id0+isnext,1);
			chn->transfertemp(thead,t,r->id0+isnext,1);	
			r->transfer(r2);
			r->copytemp(r2);
			if(t) {				 
				rot.link(t,r,r->id0,1);
				if(isnext) rot.link(r,r->next,r->next->id0,1);	
					
				if(cis) {					
					float rt=r->atm->gettorsionangle();						
					rot.rotateomega(r->atm,-rt,r->id0,1);						
				}														
			}
			float d;
			 
			d=getrmsdmore(t,r,1);
			if(d<dmin) {
				dmin=d;
				chn->transfer(txyz,t,r->id0+isnext,0);
				chn->transfertemp(thead,t,r->id0+isnext,0);				
			}
			if(TRES.logg>3) pdb->write("s0");
		}
 
		if(TRES.logg>3) cerr<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
		
		if(t) {
			chn->transfer(txyz,t,r->id0+isnext,1);
			chn->transfertemp(thead,t,r->id0+isnext,1);				 
		}
		 		
		//rotate continously to find the local minima
 		
		chn->transfer(txyz,t,r->id0+isnext,0);
		chn->transfertemp(thead,t,r->id0+isnext,0);
		 
		 
		dmin=getrmsdmore(t,r,1);

		Atm *g1,*g2,*g3;
		 		 
		re200:
		float dmin0=dmin;	
		int aa=5;
		while(1&&t&&r->id-t->id==1) {
			if(dmin<0.01) break;
			int a,b,c;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			for(c=-1;c<2;c++) {
				float e,f,g;
				e=a*aa;f=b*aa;g=c*aa;				
				chn->transfer(txyz,t,r->id0+isnext,1);
				chn->transfertemp(thead,t,r->id0+isnext,1);
				g1=t->isatmid(1);
				g2=t->isatmid(2);
				g3=r->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;				
				rot.rotate(g1,e,r->id0+isnext,1);
				rot.rotate(g2,f,r->id0+isnext,1);	
				rot.rotateomega(g3,g,r->id0+isnext,1);				 
				float d;
				d=getrmsdmore(t,r,1);
				if(d<dmin) {
					float rt=r->atm->gettorsionangle();
                                        if(fabs(rt)>160&&cis==0) {
                                                if(fabs(dmin-d)>0.01) nn++;
                                                dmin=d;
                                                chn->transfer(txyz,t,r->id0+isnext,0);
						chn->transfertemp(thead,t,r->id0+isnext,0);
                                        }
                                        else if(fabs(rt)<20&&cis==1) {
                                                if(fabs(dmin-d)>0.01) nn++;
                                                dmin=d;
                                                chn->transfer(txyz,t,r->id0+isnext,0);
						chn->transfertemp(thead,t,r->id0+isnext,0);
                                        }					
				}
			}
 
			if(nn==0) aa=aa/2;
			chn->transfer(txyz,t,r->id0+isnext,1);
			chn->transfertemp(thead,t,r->id0+isnext,1);
			if(aa==0) break;
		}

		if(TRES.logg>3) cerr<<"rotate t:"<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;		
		aa=5;
		while(1&&t&&r->id-t->id==1) {
			if(dmin<0.01) break;
			int a,b;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			/*for(c=-1;c<0;c++)*/ {
				float e,f;
				e=a*aa;f=b*aa;
				chn->transfer(txyz,t,r->id0+isnext,1);
				chn->transfertemp(thead,t,r->id0+isnext,1);
				g1=r->isatmid(1);
				g2=r->isatmid(2);
				g3=r->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;				
				rot.rotate(g1,e,r->id0+isnext,1);
				rot.rotate(g2,f,r->id0+isnext,1);
				float d;
				d=getrmsdmore(t,r,1);
				if(d<dmin) {
					if(fabs(dmin-d)>0.01) nn++;
					dmin=d;
					chn->transfer(txyz,t,r->id0+isnext,0);
					chn->transfertemp(thead,t,r->id0+isnext,0);		
				}
			}
			if(nn==0) aa=aa/2;
						 
			chn->transfer(txyz,t,r->id0+isnext,1);
			chn->transfertemp(thead,t,r->id0+isnext,1);
			if(aa==0) break;
		}
		chn->transfer(txyz,t,r->id0+isnext,1);
		chn->transfertemp(thead,t,r->id0+isnext,1);	
	 
		if(TRES.logg>3) cerr<<"rotate r: "<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
 
		if(fabs(dmin-dmin0)>0.05) {
			dmin0=dmin;
			goto re200;
		}

		chn->transfer(txyz,t,r->id0+isnext,1);
		chn->transfertemp(thead,t,r->id0+isnext,1);	
			
		//rotate the bond so that it is more closer to the exact
		//start exacting the bond to the original
		Res *r0=atmat[r->id0];
		Atm *a0,*a;
		if(exact&&r0&&r0->next&&r0->next->id-r0->id==1&&t&&r->id-t->id==1) {
			 
			a0=r0->isatmid(1);					 
			a=r->isatmid(1);
			 
			
			if(a0&&a) {
				float d=TRES.distance(a0,a);
				if(d<0.1) {
					float dx=(a0->xyz[0]-a->xyz[0]);
					float dy=(a0->xyz[1]-a->xyz[1]);
					float dz=(a0->xyz[2]-a->xyz[2]);
					Atm *aa;
					for(aa=r->atm;aa;aa=aa->next) {
						aa->xyz[0]+=dx;
                                		aa->xyz[1]+=dy;
                                		aa->xyz[2]+=dz;
					}
					if(r->next) {	
						for(aa=r->next->atm;aa;aa=aa->next) {
        						aa->xyz[0]+=dx;
        						aa->xyz[1]+=dy;
        						aa->xyz[2]+=dz;
						}	
						int i;
						for(i=0;i<3;i++) {
							r->next->temp[i*3+0]+=dx;
							r->next->temp[i*3+1]+=dy;
							r->next->temp[i*3+2]+=dz;
						}
					}							 
				}
			}
		}
		//end
 
		if(TRES.logg>3) pdb->write("ss");
				
		t=r;			
	}	
}


void PdbFix::mynewreversefix(Chn *chn,Res *rr0,int rend) {

	Rotate rot;

	chn->setlastresorder();
	
	Res *r0=chn->lastres();	
	Res *r;
	
	Chn *rotamer=TRES.findrotamer("backbone")->chn;
	if(rotamer==0) {
		cerr<<"backbone rotamer not found..."<<endl;
		exit(0);
	}

	Res *t=rr0;	
	for(r=rr0;r;r=r->last) {
		if(r->id0>rend) break;
		float end=0;
		if(t==0||t->last==0||t->last->id-t->id!=-1) {end=1; t=r; continue;}
		if(r->id-t->id!=1) {end=1; t=r; continue;}
		Res *r1=rotamer->isres(r->name,0);
		Res *r2;
		
		float dmin=100000;
		
		int cis=0;
		if(t&&atmat[t->id0]&&atmat[r->id0]) {
			Atm *ca0=atmat[t->id0]->isatmid(1);
			Atm *ca1=atmat[r->id0]->isatmid(1);
			if(ca0&&ca1) {
				float e=TRES.distance(ca0,ca1);
				if(e<3) cis=1;
			}
		}	
		 
 
		float txyz[200];
		float thead[200];

		//save t conformation
		if(r->next) rot.link(r->next,r,r->next->id0,1);		
	
		chn->transfer(txyz,t,r->id0+1,0);
		chn->transfertemp(thead,t,r->id0+1,0);				 
		float ang1=r->isatmid(1)->gettorsionangle();
		float ang2=r->isatmid(2)->gettorsionangle();
		
		for(r2=r1->more;r2;r2=r2->more) {
			chn->transfer(txyz,t,r->id0+1,1);
			chn->transfertemp(thead,t,r->id0+1,1);		
			r->transfer(r2);
			r->copytemp(r2);							 
			 
			rot.link(r,t,r->id0,1);
			if(r->next) rot.link(r->next,r,r->next->id0,1);				
			if(cis) {					
				float rt=t->atm->gettorsionangle();						
				rot.rotateomega(r->atm,-rt,r->id0+1,1);						
			}	
			float rt;													
							
			//a
			Atm *a=r->isatmid(2);
			rt=a->gettorsionangle();
			rot.rotate(a,-rt+ang2,r->next->id0,1);
			//a
			a=r->isatmid(1);
			rt=a->gettorsionangle();
			rot.rotate(a,-rt+ang1,r->last->id0,1);			 
			// 
			float d;
			if(end) d=getendrmsd(r,t,1);  
			else    d=getrmsd(r,t,1);
			if(d<dmin) {
				dmin=d;					
				chn->transfer(txyz,t,r->id0+1,0);
				chn->transfertemp(thead,t,r->id0+1,0);							
			}
		}		
			
		
		if(TRES.logg>3) cerr<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
		
		
		chn->transfer(txyz,t,r->id0+1,1);
		chn->transfertemp(thead,t,r->id0+1,1);	
		
		 		
		//rotate continously to find the local minima
 		
		if(end) dmin=getendrmsd(r,t,1);
		else	dmin=getrmsd(r,t,1);

		Atm *g1,*g2,*g3;
		 		 
		re200:
		float dmin0=dmin;	

		if(TRES.logg>3) cerr<<"rotate t:"<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;		
		int aa=5;
		while(1&&t&&r->id-t->id==1) {
			if(dmin<0.01) break;
			int a,b;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			{
				float e,f;
				e=a*aa;f=b*aa;
				chn->transfer(txyz,t,r->id0+1,1);
				chn->transfertemp(thead,t,r->id0+1,1);
				g1=r->isatmid(2);
				g2=r->isatmid(1);
				 
				if(g1==0||g2==0||g3==0) continue;				
				rot.rotate(g1,e,r->id0+1,1);
				rot.rotate(g2,f,r->id0+1,1);
				float d;
				if(end) d=getendrmsd(t,r,1); //1, all fixed atm, 0 this atom	
				else    d=getrmsd(t,r,1);
				if(d<dmin) {
					dmin=d;
					chn->transfer(txyz,t,r->id0+1,0);
					chn->transfertemp(thead,t,r->id0+1,0);
					nn++;
				}
			}
			if(nn==0) aa=aa/2;
						 
			chn->transfer(txyz,t,r->id0+1,1);
			chn->transfertemp(thead,t,r->id0+1,1);
			if(aa==0) break;
		}
		
		chn->transfer(txyz,t,r->id0+1,1);
		chn->transfertemp(thead,t,r->id0+1,1);	
		
		if(TRES.logg>3) cerr<<"rotate r: "<<r->name<<r->id0<<"  minimal rmsd: "<<dmin<<endl;
	 

		aa=5;
		while(1&&r&&r->id-t->id==1) {
			if(dmin<0.01) break;
			int a,b,c;
			int nn=0;
			for(a=-1;a<2;a++)  	
			for(b=-1;b<2;b++)  
			for(c=-1;c<2;c++) {
				float e,f,g;
				e=a*aa;f=b*aa;g=c*aa;				
				chn->transfer(txyz,t,r->id0+1,1);
				chn->transfertemp(thead,t,r->id0+1,1);
				g1=t->isatmid(2);
				g2=t->isatmid(1);
				g3=t->isatmid(0);
				if(g1==0||g2==0||g3==0) continue;	
							
				rot.rotate(g1,e,r->id0+1,1);
				rot.rotate(g2,f,r->id0+1,1);	
				rot.rotateomega(g3,g,r->id0+1,1);
				 
				float d;
				if(end) d=getendrmsd(t,r,1); //1, all fixed atm, 0 this atom	
				else    d=getrmsd(t,r,1);
				if(d<dmin) {
					dmin=d;
					chn->transfer(txyz,t,r->id0+1,0);
					chn->transfertemp(thead,t,r->id0+1,0);
					nn++;
				}
			}
 
			if(nn==0) aa=aa/2;
			chn->transfer(txyz,t,r->id0+1,1);
			chn->transfertemp(thead,t,r->id0+1,1);
			if(aa==0) break;
		}

		
		if(fabs(dmin-dmin0)>0.05) {
			dmin0=dmin;
			goto re200;
		}
		//end
		chn->transfer(txyz,t,r->id0+1,1);
		chn->transfertemp(thead,t,r->id0+1,1);	
			
		//rotate the bond so that it is more closer to the exact
		//start exacting the bond to the original
		Res *r0=atmat[r->id0];
		Atm *a0,*b0,*a,*b;
		if(exact&&r0&&r0->last&&r0->last->id-r0->id==-1&&t&&t->id-r->id==-1) {

			b0=r0->last->isatmid(1);
			a0=r0->isatmid(1);	
	
			b=r->last->isatmid(1);
			a=r->isatmid(1);			
			if(a0&&b0&&a&&b) {
				float d=TRES.distance(a0,b0);
				if(d>3.9) d=3.9;
				else if(d<3.6&&d>3.0) d=3.6;
				float d0=TRES.distance(a,b);
				float x,y,z;
                        	float dx,dy,dz;
				dx=(a->xyz[0]-b->xyz[0])/d0;
                        	dy=(a->xyz[1]-b->xyz[1])/d0;
                        	dz=(a->xyz[2]-b->xyz[2])/d0;
                        	x=(d-d0)*dx;
                        	y=(d-d0)*dy;
                        	z=(d-d0)*dz;	
				Atm *aa;
				for(aa=r->atm;aa;aa=aa->next) {
					aa->xyz[0]+=x;
                                	aa->xyz[1]+=y;
                                	aa->xyz[2]+=z;
				}
				int i;
				for(i=0;i<3;i++) {
					r->temp[i*3+0]+=x;
					r->temp[i*3+1]+=y;
					r->temp[i*3+2]+=z;
				}
				cerr<<d<<" "<<d0<<" "<<TRES.distance(a,b)<<endl;
			}			
		}
		//end


		if(r->last) rot.link(r->last,r,r->last->id0,0);
		if(t&&end==1&&r->last)  superimposesegment(r->last,t,1);
		else if(t&&end==1) superimposesegment(r,t,1);
		if(TRES.logg>3) pdb->write("ss");		
		
		t=r;			
	}	
}
 
float PdbFix::energy(Res *r,Res *t) {
	
	Res *rr;
	float d=0;
	int n=0;
	for(rr=r;rr;rr=rr->next) {
		if(rr->id0>t->id0) break;
		Res *rr0=atmat[rr->id0];
		if(rr0==0) continue;
		Atm *a;
		for(a=rr0->atm;a;a=a->next) {
			Atm *b;
			b=rr->isatmid(a->tatm->id);
			if(b==0) continue;
			d+=TRES.distsqr(a->xyz,b->xyz);
			n++;
		}
	}
	if(n==0) n=1;
	return sqrt(d/n);
}

Res **PdbFix::maxenergystem(Res *r,Res *t) {
	
	Res *rr,*tt;
	float dmin=-1000;	
	Res **stem=new Res *[2];
	stem[0]=0;
	stem[1]=0;
	for(rr=r;rr;rr=rr->next) {
		if(rr->id0>t->id0) break;
		tt=rr->chn->isres0(rr->id0+10);
		if(tt==0||tt->id0>t->id0) tt=t;
		float d=energy(rr,tt);
		if(d>dmin) {
			dmin=d;
			stem[0]=rr;
			stem[1]=tt;
		} 
		if(tt->id0>=t->id0) break; 
	}
	return stem;	
}
void PdbFix::minimize0(Chn *chn) {

        ProFix segen;
        segen.pdb=pdb;
        segen.flex=1;
        segen.arbt=100;
        segen.cid=chn->id;
        segen.randcoil=2;
        segen.smoothclash=1;
        segen.part=0;
        segen.pdborg=pdborg;
        segen.onlyenergy=0;
	segen.cycle=1000;
	segen.flat=111;
	segen.step=0.17;
	segen.outlet=0.4;
        Res *r;
        for(r=chn->res;r;r=r->next) {
                Res **stem=findsegment(r,5);
                if(stem==0) continue;
		float d=energy(stem[0],stem[1]);
		if(d<0.5) {
			if(stem) delete [] stem;stem=0;
			continue;
		}
                float *xyzorg=chn->gettransfer(stem[0],stem[1]->id0);
		
                float *xyzout=segen.fixpdbsegment(stem);
		
                if(xyzout) {
                        chn->transfer(xyzout,stem[0],stem[1]->id0,1);
                }
                else {
                        chn->transfer(xyzorg,stem[0],stem[1]->id0,1);
                }
		
                if(stem) delete [] stem;stem=0;		
        }
        segen.pdb=0;
}


int PdbFix::isbreaker(Res *r1,Res *r2){
	Res *rr;
 
	for(rr=r1;rr;rr=rr->next) {
		if(rr->id0>r2->id0) break;
		if(rr->next&&rr->next->id-rr->id!=1) return 1;
		if(rr->last&&rr->id-rr->last->id!=1) return 1;			
	}
	return 0;
}
int PdbFix::islinked(Res *r1,Res *r2){
	Res *rr;
 
	for(rr=r1;rr;rr=rr->next) {
		if(rr->id0>r2->id0) break;
		if(rr->next) {
			int n=rr->chn->islinked(rr,rr->next);
			if(n==0) return 0;		
		}
		else if(rr->last) {
			int n=rr->chn->islinked(rr->last,rr);
			if(n==0) return 0;		
		}		
	}
	return 1;
}
void PdbFix::minimize(Chn *chn) {

	ProFix segen;
	segen.pdb=pdb;	
        segen.flex=1;
        segen.arbt=100;
        segen.cid=chn->id;
        segen.randcoil=2;
        segen.smoothclash=1;
        segen.part=0;
	segen.pdborg=pdborg;
	segen.onlyenergy=0;
	segen.cycle=1000;
	segen.flat=111;
	segen.step=0.17;
	segen.outlet=0.4;
	 
	Res *r;
	for(r=chn->res;r;r=r->next) {
		Res **stem=findsegment(r,4);
		
		if(stem==0) continue;
		Res *r1,*r2;
		r1=stem[0];r2=stem[1];
		if(isbreaker(r1,r2)==1) continue;
		float d=energy(stem[0],stem[1]);
		if(islinked(r1,r2)==1) {
			segen.addoriginal=1;
		}
		else {
			segen.addoriginal=0;
		}
		if(d<0.2&&segen.addoriginal) {
			float *xyzorg=chn->gettransfer(stem[0],stem[1]->id0);
			segen.readyminfix(stem[0],stem[1]->id0);
			strcpy(segen.boundtag,"P");
			strcpy(segen.dforce,"Ddu");
			segen.setpdbfixbound(2);
			segen.minimizefix(xyzorg,1);
			if(segen.bound) delete segen.bound;segen.bound=0;
			chn->transfer(xyzorg,stem[0],stem[1]->id0,1);				
			if(stem) delete [] stem;stem=0;	
			r1=chn->isres0((r1->id0+r2->id0)/2);
			if(r1) r=r1;	
		}
		else {
			float *xyzorg=chn->gettransfer(stem[0],stem[1]->id0);		
                	float *xyzout=segen.fixpdbsegment(stem);		
                	if(xyzout) {
                        	chn->transfer(xyzout,stem[0],stem[1]->id0,1);
               	 	}
                	else {
                        	chn->transfer(xyzorg,stem[0],stem[1]->id0,1);
                	}		
                	if(stem) delete [] stem;stem=0;	
		}
		chn->header(r1->id0-3,r2->id0+3);	
	}
	segen.pdb=0;	 
}


void PdbFix::minimize(Res *first,Res *last) {
	Chn *chn=first->chn;
	ProFix segen;
	segen.pdb=pdb;	
        segen.flex=1;
        segen.arbt=150;
        segen.cid=chn->id;
        segen.randcoil=2;
        segen.smoothclash=1;
        segen.part=0;
	segen.pdborg=pdborg;
	segen.onlyenergy=0;
	segen.cycle=1000;
	segen.flat=111;
	segen.step=0.17;
	segen.outlet=0.4;
	 
	Res *r;
	for(r=first;r;r=r->next) {		
		if(r->id0>last->id0) break;
		Res **stem=findsegment(r,3);
		if(stem==0) continue;
		if(stem[0]->id0<first->id0) stem[0]=first;
		if(stem[1]->id0>last->id0) stem[1]=last;
		Res *r1,*r2;
		r1=stem[0];r2=stem[1];
		
		if(isbreaker(r1,r2)==1) {
			if(stem) delete [] stem;stem=0;
			continue;
		}
		float d=energy(stem[0],stem[1]);
		if(islinked(r1,r2)==1) {
			segen.addoriginal=1;
		}
		else {
			segen.addoriginal=0;
		}	
			
		if(d<0.2&&segen.addoriginal) {
			float *xyzorg=chn->gettransfer(stem[0],stem[1]->id0);
			segen.readyminfix(stem[0],stem[1]->id0);
			strcpy(segen.boundtag,"P");
			strcpy(segen.dforce,"Ddu");
			segen.setpdbfixbound(2);
			segen.minimizefix(xyzorg,1);
			if(segen.bound) delete segen.bound;segen.bound=0;
			chn->transfer(xyzorg,stem[0],stem[1]->id0,1);				
			if(stem) delete [] stem;stem=0;	
			r1=chn->isres0((r1->id0+r2->id0)/2);
			if(r1) r=r1;
			if(xyzorg) delete [] xyzorg;xyzorg=0;
		}
		else {
			float *xyzorg=chn->gettransfer(stem[0],stem[1]->id0);	
			float *xyzorgtemp=chn->gettransfertemp(stem[0],stem[1]->id0);		
                	float *xyzout=segen.fixpdbsegment(stem);		
                	if(xyzout) {
                        	chn->transfer(xyzout,stem[0],stem[1]->id0,1);
				chn->header(r1,r2->id0);
               	 	}
                	else {
                        	chn->transfer(xyzorg,stem[0],stem[1]->id0,1);
				chn->transfertemp(xyzorgtemp,stem[0],stem[1]->id0,1);
                	}		
			if(xyzorg) delete [] xyzorg;xyzorg=0;
			if(xyzorgtemp) delete [] xyzorgtemp;xyzorgtemp=0;
			if(xyzout) delete [] xyzout;xyzout=0;
                	if(stem) delete [] stem;stem=0;	
		}
		if(stem) delete [] stem;stem=0;	
	}
	segen.pdb=0;	 
}

int PdbFix::directlink(Res *first,Res *last,int fast0) {

	Chn *chn=first->chn;
	ProFix segen;
	segen.pdb=pdb;	
        segen.flex=1;
        segen.arbt=150;
        segen.cid=chn->id;
        segen.randcoil=3;
        segen.smoothclash=1;
        segen.part=0;
	segen.pdborg=pdborg;
	segen.onlyenergy=0;
	segen.cycle=1000;
	segen.flat=111;
	segen.step=0.17;
	segen.outlet=0.4;	 
	segen.addoriginal=0;
	segen.dconst=0;
	segen.rbackbone=1;
	if(fast0==1) segen.arbt=50;
	Res *stem[2];

	stem[0]=first;stem[1]=last;
	float *xyzorg=chn->gettransfer(stem[0],stem[1]->id0);	
	float *xyzorgtemp=chn->gettransfertemp(stem[0],stem[1]->id0);		
    	float *xyzout=segen.fixpdbsegment(stem);	

	int nt=1;
	if(xyzout) {
              	chn->transfer(xyzout,stem[0],stem[1]->id0,1);
		chn->header(first,last->id0);
      	}
    	else {
         	chn->transfer(xyzorg,stem[0],stem[1]->id0,1);
		chn->transfertemp(xyzorgtemp,stem[0],stem[1]->id0,1);	
		nt=0;
    	}		
		
	
	segen.pdb=0;	 
	if(xyzout) delete [] xyzout; xyzout=0;
	if(xyzorg) delete [] xyzorg; xyzorg=0;
	if(xyzorgtemp) delete [] xyzorgtemp; xyzorgtemp=0;
	return nt;
}



void PdbFix::linkallres(Res *first,Res *last) {
	Chn *chn=first->chn;
	ProFix segen; 
	segen.pdb=pdb;	
        segen.flex=1;
        segen.arbt=150;
        segen.cid=chn->id;
        segen.randcoil=3;
        segen.smoothclash=1;
        segen.part=0;
	segen.pdborg=pdborg;
	segen.onlyenergy=0;
	segen.cycle=1000;
	segen.flat=111;
	segen.step=0.17;
	segen.outlet=0.4;
	segen.rbackbone=1;
	segen.addoriginal=0;
	segen.dconst=0;
	 
	Res *r;
	for(r=first;r;r=r->next) {		
		if(r->id0>last->id0) break;
		if(r->next==0) continue;
		if(r->next->id-r->id!=1) continue; //breaker	
		if(chn->isdatabaselinked(r,r->next)==1) continue; //not linked
		//if(r->last&&r->id-r->last->id!=1) continue;
		//if(r->next&&r->next->id-r->id!=1) continue; 
		Res **stem=new Res*[2];
		stem[0]=r;
		stem[1]=r->next;
		
		while(1) {
			//segen.addoriginal=0;	
			float *xyzorg=chn->gettransfer(stem[0],stem[1]->id0);	
			float *xyzorgtemp=chn->gettransfertemp(stem[0],stem[1]->id0);		
                	float *xyzout=segen.fixpdbsegment(stem);		
                	if(xyzout) {
                        	chn->transfer(xyzout,stem[0],stem[1]->id0,1);				
				if(xyzorg) delete [] xyzorg;xyzorg=0;
				if(xyzout) delete [] xyzout;xyzout=0;
				if(xyzorgtemp) delete [] xyzorgtemp;xyzorgtemp=0;	
				chn->header(stem[0],stem[1]->id0);							 
				break;
               	 	}
                	else if(stem[1]->id-stem[0]->id<15) {
                        	chn->transfer(xyzorg,stem[0],stem[1]->id0,1);
				chn->transfertemp(xyzorgtemp,stem[0],stem[1]->id0,1);
				if(stem[0]->last) {
					stem[0]=stem[0]->last;
				} 
				if(stem[1]->next) {
					stem[1]=stem[1]->next;
				} 	
				if(xyzorg) delete [] xyzorg;xyzorg=0;	
				if(xyzorgtemp) delete [] xyzorgtemp;xyzorgtemp=0;
                	}	
			else {
				cerr<<"could not link the residues: "<<r->name<<r->oid<<"--"<<r->next->name<<r->next->oid<<endl;
				chn->transfer(xyzorg,stem[0],stem[1]->id0,1);
				chn->transfertemp(xyzorgtemp,stem[0],stem[1]->id0,1);
				if(xyzorg) delete [] xyzorg;xyzorg=0;
				if(xyzorgtemp) delete [] xyzorgtemp;xyzorgtemp=0;
				break;
			}	
		}
		if(stem) delete [] stem;stem=0;						
	}
	segen.pdb=0;	 
}

float PdbFix::getendrmsd(Res *r,Res *t,int flg) {

	Res *r1,*t1,*r2;
	Chn *chn0;
	chn0=pdborg->ischain(r->chn->id);
	if(chn0==0) {
		cerr<<"unexpected errors!"<<endl;
		return 0;
	}
	r1=0;t1=0;r2=0;
	//if(r)  r1=chn0->isresoid(r->oid);
	//if(t)  t1=chn0->isresoid(t->oid);
	if(r) r1=atmat[r->id0];
	if(t) t1=atmat[t->id0];	 	 
	if(r1) r2=r1->last;
	if(r1==0) {
		cerr<<"unexpected errors!"<<endl;
		return 0;
	}
	
	if(r2&&r1->id-r2->id!=1) r2=0;
	
	
	float rmsd=0;
	int   ntot=0;

	
	Atm *a,*b;
	Atm *a0,*b0;

	//distance within r1
	for(a0=r1->atm;a0;a0=a0->next) {
		if(flg==0&&a0->tatm->id!=1) continue;
		if(flg==1&&a0->tatm->name[1]!='H'&&a0->tatm->id>4) continue;
		if(flg==1&&a0->tatm->name[1]=='H'&&a0->tatm->bond[0]->id>0) continue;
		a=r->isatmid(a0->tatm->id);
		if(a==0) continue;
		for(b0=a0->next;b0;b0=b0->next) {			
			//only ca			
			if(flg==0&&b0->tatm->id!=1) continue;
			//only backbone including cb			
			if(flg==1&&b0->tatm->name[1]!='H'&&b0->tatm->id>4) continue;
			if(flg==1&&b0->tatm->name[1]=='H'&&b0->tatm->bond[0]->id>0) continue;
			//all atoms;
			b=r->isatmid(b0->tatm->id);
			if(b==0) continue;
			float d1=TRES.distance(a0,b0);
			float d2=TRES.distance(a,b);
			rmsd+=(d1-d2)*(d1-d2);
			ntot++;
		}
	}
	
	//distance between r1 and t1
	for(a0=r1->atm;a0;a0=a0->next) {
		if(t1==0) break;	
		if(flg==0&&a0->tatm->id!=1) continue;
		if(flg==1&&a0->tatm->name[1]!='H'&&a0->tatm->id>4) continue;
		if(flg==1&&a0->tatm->name[1]=='H'&&a0->tatm->bond[0]->id>0) continue;
		a=r->isatmid(a0->tatm->id);
		if(a==0) continue;
		for(b0=t1->atm;b0;b0=b0->next) {			
			//only ca			
			if(flg==0&&b0->tatm->id!=1) continue;
			//only backbone including cb			
			if(flg==1&&b0->tatm->name[1]!='H'&&b0->tatm->id>4) continue;
			if(flg==1&&b0->tatm->name[1]=='H'&&b0->tatm->bond[0]->id>0) continue;
			//all atoms;
			b=t->isatmid(b0->tatm->id);
			if(b==0) continue;
			float d1=TRES.distance(a0,b0);
			float d2=TRES.distance(a,b);
			rmsd+=(d1-d2)*(d1-d2);
			ntot++;
		}
	}

	if(r2==0) {
		if(ntot==0)ntot=1;
		return sqrt(rmsd/ntot);
	}

	//distance between r2 tail and r1
	for(a0=r2->atm;a0;a0=a0->next) {
		if(flg==0&&a0->tatm->id!=1) continue;
		if(flg==1&&a0->tatm->name[1]!='H'&&a0->tatm->id>4) continue;
		if(flg==1&&a0->tatm->name[1]=='H'&&a0->tatm->bond[0]->id>0) continue;
		if(a0->tatm->id<1||a0->tatm->id>3) continue;
		float *xyz;
		if(a0->tatm->id==1) xyz=r->temp;
		else if(a0->tatm->id==2) xyz=r->temp+3;
		else if(a0->tatm->id==3) xyz=r->temp+6;
		else continue;
		//a=r->isatmid(a0->tatm->id);
		//if(a==0) continue;
		for(b0=r1->atm;b0;b0=b0->next) {			
			//only ca			
			if(flg==0&&b0->tatm->id!=1) continue;
			//only backbone including cb			
			if(flg==1&&b0->tatm->name[1]!='H'&&b0->tatm->id>4) continue;
			if(flg==1&&b0->tatm->name[1]=='H'&&b0->tatm->bond[0]->id>0) continue;
			//all atoms;
			b=r->isatmid(b0->tatm->id);
			if(b==0) continue;
			float d1=TRES.distance(a0,b0);
			float d2=TRES.distance(xyz,b->xyz);
			rmsd+=(d1-d2)*(d1-d2);
			ntot++;
		}
	}
	if(t1==0||t==0) {
		if(ntot==0)ntot=1;
		return sqrt(rmsd/ntot);
	}
	//distance between r2 tail and t1;
	for(a0=r2->atm;a0;a0=a0->next) {
		if(flg==0&&a0->tatm->id!=1) continue;
		if(flg==1&&a0->tatm->name[1]!='H'&&a0->tatm->id>4) continue;
		if(flg==1&&a0->tatm->name[1]=='H'&&a0->tatm->bond[0]->id>0) continue;
		if(a0->tatm->id<1||a0->tatm->id>3) continue;
		float *xyz;
		if(a0->tatm->id==1) xyz=r->temp;
		else if(a0->tatm->id==2) xyz=r->temp+3;
		else if(a0->tatm->id==3) xyz=r->temp+6;
		//a=r->isatmid(a0->tatm->id);
		//if(a==0) continue;
		for(b0=t1->atm;b0;b0=b0->next) {			
			//only ca			
			if(flg==0&&b0->tatm->id!=1) continue;
			//only backbone including cb			
			if(flg==1&&b0->tatm->name[1]!='H'&&b0->tatm->id>4) continue;
			if(flg==1&&b0->tatm->name[1]=='H'&&b0->tatm->bond[0]->id>0) continue;
			//all atoms;
			b=t->isatmid(b0->tatm->id);
			if(b==0) continue;
			float d1=TRES.distance(a0,b0);
			float d2=TRES.distance(xyz,b->xyz);
			rmsd+=(d1-d2)*(d1-d2);
			ntot++;
		}
	}

	if(ntot==0)ntot=1;
	return sqrt(rmsd/ntot);	 
}

float PdbFix::getrmsd(Res *r,Res *t,int flg) {

	Res *r1,*t1,*r2;
	Chn *chn0;
	chn0=pdborg->ischain(r->chn->id);
	if(chn0==0) {
		cerr<<"unexpected errors!"<<endl;
		return -1;
	}
	r1=atmat[r->id0];
	t1=atmat[t->id0];
	//r1=chn0->isresoid(r->oid);
	//t1=chn0->isresoid(t->oid);	 
	r2=r1->last;
	if(r2&&r1->id-r2->id!=1) r2=0;
	if(r1==0||t1==0) {
		cerr<<"unexpected errors!"<<endl;
		return 0;
	}
	
	float rmsd=0;
	int   ntot=0;

	//r1
	Atm *a,*b;

	for(a=r1->atm;a;a=a->next) {
		if(flg==0&&a->tatm->id!=1) continue;
		if(flg==1&&a->tatm->name[1]!='H'&&a->tatm->id>4) continue;
		if(flg==1&&a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>1) continue; 
		b=r->isatmid(a->tatm->id);
		if(b==0) continue;
		rmsd+=TRES.distsqr(b->xyz,a->xyz);
		ntot++;
	}

	if(r2==0) {
		if(ntot==0) ntot=1;
		return sqrt(rmsd/ntot);
	} 
	
	//r2
	for(a=r2->atm;a;a=a->next) {
		if(flg==0&&a->tatm->id!=1) continue;
		if(flg==1&&a->tatm->name[1]!='H'&&a->tatm->id>4) continue;
		if(flg==1&&a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>1) continue; 
		if(a->tatm->id<1||a->tatm->id>3) continue;
		float *xyz;
		if(a->tatm->id==1) xyz=r->temp;
		else if(a->tatm->id==2) xyz=r->temp+3;
		else if(a->tatm->id==3) xyz=r->temp+6;
		else continue;			
		rmsd+=TRES.distsqr(xyz,a->xyz);
		ntot++;
	}

	if(ntot==0) ntot=1;
	return sqrt(rmsd/ntot);
}

float PdbFix::getrmsdmore(Res *t,Res *r,int flg) {

	Res *r1,*t1,*r2;
	Chn *chn0;
	chn0=pdborg->ischain(r->chn->id);
	if(chn0==0) {
		cerr<<"unexpected errors!"<<endl;
		return -1;
	}
	
	t1=atmat[t->id0];
	r1=atmat[r->id0];
	 
 
	r2=r1->next;
	if(!(r2&&r2->id-r1->id==1))  r2=0;
 

	if(r1==0||t1==0) {
		cerr<<"unexpected errors!"<<endl;
		return 0;
	}
	
	float rmsd=0;
	int   ntot=0;

	//r1
	Atm *a,*b;

	for(a=r1->atm;a;a=a->next) {
		if(flg==0&&a->tatm->id!=1) continue;
		if(flg==1&&a->tatm->name[1]!='H'&&a->tatm->id>4) continue;
		if(flg==1&&a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>1) continue; 
		b=r->isatmid(a->tatm->id);
		if(b==0) continue;
		rmsd+=TRES.distsqr(b->xyz,a->xyz);
		ntot++;
	}

	if(r2==0) {
		if(ntot==0) ntot=1;
		return sqrt(rmsd/ntot);
	} 
	
	//r2
	for(a=r2->atm;a;a=a->next) {
		if(flg==0&&a->tatm->id!=1) continue;
		if(flg==1&&a->tatm->name[1]!='H'&&a->tatm->id>1) continue;
		if(flg==1&&a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>1) continue; 	
		if(r->next==0) continue;
		if(r->next->id-r->id!=1) continue;	 
		b=r->next->isatmid(a->tatm->id);
		if(b==0) continue; 			
		rmsd+=TRES.distsqr(b->xyz,a->xyz);
		ntot++;
	}

	if(ntot==0) ntot=1;
	return sqrt(rmsd/ntot);
}


void PdbFix::superimposesegment(Res *r1,Res *r2,int flg) {
	 
	Res *r,*rr;
	Stralg algn;
	algn.flag=1;
	//algn.onlyaligned=1;
	Res **alga=new Res*[r2->id0-r1->id0+100];
	Res **algb=new Res*[r2->id0-r1->id0+100];
	
	Chn *chn0=pdborg->ischain(r1->chn->id);
	if(chn0==0) {
		cerr<<"unexpected error"<<endl;
		exit(0);
	}
	int n=0;
	int jj=0;
	for(r=r1;r;r=r->next) {	
		if(r->id0>r2->id0) break;	
		alga[n]=r;
		rr=chn0->isresoid(r->oid);
		if(rr==0) {
			cerr<<"program detects error in constructing missing atoms..."<<endl;
			exit(0);
		}
		algb[n]=rr;
		Atm *aa;
		for(aa=rr->atm;aa;aa=aa->next){
			if(flg==1&&aa->tatm->name[1]!='H'&&aa->tatm->id>=4) continue;
			if(flg==1&&aa->tatm->name[1]=='H'&&aa->tatm->bond[0]->id>=3) continue;
			jj++;
		}
		n++;			
	}
	if(jj<3) algn.flag=1000;
	alga[n]=0;
	algb[n]=0;
	algn.alga=alga;
	algn.algb=algb;
	if(TRES.logg>3)pdb->write("sb");
	float dd=algn.superimpose(flg);
	algn.alga=0;
	algn.algb=0;
	if(alga) delete [] alga; alga=0;
	if(algb) delete [] algb; algb=0;
	//write
	if(TRES.logg>3) {
		char cid=r1->chn->id;
		r1->chn->id='A';
		chn0->id='B';
		FILE *fp=fopen("s","w");
		r1->chn->write(fp);
		chn0->write(fp);
		r1->chn->id=cid;
		chn0->id=cid;
		fclose(fp);
	}
	 
}

Res **PdbFix::findsegment(Res *r,int n) {

	Res *r1;
 
	for(r1=r;r1;r1=r1->next) {
		if(r1->next==0) break;
		if(r1->next->id-r1->id!=1) break;
		if(r1->next->id0-r->id0>=n) break;  
	}
	 
	Res **stem=new Res*[2];
	stem[0]=r;
	stem[1]=r1;
	return stem;
}

void PdbFix::myfix(Chn *chn) {

	chn->setlastresorder();
	
	Res *r;
	for(r=chn->res;r;r=r->next) {
		Res **stem=findsegment(r,4);
		if(stem==0) continue;
		Res *r1=stem[0];
		Res *r2=stem[1];
		cerr<<r1->name<<r1->id0<<" "<<r2->name<<r2->id0<<endl;
		fixsegment(stem);
		chn->start=1;
		chn->write("sa.pdb");	
		delete [] stem;stem=0;
		r=r2;
	}
}

void PdbFix::fixsegment(Res **stem) {
 
	//create chain
	
	Res *r1,*r2;
	Res *r,*rr;
	Atm *a,*aa;
	r1=stem[0];
	r2=stem[1];	
	Chn *s=new Chn;
	s->create(r1,r2->id0);
	s->id=r1->chn->id;
	
	Pdb temp;
	temp.chn=s;
	s->pdb=&temp;
	s->configure();
	s->start=r1->chn->start;
	Chn *chn0=pdborg->ischain(s->id);
	Chn *chn=r1->chn;

	Res *t;
	for(t=s->res;t;t=t->next) {	
		r=chn0->isresoid(t->oid);
		t->id=r->id;
		t->id0=r->id0;
	}
	//align the segment
	Stralg algn;
	algn.flag=1;
	Res **alga=new Res*[r2->id0-r1->id0+100];
	Res **algb=new Res*[r2->id0-r1->id0+100];

	int j;
	for(j=0;j<r2->id0-r1->id0+100;j++) {
		alga[j]=0;
		algb[j]=0;
	}	
	
	

	
	int n=0;
	for(r=s->res;r;r=r->next) {		
		alga[n]=r;
		rr=chn0->isresoid(r->oid);
		if(rr==0) {
			cerr<<"program detects error in constructing missing atoms..."<<endl;
			exit(0);
		}
		algb[n]=rr;
		n++;
		if(n==2) break;
	}
	alga[n]=0;
	algb[n]=0;
	algn.alga=alga;
	algn.algb=algb;
	int i;
	int m=0;
	
	for(i=0;i<n;i++) {
		a=alga[i]->isatmid(1);
		aa=algb[i]->isatmid(1);	
		if(a&&aa) {
			m++;
		}
	}
	if(m<2) m=1;
	else	m=0;
	if(TRES.logg>3) s->write("sb");
	float dd=algn.superimpose(m);
	algn.alga=0;
	algn.algb=0;
	if(alga) delete [] alga; alga=0;
	if(algb) delete [] algb; algb=0;
	//write
	if(TRES.logg>3) {
		char cid=s->id;
		s->id='A';
		chn0->id='B';
		FILE *fp=fopen("s","w");
		s->write(fp);
		chn0->write(fp);
		s->id=cid;
		chn0->id=cid;
		fclose(fp);
	}
	if(r2->id0-r1->id0+1<=2) {
		float *xyz=s->gettransfer(s->res,r2->id0);
		r1->chn->transfer(xyz,r1,r2->id0,1);
		if(xyz) delete [] xyz;xyz=0;
		return;
	} 	

	ProFix segen;
	segen.pdb=&temp;	
        segen.flex=1;
        segen.arbt=100;
        segen.cid=s->id;
        segen.randcoil=2;
        segen.smoothclash=1;
        segen.part=0;
	segen.pdborg=pdborg;
	segen.onlyenergy=0;
	temp.header();
	stem[0]=s->res->next;
	stem[1]=s->lastres();
	float *xyzout=segen.fixpdbsegment(stem);
	if(xyzout) {			
		r1=s->isresoid(stem[0]->oid);
		r2=s->isresoid(stem[1]->oid);
		s->transfer(xyzout,r1,r2->id0,1);
		if(xyzout) delete [] xyzout;xyzout=0;
		xyzout=s->gettransfer(s->res,10000);
		r1=chn->isresoid(s->res->oid);
		r2=chn->isresoid(s->lastres()->oid);
		chn->transfer(xyzout,r1,r2->id0,1);
		if(xyzout) delete [] xyzout;xyzout=0;
	}
	segen.pdb=0;
}


void PdbFix::myfix0(Chn *chn) {

	chn->setlastresorder();
	
	Res *r;
	for(r=chn->res;r;r=r->next) {
		Res **stem=findsegment(r,5);
		if(stem==0) continue;
		Res *r1=stem[0];
		Res *r2=stem[1];
		cerr<<r1->name<<r1->id0<<" "<<r2->name<<r2->id0<<endl;
		fixsegment0(stem);
		chn->start=1;
		chn->write("sa.pdb");	
		delete [] stem;stem=0;
		r=r2;
	}
}
 
void PdbFix::fixsegment0(Res **stem) {
	
	//create chain
	
	Res *r1,*r2;
	Res *r,*rr;
	Atm *a,*aa;
	r1=stem[0];
	r2=stem[1];	
	Chn *s=new Chn;
	s->create(r1,r2->id0);
	s->id=r1->chn->id;
	
	Pdb temp;
	temp.chn=s;
	s->pdb=&temp;
	s->configure();
	s->start=r1->chn->start;
	Chn *chn0=pdborg->ischain(s->id);
	Chn *chn=r1->chn;

	Res *t;
	for(t=s->res;t;t=t->next) {	
		r=chn0->isresoid(t->oid);
		t->id=r->id;
		t->id0=r->id0;
	}
	//align the segment
	Stralg algn;
	algn.flag=1;
	Res **alga=new Res*[r2->id0-r1->id0+100];
	Res **algb=new Res*[r2->id0-r1->id0+100];

	int j;
	for(j=0;j<r2->id0-r1->id0+100;j++) {
		alga[j]=0;
		algb[j]=0;
	}	
	
	

	
	int n=0;
	for(r=s->res;r;r=r->next) {		
		alga[n]=r;
		rr=chn0->isresoid(r->oid);
		if(rr==0) {
			cerr<<"program detects error in constructing missing atoms..."<<endl;
			exit(0);
		}
		algb[n]=rr;
		n++;
		if(n==2) break;
	}
	alga[n]=0;
	algb[n]=0;
	algn.alga=alga;
	algn.algb=algb;
	int i;
	int m=0;
	
	for(i=0;i<n;i++) {
		a=alga[i]->isatmid(1);
		aa=algb[i]->isatmid(1);	
		if(a&&aa) {
			m++;
		}
	}
	if(m<2) m=1;
	else	m=0;
	if(TRES.logg>3) s->write("sb");
	float dd=algn.superimpose(m);
	algn.alga=0;
	algn.algb=0;
	if(alga) delete [] alga; alga=0;
	if(algb) delete [] algb; algb=0;
	//write
	if(TRES.logg>3) {
		char cid=s->id;
		s->id='A';
		chn0->id='B';
		FILE *fp=fopen("s","w");
		s->write(fp);
		chn0->write(fp);
		s->id=cid;
		chn0->id=cid;
		fclose(fp);
	}
	if(r2->id0-r1->id0+1<=2) {
		float *xyz=s->gettransfer(s->res,r2->id0);
		r1->chn->transfer(xyz,r1,r2->id0,1);
		if(xyz) delete [] xyz;xyz=0;
		return;
	} 	

	ProFix segen;
	segen.pdb=&temp;	
        segen.flex=1;
        segen.arbt=100;
        segen.cid=s->id;
        segen.randcoil=2;
        segen.smoothclash=1;
        segen.part=0;
	segen.pdborg=pdborg;
	segen.onlyenergy=0;
	temp.header();
	stem[0]=s->res->next;
	stem[1]=s->lastres();
	float *xyzout=segen.fixpdbsegment(stem);
	if(xyzout) {			
		r1=s->isresoid(stem[0]->oid);
		r2=s->isresoid(stem[1]->oid);
		s->transfer(xyzout,r1,r2->id0,1);
		if(xyzout) delete [] xyzout;xyzout=0;
		xyzout=s->gettransfer(s->res,10000);
		r1=chn->isresoid(s->res->oid);
		r2=chn->isresoid(s->lastres()->oid);
		chn->transfer(xyzout,r1,r2->id0,1);
		if(xyzout) delete [] xyzout;xyzout=0;
	}
	segen.pdb=0;
}




//temp
void PdbFix::fixsave() {
  
	Chn *c,*c1;
	 
	c1=0;
	for(c=pdb->chn;c;c=c->next) {   
		Chn *c0=createchain(c);
		c0->transform(3);
		c0->setallnear(c0->res,10000);
		c0->pdb=pdb;
		//
		c0->next=c->next;
		c->next=0;		
		delete c;c=0;
		if(c1==0) pdb->chn=c0;
		else 	  c1->next=c0;
		pdb->configure(); 	
		c=c0;	 		 
		if(TRES.logg>3) pdb->write("temp.pdb");
		fix(c);
		c1=c; 
	}	
}

void PdbFix::fixold() {
  
	Chn *c,*c1;
	 
	c1=0;
	for(c=pdb->chn;c;c=c->next) {   
		Chn *c0=createchain(c);
		c0->transform(3);
		c0->setallnear(c0->res,10000);
		c0->pdb=pdb;
		//
		c0->next=c->next;
		c->next=0;		
		delete c;c=0;
		if(c1==0) pdb->chn=c0;
		else 	  c1->next=c0;
		pdb->configure(); 	
		c=c0;	 		 
		if(TRES.logg>3) pdb->write("temp.pdb");
		fix(c);
		c1=c; 
	}	
}
 
void PdbFix::setfixedatm(Chn *chnorg,Bound *bound){

	Res *r;
	Atm *a;

	for(r=chnorg->res;r;r=r->next)  
	for(a=r->atm;a;a=a->next) {
		int i;
		for(i=0;i<bound->natom;i++){
			Atm *b;
			b=bound->atoms[i];
			if(b->res->id0!=a->res->id0) continue;
			if(b->tatm->id!=a->tatm->id) continue;
			//b->transfer(a);
			bound->predt[i]=0;
			break; 
		}
	}
}

void PdbFix::setfixedatmxyz(Chn *chnorg,Bound *bound){

	Res *r;
	Atm *a;

	for(r=chnorg->res;r;r=r->next)  
	for(a=r->atm;a;a=a->next) {
		int i;
		for(i=0;i<bound->natom;i++){
			Atm *b;
			b=bound->atoms[i];
			if(b->res->id0!=a->res->id0) continue;
			if(b->tatm->id!=a->tatm->id) continue;
			b->transfer(a);
			//bound->predt[i]=0;
			break; 
		}
	}
}



void PdbFix::fix(Chn *chn) {

	Bound bound;
	Pdb tmp;
	tmp.chn=chn;
	chn->pdb=&tmp;
	tmp.configure();
	bound.model=&tmp;
	bound.flag|=TRES.constant->pdbundeletable;
	bound.ready();
 
	//check real constraints
	Chn *chnorg=pdborg->ischain(chn->id);
	if(chnorg==0) {
		cerr<<"error in executing..."<<endl;
		cerr<<"chain with id: "<<chn->id<<" not found in pdb copy.."<<endl;
		exit(0);
	}
	bound.setatmat(chnorg); 
	//
	 	 

	//bound.setpredt(0);
        //bound.setpredt(begin,end->id0+1,1);
	//bound.setcloseresidue(); 
	
	bound.setpdbfixbound(2);
	setfixedatm(chnorg,&bound);

	bound.maxnit=500;
        bound.adjust();
        bound.reorder();
        bound.writeatmbound("out_0");
        bound.clearredundancy();
        bound.setimproper();
        bound.settag();
        //bound.boundsmooth(5);

	pdb->write("dist.pdb");
        bound.writeoutdistbound("dist.dat");
        bound.setbuffertag();
 
        bound.setbond();
	//bound.setsmallbound();
	bound.setrange();
	
	
        // 	 
        bound.writeatmbound("out_1");
        bound.boundsquare(2);
        bound.printoutimproper();
        //bound->setbuddy();
	setfixedatmxyz(chnorg,&bound);
        bound.searchstructure(10);
}

void PdbFix::fixold(Chn *chn) {

	Bound bound;
	Pdb tmp;
	tmp.chn=chn;
	chn->pdb=&tmp;
	tmp.configure();
	bound.model=&tmp;
	bound.flag|=TRES.constant->pdbundeletable;
	bound.ready();
	 
	//check real constraints
	Chn *chnorg=pdborg->ischain(chn->id);
	if(chnorg==0) {
		cerr<<"error in executing..."<<endl;
		cerr<<"chain with id: "<<chn->id<<" not found in pdb copy.."<<endl;
		exit(0);
	}
	setfixedatm(chnorg,&bound);
	//
	 
	//bound.setpredt(0);
        //bound.setpredt(begin,end->id0+1,1);
	bound.setcloseresidue(); 
	bound.maxnit=500;
        bound.adjust();
        bound.reorder();
        bound.clearredundancy();
        bound.setimproper();
        bound.settag();
        bound.boundsmooth(5);
        bound.writeoutdistbound("out0");
        bound.setbuffertag();
 
        bound.setbond();
	bound.setsmallbound();
	bound.setrange();
	
	
	bound.writeoutdistbound("out");
        // 	 
        bound.writeatmbound("out_1");
        bound.boundsquare(2);
        bound.printoutimproper();
        //bound->setbuddy();
        bound.searchstructure(10);
}
 
void PdbFix::fix() {
  
	Chn *c,*c1;
	 
	c1=0;
	for(c=pdb->chn;c;c=c->next) {   
		Chn *c0=createchain(c);
		c0->transform(3);
		c0->setallnear(c0->res,10000);
		c0->pdb=pdb;
		//
		c0->next=c->next;
		c->next=0;		
		delete c;c=0;
		if(c1==0) pdb->chn=c0;
		else 	  c1->next=c0;
		pdb->configure(); 	
		c=c0;	 		 
		if(TRES.logg>3) pdb->write("temp.pdb");
		fix(c);
		c1=c; 
		c->write("ss1");
	}	
}


void PdbFix::fixtemp() {
  
	Chn *c;
	
	for(c=pdb->chn;c;c=c->next) {   
		c->addmissingatms(c->res,100000);
		pdb->configure();
		if(TRES.logg>3) pdb->write("temp.pdb");
		fix(c);
	}	
}

Chn *PdbFix::createchain0(Chn *c) {

	Rotate rot;
	//create a new chain with the same sequecen as the old

	char *seqn=c->getseqn();
	Chn *c0=new Chn();
	c0->create(seqn);
	
	if(TRES.logg>3) c0->write("s");

	//configure them
	c0->configure();
	c->configure();
	c0->id=c->id;
	c0->start=c->start;
	
	
	//assign old id to new chains
	Res *r,*r0;	 
	for(r=c->res;r;r=r->next) {
		r0=c0->isres0(r->id0);
		if(r0==0||r0->name!=r->name) {
			cerr<<"program detected errors when trying to fix missing atoms..."<<endl;
			exit(0);
		}		 
		r0->oid=r->oid; //old id
		r0->id=r->id;   //new id
	}
	
	//check trans or cis
	for(r=c->res;r;r=r->next) {
		if(r->next==0) continue;
		if(r->next->id-r->id!=1) continue;
		Atm *ca,*ca0;
		ca=r->isatmid(1);
		ca0=r->next->isatmid(1);	
		if(ca==0|ca0==0) continue;
		//
		r0=c0->isres0(r->id0);
		if(r0==0) {
			cerr<<"program detected errors when trying to fix missing atoms..."<<endl;
			exit(0);
		}
		Atm *fa,*fa0;
		fa=r0->isatmid(1);
		fa0=r0->next->isatmid(1);	
		if(fa==0|fa0==0) continue;

		//cis case
		float d=TRES.distance(ca,ca0);	
		float rt;		
		if(d<3.0) {//cis 
			cerr<<"the CA distance between "<<r->name<<r->oid<<"--"<<r->next->name<<r->next->oid<<":"<<d<<endl;
			cerr<<"the distance is smaller than 3.0"<<endl;
			cerr<<"the peptide between the two residues considered as cis..."<<endl;	
			rt=r0->next->atm->gettorsionangle();					
			rot.rotateomega(r0->next->atm,-rt,100000,1);
		}	
		else {
			rt=r0->next->atm->gettorsionangle();
                        rot.rotateomega(r0->next->atm,180-rt,100000,1);  
		}
		
		if(d>4.5) {
			cerr<<"the CA distance between "<<r->name<<r->oid<<"--"<<r->next->name<<r->next->oid<<":"<<d<<endl;
			cerr<<"the distance is unusually large..."<<endl;
			cerr<<"the two residues should not be continous"<<endl;
			cerr<<"however their residue ids are sequential so they are treated as continuous anyway.."<<endl;
		}	
		//rotate omega to find better CA
		d=TRES.distance(ca,ca0);	
		float d0=TRES.distance(fa,fa0);
		
		if(TRES.logg>3) {
			rt=r0->next->atm->gettorsionangle();
			cerr<<"the distance new and old  "<<rt<<" "<<fa->res->id<<"--"<<fa0->res->id0<<":"<<d0<<" "<<d<<endl;
		}
		
		//if(d<3.5||d>4.0||exact==0) continue;
		if(exact==1)  {

			float x,y,z;
			float dx,dy,dz;
			dx=(fa0->xyz[0]-fa->xyz[0])/d0;
			dy=(fa0->xyz[1]-fa->xyz[1])/d0;
			dz=(fa0->xyz[2]-fa->xyz[2])/d0;
			x=(d-d0)/2*dx;
			y=(d-d0)/2*dy; 
			z=(d-d0)/2*dz; 
			Res *r1;
			Atm *a1;
			for(r1=r0->next;r1;r1=r1->next)
			for(a1=r1->atm;a1;a1=a1->next) {
				a1->xyz[0]+=x;
				a1->xyz[1]+=y;
				a1->xyz[2]+=z;
			}
			for(r1=c0->res;r1;r1=r1->next)
                	for(a1=r1->atm;a1;a1=a1->next) {
				if(r1->id0>r0->id0) break;
                        	a1->xyz[0]-=x;
                        	a1->xyz[1]-=y;
                        	a1->xyz[2]-=z;
                	}    
			if(TRES.logg>3) {
				d0=TRES.distance(fa,fa0);
                		rt=r0->next->atm->gettorsionangle();
           			cerr<<"**** the distance new and old  "<<rt<<" "<<fa->res->id<<"--"<<fa0->res->id0<<":"<<d0<<" "<<d<<endl;     
			}
		}
	}
	return c0;
}



Res **PdbFix::findbackbonelostsegment(Res *rr) {
	
	Res *r1,*r2;
	Res *r;	
	
	r1=rr->last;
	if(r1==0) r1=rr->chn->res;

	for(r=rr;r;r=r->next) {
		Res *r0=atmat[r->id0];	
		if(r0==0) continue;
		if(r0->isbackbonecomplete()==1) break; 
	}
			
	if(r==0) r=rr->chn->lastres();

	r2=r;
		
	Res **stem=new Res*[2];

	stem[0]=r1;

	stem[1]=r2;

	return stem;
}

void PdbFix::printhelp(){
fprintf(stderr,"*****************************************************\n");
fprintf(stderr,"ctrip is to build complete structure based on CA only\n");
fprintf(stderr,"comments? send to Jason Xiang at jsxzx@yahoo.com\n");		
fprintf(stderr,"*****************************************************\n");
fprintf(stderr,"\n");
fprintf(stderr,"USAGE: ctrip -prm num -k num -fast num file.pdb\n");
fprintf(stderr,"-prm select  atom model.default is 1 \n");
fprintf(stderr,"     1,all atom model\n");
fprintf(stderr,"     2,heavy atom model\n"); 
fprintf(stderr,"     3,backbone atom model\n");
fprintf(stderr,"-k   is original atoms kept beyond CA? default is 1.\n");
fprintf(stderr,"     0, delete all original atoms except CA\n");
fprintf(stderr,"     1, keep all original atoms\n");
//fprintf(stderr,"-fast  fast mode.from 0-5. default is 4. with 5 fastest\n");
}
void PdbFix::printhelpfix(){

fprintf(stderr,"*******************************************************************\n");
fprintf(stderr,"profix is to fix missing atoms and residues in pdb file\n");
fprintf(stderr,"the missing residues decided by comparing SEQRES card and ATOM card\n");
fprintf(stderr,"comments? send to Jason Xiang at jsxzx@yahoo.com\n");	
fprintf(stderr,"*******************************************************************\n");
fprintf(stderr,"\n");
fprintf(stderr,"USAGE: pdbfix -prm num -fix num -len num file.pdb\n");
fprintf(stderr,"-prm select  atom model.default is 1\n");
fprintf(stderr,"     1,all atom model\n");
fprintf(stderr,"     2,heavy atom model\n"); 
fprintf(stderr,"-len the maximum number of missing residues to be fixed. default is no limit\n");
fprintf(stderr,"-fix what to be fixed? default is 1.\n");
fprintf(stderr,"     0, fix missing atoms including missing atoms in backbone and sidechain\n");
fprintf(stderr,"     1, besides option 0, also fix missing residues by reading PDB SEQRES card\n");
//fprintf(stderr,"-fast  fast mode.from 0-5. default is 4. with 5 fastest\n");
}
