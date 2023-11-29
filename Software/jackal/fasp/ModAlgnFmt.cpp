#include "source.h"

ModAlgnFmt::ModAlgnFmt() {

	fileName=0;
	code=0;
	start=-1;
	end=-1;
	cid='0';
	seqn=0;
	seqngap=0;
    	match=0;
    	//score=0;
	compare=0;
	source=' ';
	pdb=0;
	next=0;
	more=0;
	up=0;
	last=0;
	token=0;
	flag=0;
	mutatepdb=0;
	sitescore=0;
}

ModAlgnFmt::~ModAlgnFmt() {

	Strhandler cc;
	fileName=cc.strdel(fileName);
	code=cc.strdel(code);
	seqn=cc.strdel(seqn);
	seqngap=cc.strdel(seqngap);
    	if(match) delete [] match;
    	//if(score) delete [] score;
	if(compare) delete [] compare;
	if(pdb&&(flag&TRES.constant->pdbundeletable)==0) delete pdb;
	if(next) delete next;
	if(more) delete more;
	if(mutatepdb) delete mutatepdb;
	if(sitescore) delete [] sitescore;
	last=0;up=0;
}

void ModAlgnFmt::setmutateresn() {


	ModAlgnFmt *t,*t0;

	for(t=this;t;t=t->next)
	for(t0=t;t0;t0=t0->more) {
		
		if(strcmp(t0->token,"structure")!=0) continue;

		int n=strlen(t0->seqn);

		for(int i=0;i<n;i++) {
			Res *r=resn[i];
			r=mutatepdb->chn->isres(r->id);
			resn[i]=r;
		}		
	}

}

void ModAlgnFmt::setresn() {

	if(seqn==0) return;

	int  n=strlen(seqn);

	resn=new Res *[n];

	int i;

	for(i=0;i<n;i++) resn[i]=0;

	if(start!=-1&&end!=-1) {

		start=start-pdb->chn->start;
		end=end-pdb->chn->start;
		Res *rr;
        	int nn=0;
        	for(i=0;i<n;i++) {
                    rr=pdb->chn->isres(i+start);
                    if(rr&&rr->name==seqn[i]) nn++;
        	}

        	if(nn>n/2) {
			for(i=0;i<n;i++) {
                        	rr=pdb->chn->isres(i+start);
                        	if(rr&&rr->name==seqn[i]) {
					resn[i]=rr;
				}
                	}
			return;
		}

	}
	else if(start!=-1) {
		start=start-pdb->chn->start;
		end=start+n-1;
	}
	else if(end!=-1) {
		end=end-pdb->chn->start;
		start=end-n+1;
	}
	else {
		start=0;
		end=pdb->chn->lastres()->id;
	}


	Algn algn;
	if(pdb->chn->seqcard==0) pdb->chn->setseqcard();
	algn.setsequence(pdb->chn->seqcard,seqn);
	algn.defalgn();

	int n1,n2;
	n1=0;
	n2=0;
	int slen=algn.getroutelength();
	char **result=algn.output(0);
	for(i=0;i<slen;i++) {
		if(result[0][i]!='-') n1++;
		if(result[1][i]!='-') n2++;
		if(result[0][i]==result[1][i]&&result[0][i]!='-') {
			resn[n2-1]=pdb->chn->isres0(n1-1);
		}
	}
}


void ModAlgnFmt::readModelAlgnFmt(char *f) {

	Strhandler cc;

	if(f==0) return;

	fileName=strdup(f);

	//char **lines = cc.opnfilebylinesimple(f);

	//if(lines==0) return;

	//int n= cc.gettotnum(lines);
	char line[1000],*lines;

	int i=0;
	//int k=0;

	//clear lines
	/*
	for(i=0;i<n;i++) {
		lines[i]=cc.clearfirstchar(lines[i]," \n\r");
	}

	k=0;
	for(i=0;i<n;i++) {
		if(lines[i]==0) continue;
		lines[k++]=lines[i];
	}

	lines[k]=0;

	n= cc.gettotnum(lines);
	*/

	ModAlgnFmt *tmp,*xmp;

	tmp=0;xmp=0;

	int num=0;

	int nl=strlen("#start new alignment");

	FILE *fp=fopen(f,"r");
	if(fp==0) {
		cerr<<"could not open file:"<<f<<endl;
		exit(0);
	}

	//for(i=0;i<n;i++) {
	while(fgets(line,1000,fp)!=NULL) {
		lines=strdup(line);
		if(lines==0) continue;
		//cc.clearendemptyspace(lines[i]);
		lines=cc.clearendchar(lines," \n\r");
		if(lines==0) continue;
		if(strncmp(lines,"#start new alignment",nl)==0) {
			if(xmp==0) {
				xmp=this;
			}
			else {
				xmp->next=new ModAlgnFmt();
				xmp->next->up=xmp;
				xmp=xmp->next;
			}
			tmp=0;
			continue;
		}
		if(lines[0]=='!') continue;
		if(strncasecmp(lines,">P1;",4)==0) {
			num=0;
			if(xmp==0) {
				tmp=this;
				xmp=this;
				num=0;
			}
			else if(tmp==0) {
				tmp=xmp;
				num=0;
			}
			else {
				tmp->more=new ModAlgnFmt();
				tmp->more->last=tmp;
				tmp=tmp->more;
				num=0;
			}
			num++;
		}
		else if(num==1) {
			num++;
			if(strncasecmp(lines,"structureX:",strlen("structureX:"))==0) {
				tmp->source='X';
				tmp->token=strdup("structure");
			}
			else if(strncasecmp(lines,"structureN:",strlen("structureN:"))==0) {
				tmp->source='N';
				tmp->token=strdup("structure");
			}
			else if(strncasecmp(lines,"structureM:",strlen("structureM:"))==0) {
				tmp->source='M';
				tmp->token=strdup("structure");
			}
			else if(strncasecmp(lines,"sequence:",strlen("sequence:"))==0) {
                                tmp->source='U';
                                tmp->token=strdup("sequence");
                        }
			else if(strncasecmp(lines,"squence:",strlen("sequence:"))==0) {
                                tmp->source='U';
                                tmp->token=strdup("sequence");
                        }
			else if(strncasecmp(lines,"target:",strlen("target:"))==0) {
                                tmp->source='U';
                                tmp->token=strdup("target");
                        }

			char **tt=cc.pairbytokensimple(lines,":");

			int to=cc.gettotnum(tt);

			for(int j=0;j<to;j++) {
				if(j==1) {
					tmp->code=strdup(tt[j]);
					tmp->code=cc.clearendchar(code,"\r\t\n ");
				}
				else if(j==2) {
					tmp->start=atoi(tt[j]);
				}
				else if(j==3) {
					if(tt[j]&&strlen(tt[j])>0) {
						tmp->cid=tt[j][0];
					}
				}
				else if(j==4) {
					tmp->end=atoi(tt[j]);
				}
				else if(j==5) {
					if(tt[i]&&strlen(tt[j])>0) {
                                                tmp->cid=tt[j][0];
                                        }
				}
			}
			tt=cc.strdel(tt);
		}
		else {
			tmp->seqngap=cc.straddup(tmp->seqngap,lines);
		}
		if(lines) delete [] lines;lines=0;
	}

	ModAlgnFmt *c,*c0;

        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
	    c0->seqngap=cc.clearchar(c0->seqngap," *\n\t\r");
	    if(c0->seqngap&&c0->seqngap[0]=='\0') {
		delete [] c0->seqngap; c0->seqngap=0;
	    }
	    if(c0->seqngap) cerr<<c0->seqngap<<endl;
        }

	removeallnonstandard();
	clearemptyseq();
	setalldefaultseq();
	printoutallnoseqmesg();
}

void ModAlgnFmt::printoutallnoseqmesg()
{
	ModAlgnFmt *c;
        for(c=this;c;c=c->next) c->printoutnoseqmesg();
}

void ModAlgnFmt::printoutnoseqmesg()
{
        ModAlgnFmt *c;

        for(c=this;c;c=c->more) {
		if(strcmp(c->token,"sequence")==0) return;
	}

	cerr<<"Warning! some alignments does not have query"<<endl;
	cerr<<"guess you have noticed it"<<endl;
}


void ModAlgnFmt::clearemptyseq() {

        ModAlgnFmt *c,*c0;

        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
		if(c0->seqngap==0) continue;
		int n=strlen(c0->seqngap);
		int m=0;
		for(int i=0;i<n;i++) {
			if(c0->seqngap[i]!='-') {
				m=1;
				break;
			}
		}
		if(m==0) {
			delete [] c0->seqngap;
			c0->seqngap=0;
		}
        }
}


void ModAlgnFmt::deletegap() {

	if(seqn) return;

	if(seqngap==0) return ;

	int n=strlen(seqngap);

	seqn=new char[n+1];

	match=new int[n+1];
	compare=new int[n+1];

	int i;

	for(i=0;i<n+1;i++) match[i]=-1;

	for(i=0;i<n+1;i++) compare[i]=-1;

	int m=0;

	for(i=0;i<n;i++) {
		if(seqngap[i]=='-') continue;
		seqn[m]=seqngap[i];
		match[i]=m;
		compare[m]=i;
		m++;
	}

	seqn[m]='\0';

	return;
}

void ModAlgnFmt::deleteallgap() {

	ModAlgnFmt *t;
	ModAlgnFmt *tt;

	for(tt=this;tt;tt=tt->next)
	for(t=tt;t;t=t->more) {
		t->deletegap();
	}
}

void ModAlgnFmt::takeofzero() {

        ModAlgnFmt *t;
        ModAlgnFmt *tt;

        for(tt=this;tt;tt=tt->next)
        for(t=tt;t;t=t->more) {
		if(t->token==0) token=strdup("UNKNOWN");
        }
}

ModAlgnFmt *ModAlgnFmt::findfirstsequencefmt() {

        ModAlgnFmt *t;
        ModAlgnFmt *tt;

        for(tt=this;tt;tt=tt->next)
        for(t=tt;t;t=t->more) {
                if(strcmp(t->token,"sequence")==0) return t;
        }

	return 0;
}

void ModAlgnFmt::removeallnonstandard() {

	ModAlgnFmt *c,*c0;

        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
        	c0->removenonstandard();
        }
}

void ModAlgnFmt::removenonstandard() {

	//if(s==0) return;

	if(seqngap==0) return;

	int n=strlen(seqngap);

	Tres *tt;

	for(int i=0;i<n;i++) {
		if(seqngap[i]=='-') continue;
		tt=TRES[seqngap[i]];
		if(tt==0)   {
			cerr<<"the residue :"<<i<<" "<<seqngap[i]<<" is not standard! treated as gap!"<<endl;
			seqngap[i]='-';
		}
	}
	return;
}

void ModAlgnFmt::createstructure(char *s) {

        pdb=new Pdb();
        pdb->chn=new Chn();
        pdb->chn->create(s);
        pdb->configure();
	pdb->chn->header();
}

void ModAlgnFmt::createstructure() {

	createstructure(seqn);
}

void ModAlgnFmt::setseqnstructure() {

	ModAlgnFmt *c,*c0,*c1;

	int n;
	int nn=0;
	c1=0;
        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
               if(strcmp(c0->token,"sequence")!=0) continue;
	       if(seqn==0) continue;
	       n=strlen(seqn);
	       if(nn<n) {
		     nn=n;
	             c1=c0;
	       }
        }

	if(c1==0||c1->seqn==0||strlen(c1->seqn)==0) {

		cerr<<"the sequence of the protein to be built does not exist!"<<endl;
		exit(0);
	}

	c1->createstructure(c1->seqn);
	c1->pdb->chn->start=c1->start;
	c1->setresn();

	for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
		if(strcmp(c0->token,"sequence")!=0) continue;
		if(c0==c1) continue;
		//c0->createstructure(c1->seqn);
		//c0->pdb->chn->start=c1->start;
		c0->pdb=c1->pdb;
		c0->setresn();
		c0->flag|=TRES.constant->pdbundeletable;
        }
}

ModAlgnFmt *ModAlgnFmt::getlongestsequence() {

	ModAlgnFmt *c,*c0,*c1;

        int n;
        int nn=0;
        c1=0;

        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
               if(strcmp(c0->token,"sequence")!=0) continue;
               if(seqn==0) continue;
               n=strlen(seqn);
               if(nn<n) {
                     nn=n;
                     c1=c0;
               }
        }

	return c1;
}


void ModAlgnFmt::setseqnstructure(Pdb *p) {

	ModAlgnFmt *c,*c0;

	for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
                if(strcmp(c0->token,"sequence")!=0) continue;
                //c0->createstructure(c1->seqn);
                //c0->pdb->chn->start=c1->start;
		c0->pdb=p;
		c0->flag |= TRES.constant->pdbundeletable;
                c0->setresn();
        }
}

void ModAlgnFmt::setalldefaultseq() {

	ModAlgnFmt *c,*c0;

        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
                if(strcmp(c0->token,"structure")==0) c0->setdefaultstructseq();
		else if(strcmp(c0->token,"sequence")==0) c0->setdefaultqueryseq();
        }
}

void ModAlgnFmt::setdefaultqueryseq() {

	ModAlgnFmt *c=getrootModAlgnFmt();

	int n=-1000;

        ModAlgnFmt *t,*t0,*g=0;

        for(t=c;t;t=t->next)
	for(t0=t;t0;t0=t0->more) {
		if(strcmp(t0->token,"sequence"))continue;
		if(t0->seqngap==0) continue;
		int m=strlen(t0->seqngap);
                if(m>n) {
                        n=m;
                        g=t0;
                }
        }

        if(g&&g->seqngap) {
                seqngap=strdup(g->seqngap);
        }
}

void ModAlgnFmt::setdefaultstructseq() {

	if(seqngap) return;

	if(token==0) {
		cerr<<"the code does not exist for the alignment:"<<endl;
		cerr<<seqngap<<endl;
		return;
	}

	if(strcmp(token,"structure")!=0) return;

	Pdb ppdb;
	ppdb.read(code,cid);
	if(ppdb.chn==0) {
		return;
	}

	int n=-10000;
	Chn *t=0;
	for(Chn* c=ppdb.chn;c;c=c->next) {
		char *s=c->getseqn();
		if(s==0) continue;
		int m=strlen(s);
		if(m>n) {
			n=m;
			t=c;
		}
		delete [] s;s=0;
	}

	if(t) {
		seqngap=t->getseqn();
	}

}
void ModAlgnFmt::getstructure() {

	pdb=new Pdb();
	pdb->read(code,cid);
	if(pdb->chn==0) {
		cerr<<"the protein structure "<<code<<" "<<cid<<" does not exist"<<endl;
		exit(0);
	}
	deletechain();
	setresn();
	/*
	for(Chn *chn=pdb->chn;chn;chn=chn->next) {
		chn->header();
		chn->buildhbond();
		chn->setdsspstr();
		chn->write("sss");
		for(Res *r=chn->res;r;r=r->next) cerr<<r->id<<r->name<<" "<<r->sec<<endl;
		for(r=chn->res;r;r=r->next) {
			for(HBondList *h=r->hbond;h;h=h->next) {
				cerr<<h->donor->res->name<<" "<<h->donor->res->id<<" "<<h->donor->name<<" ";
				cerr<<h->acceptor->res->name<<" "<<h->acceptor->res->id<<" "<<h->acceptor->name<<endl;
			}
		}
		
	}*/
}


void ModAlgnFmt::getallstructure() {

	ModAlgnFmt *c,*c0;

	for(c=this;c;c=c->next)
	for(c0=c;c0;c0=c0->more) {

		if(strcmp(c0->token,"structure")==0) {
			c0->getstructure();
		}
	}
}


void ModAlgnFmt::deleteallchain() {

	ModAlgnFmt *c,*c0;

        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {

                if(strcmp(c0->token,"structure")==0) {
                        c0->deletechain();
                }
        }
}
void ModAlgnFmt::deletechain() {

	if(pdb==0||pdb->chn==0||pdb->chn->next==0) return;

	Chn *c;
	Chn *c0;
	float dmax=0;

	//delete chains have the same id iteratively

	while(1) {

		int n=0;
		char cn=' ';

		for(c=pdb->chn;c;c=c->next) {
			cn=c->id;
			for(c0=c->next;c0;c0=c0->next) {
				if(c0==c) continue;
				if(c0->id==cn) {
					n=1;
					pdb->removechain(c0);
					break;
				}
			}
			if(n) break;
		}
		if(n==0) break;
	}

	//delete chains have different id

	if(pdb==0||pdb->chn==0||pdb->chn->next==0) return;

	for(c=pdb->chn;c;c=c->next) {
		c->setseqcard();
		Algn aa;
		aa.setsequence(c->seqcard,seqn);
		aa.defalgn();
		float d=aa.calctotalscore();
		if(c==pdb->chn) {
			dmax=d;
			c0=c;
		}
		else if(d>dmax) {
			dmax=d;
			c0=c;
		}
	}

	while(1) {
		int n=0;
		for(c=pdb->chn;c;c=c->next) {
			if(c!=c0) {
				pdb->removechain(c);
				n=1;
				break;
			}
		}
		if(n==0) break;
	}
	pdb->configure();
	if(pdb->chn) cid=pdb->chn->id;
}

ModAlgnFmt * ModAlgnFmt ::findsequencefmt() {

	ModAlgnFmt *c;

	for(c=this;c;c=c->more) {

		if(c->token&&strcmp(c->token,"sequence")==0) return c;

	}

	return 0;
}
void ModAlgnFmt::setmutatesec() {
	
	ModAlgnFmt *c;

        for(c=this;c;c=c->next) {
                if(strcmp(c->token,"structure")) continue;
                Chn *cc;
		for(cc=c->mutatepdb->chn;cc;cc=cc->next) {
			cc->header();
			cc->buildhbond();
			cc->setdsspstr();
			
			cc->write("sss");
			Res *r;
                	for(r=cc->res;r;r=r->next) cerr<<r->id<<r->name<<" "<<r->sec<<endl;
                	for(r=cc->res;r;r=r->next) {
                        	for(HBondList *h=r->hbond;h;h=h->next) {
                                	cerr<<h->donor->res->name<<" "<<h->donor->res->id<<" "<<h->donor->name<<" ";
                                	cerr<<h->acceptor->res->name<<" "<<h->acceptor->res->id<<" "<<h->acceptor->name<<endl;
				}
                        }
		}
        }
}

void ModAlgnFmt::setmutatesidepdb() {

	ModAlgnFmt *c;

	for(c=this;c;c=c->next) {
		if(strcmp(c->token,"structure")) continue;
		c->setmutatesidechainpdb();
	}
}

ModAlgnFmt *ModAlgnFmt::getparentModAlgnFmt() {

	ModAlgnFmt *c;

        for(c=this;c->last;c=c->last) ;

	return c;
}

ModAlgnFmt *ModAlgnFmt::getrootModAlgnFmt() {

        ModAlgnFmt *c,*c0;
        for(c=this;c->last;c=c->last);
	for(c0=c;c0->up;c0=c0->up);
        return c0;
}


ModAlgnFmt *ModAlgnFmt::getsequenceModAlgnFmt() {

        ModAlgnFmt *c,*a;

	c=getparentModAlgnFmt();

	for(a=c;a;a=a->more) {
		if(strcmp(a->token,"sequence")==0) return a;
	}
        return 0;
}

void ModAlgnFmt::checklength() {

	ModAlgnFmt *c;

        for(c=this;c;c=c->next) {
                c->checkalignmentlength();
        }
}



void ModAlgnFmt::checkalignmentlength() {

	ModAlgnFmt *c;

	int n=-1000;
        for(c=this;c;c=c->more) {
               n=max(n,strlen(c->seqngap));
        }

	for(c=this;c;c=c->more) {
		int i=strlen(c->seqngap);
		c->seqngap=(char *)realloc(c->seqngap,sizeof(char)*n);
		for(int j=i;j<n;j++) seqngap[j]='-';
	}
}


void ModAlgnFmt::setmutatesidechainpdb() {

	if(strcmp(token,"structure")) return;
	if(mutatepdb==0) {
		mutatepdb=new Pdb(pdb);
		mutatepdb->configure();
	}
	mutatepdb->setflgr(-99999);

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
	scprd.arbt=0;
	scprd.resultout=0;
	scprd.pdb=mutatepdb;

	ModAlgnFmt *parent=getparentModAlgnFmt();
	ModAlgnFmt *se = parent->getsequenceModAlgnFmt();
	if(se==0) return;

	int nlen=strlen(se->seqngap);

	for(int i=0;i<nlen;i++) {
		if(match[i]==-1||se->match[i]==-1) continue;
		//if(match[i]==se->match[i]) continue;
		int ia=match[i];
		int ib=se->match[i];
		if(seqn[ia]==se->seqn[ib]) continue;
		Res *rr=resn[ia];
		if(rr==0) continue;
		rr=mutatepdb->isres(rr->id0);
		if(rr==0) continue;

		if(TRES[se->seqn[ib]]==0) continue;
		rr->mutateResidue(se->seqn[ib]);
		if(rr->name=='G') rr->flag=-99999;
		else rr->flag=10000;
	}
	mutatepdb->configure();
	mutatepdb->chn->start=pdb->chn->start;
	mutatepdb->write("mutated.out1");
	int smt0=TRES.smt;
	TRES.smt=TRES.nsmt-1;
	//mutatepdb->transform(5,-99999);
	//scprd.hookside();
	pdb->setflgr('A',-99999);
	pdb->setflgr('P',-99999);
	/*
	if(TRES.isatm(" HA ")&&TRES.isatm(" HN ")) {
 		scprd.addbackboneh();
	}
	*/
	scprd.scpred(" "," ",1);
	TRES.smt=smt0;
	scprd.pdb->write("mutate.out");
	pdb->write("mutate.out2");
}

void ModAlgnFmt::prepare() {
	deleteallgap();
	getallstructure();
	setseqnstructure();
	setmutatesidepdb();
	setmutateresn();
	setmutatesec();
	setsitescore();
}


float *ModAlgnFmt::gethbondbounds(Atm *a,Atm *a0){

	Res *r=a->res;
	Res *r0=a0->res;

    	int n=compare[r->id0];
    	int n0=compare[r0->id0];

    	if(seqngap[n]!=r->name||seqngap[n0]!=r0->name) {
        	cerr<<"Warning! the residue id and the name does not match";
    	}

    	ModAlgnFmt *t,*t0;

    	//float *dst=0;

	float tt[1000];
	ModAlgnFmt *mtt[1000];

	//get all constraints

	int ntt=0;
    	for(t=this;t;t=t->next) {
        	for(t0=t;t0;t0=t0->more) {
           		if(strcmp(t0->token,"sequence")==0) continue;
			if(strcmp(t0->token,"structure")!=0) continue;

			int i=t->match[n];
			int i0=t0->match[n0];
			Res *rr=resn[i];
			Res *rr0=resn[i0];
			if(rr==0||rr0==0) continue;

			if(rr->tres!=r->tres||rr0->tres!=r0->tres) continue;

			Atm *aa,*aa0;

			aa=rr->isatmid(a->tatm->id);
			aa0=rr0->isatmid(a0->tatm->id);

			if(aa==0||aa0==0) continue;

			float dis=TRES.distance(aa,aa0);

			tt[ntt]=dis;
			mtt[ntt]=t0;
			ntt++;
       		}
    	}

	//find the most conserved constraint

	int i=0;

	DistPopular *dst=TRES.popbin->getfaraway();
	DistPopular *p=dst->getDistPopular(a->res->tres,a->name,a0->name,0);

	float g1,g2,g3;
	float s1;
	g1=g2=g3=0;
	s1=0;
	for(i=0;i<ntt;i++) {
		t=mtt[i];
		float d=t->getalgnweight(i);
		float a1=tt[i]-p->dist[0];
		float a2=p->dist[1]-tt[i];
		g1+=tt[i]*d;
		s1+=d;
		g2+=a1*d;
		g3+=a2*d;
	}
	g1=g1/s1;
	g2=g2/s1;
	g3=g3/s1;

	g2=max(p->dist[0],g1-g2);
	g3=min(p->dist[1],g1+g3);

	float *dp=new float[8];
	for(i=0;i<8;i++) dp[i]=0;
	dp[0]=g1;
	dp[1]=g2;
	dp[2]=g3;
	dp[3]=10;

	return dp;
}

float *ModAlgnFmt::getfardistbounds(Atm *a,Atm *a0) {

	//float dp[1000];

	//int ndp=0;

	Res *r=a->res;
    	Res *r0=a0->res;

    	int n=compare[r->id0];
    	int n0=compare[r0->id0];

    	if(seqngap[n]!=r->name||seqngap[n0]!=r0->name) {
        	cerr<<"Warning! the residue id and the name does not match";
    	}

    	ModAlgnFmt *t,*t0;

    	//float *dst=0;

	float tt[1000];
	ModAlgnFmt *mtt[1000];

	//get all constraints

	int ntt=0;
    	for(t=this;t;t=t->next) {
        	for(t0=t;t0;t0=t0->more) {
            		if(strcmp(t0->token,"sequence")==0) continue;
			if(strcmp(t0->token,"structure")!=0) continue;

			int i=t->match[n];
			int i0=t0->match[n0];
			Res *rr=resn[i];
			Res *rr0=resn[i0];
			if(rr==0||rr0==0) continue;

			if(rr->tres!=r->tres||rr0->tres!=r0->tres) continue;

			Atm *aa,*aa0;

			aa=rr->isatmid(a->tatm->id);
			aa0=rr0->isatmid(a0->tatm->id);

			if(aa==0||aa0==0) continue;

			float dis=TRES.distance(aa,aa0);

			tt[ntt]=dis;
			mtt[ntt]=t0;
			ntt++;
       		}
    	}

	//find the most conserved constraint

	int i=0;

	DistPopular *dst=TRES.popbin->getnearnext();
	DistPopular *p=dst->getDistPopular(a->res->tres,a->name,a0->name,a0->res->id-a->res->id);

	float g1,g2,g3;
	float s1;
	g1=g2=g3=0;
	s1=0;
	for(i=0;i<ntt;i++) {
		t=mtt[i];
		float d=t->getalgnweight(i);
		float a1=tt[i]-p->dist[0];
		float a2=p->dist[1]-tt[i];
		g1+=tt[i]*d;
		s1+=d;
		g2+=a1*d;
		g3+=a2*d;
	}
	g1=g1/s1;
	g2=g2/s1;
	g3=g3/s1;

	g2=max(p->dist[0],g2);
	g3=min(p->dist[1],g3);

	float *dp=new float[8];
	for(i=0;i<8;i++) dp[i]=0;
	dp[0]=g1;
	dp[1]=g2;
	dp[2]=g3;
	dp[3]=10;

	return dp;
}


float *ModAlgnFmt::getdistbounds(Atm *a,Atm *a0) {


	Res *r=a->res;
    	Res *r0=a0->res;

    	ModAlgnFmt *t,*t0;


	float tt[1000];
	ModAlgnFmt *mtt[1000];

	//get all constraints

	int ntt=0;
    	for(t=this;t;t=t->next) {
		ModAlgnFmt *s=t->findsequencefmt();
		int n=s->compare[r->id0];
		int n0=s->compare[r0->id0];
        	for(t0=t;t0;t0=t0->more) {
            		if(strcmp(t0->token,"sequence")==0) continue;
			if(strcmp(t0->token,"structure")!=0) continue;

			int i=t0->match[n];
			int i0=t0->match[n0];
			Res *rr=t0->resn[i];
			Res *rr0=t0->resn[i0];
			if(rr==0||rr0==0) continue;

			if(rr->tres!=r->tres||rr0->tres!=r0->tres) {
				cerr<<"should be equal:"<<n<<" "<<n0<<endl;
				continue;
			}
			Atm *aa,*aa0;

			aa=rr->isatmid(a->tatm->id);
			aa0=rr0->isatmid(a0->tatm->id);

			if(aa==0||aa0==0) continue;

			float dis=TRES.distance(aa,aa0);

			tt[ntt]=dis;
			mtt[ntt]=t0;
			ntt++;
       		}
    	}

	//find the most conserved constraint

	int i=0;

	DistPopular *dst=TRES.popbin->getnearnext();
	DistPopular *p=dst->getDistPopular(a->res->tres,a->name,a0->name,a0->res->id-a->res->id);

	float g1,g2,g3;
	float s1;
	g1=g2=g3=0;
	s1=0;
	for(i=0;i<ntt;i++) {
		t=mtt[i];
		float c=t->getalgnscore();
		float d=t->getalgnweight(i);
		float a1=tt[i]-p->dist[0];
		float a2=p->dist[1]-tt[i];
		g1+=tt[i]*d;
		s1+=d;
		g2+=a1*d;
		g3+=a2*d;
	}
	g1=g1/s1;
	g2=g2/s1;
	g3=g3/s1;

	g2=max(p->dist[0],g1-g2);
	g3=min(p->dist[1],g1+g3);

	float *dp=new float[8];
	for(i=0;i<8;i++) dp[i]=0;
	dp[0]=g1;
	dp[1]=g2;
	dp[2]=g3;
	dp[3]=10;

	return dp;
}

void ModAlgnFmt::setsitescore() {

	ModAlgnFmt *t,*t0;
        for(t0=this;t0;t0=t0->next) {

                for(t=t0;t;t=t->more) {
			if(strcmp(t->token,"structure")) continue;
			
			if(t->sitescore) delete [] t->sitescore;
			int n=strlen(t->seqngap);
			t->sitescore=new float[n+1];
			int i=0;
			for(i=0;i<n+1;i++) t->sitescore[i]=0;		
					
			for(i=0;i<n;i++) t->sitescore[i]=t->getalgnweight(i);
		}
	}
}
float *ModAlgnFmt::getdistbounds(Atm *a,Atm *a0,float low,float high) {

	 
	Res *r=a->res;
    	Res *r0=a0->res;
 
    	ModAlgnFmt *t,*t0;
 
	float tt[200];
	ModAlgnFmt *mtt[200];

	//get all constraints

	int ntt=0;
	int men=0;
    	for(t=this;t;t=t->next) {
		ModAlgnFmt *s=t->findsequencefmt();
		int n=s->compare[r->id0];
		int n0=s->compare[r0->id0];
		 
        	for(t0=t;t0;t0=t0->more) {
            		if(strcmp(t0->token,"sequence")==0) continue;
			if(strcmp(t0->token,"structure")!=0) continue;

			int i=t0->match[n];
			int i0=t0->match[n0];
			Res *rr=t0->resn[i];
			Res *rr0=t0->resn[i0];
			if(rr==0||rr0==0) continue;

			if(rr->tres!=r->tres||rr0->tres!=r0->tres) {
				cerr<<"should be equal:"<<n<<" "<<n0<<endl;
				continue;
			}

			Atm *aa,*aa0;

			aa=rr->isatmid(a->tatm->id);
			aa0=rr0->isatmid(a0->tatm->id);

			if(aa==0||aa0==0) continue;

			float dis=TRES.distance(aa,aa0);

			if(men++>190) break;
			tt[ntt]=dis;
			mtt[ntt]=t0;
			ntt++;
       		}
    	}

	//find the most conserved constraint

	int i;

	for(i=0;i<ntt;i++) {
		low=min(low,tt[i]-0.1);
		high=max(high,tt[i]+0.1);
	}
	
	float g1,s1;

	s1=g1=0;

	for(i=0;i<ntt;i++) {
		t=mtt[i];
		float d=t->getalgnweight(i);
		g1+=tt[i]*d;
		s1+=d;
	}
	if(s1==0) s1=1;
	if(g1==0) g1=0;
	g1=g1/s1;

	float *dp=new float[8];
	for(i=0;i<8;i++) dp[i]=0;
	dp[0]=g1;
	dp[1]=low;
	dp[2]=high;
	if(high-low<2)  dp[3]=1000;
	else 		dp[3]=100;

	return dp;
}

float ModAlgnFmt::getalgnweight(int n)  {

	int len=strlen(seqngap);

        ModAlgnFmt *t=getparentModAlgnFmt();
        ModAlgnFmt *s0=t->findsequencefmt();

	Algn algn;

	//return 0.5;
	float *blosum65mt=algn.getscoretable("blosum65mt");

        float a=0;
        float am=0;
        float bm=0;
        for(int i=0;i<len;i++)  {
                if(seqngap[i]=='-'&&s0->seqngap[i]=='-') continue;
                int n1= TRES.getresid(seqngap[i]);
                int n2= TRES.getresid(s0->seqngap[i]);
                if(n1==-1&&n2==-1) continue;
		float x=1+abs(i-n);
                if(n1==-1) {
                        bm+=blosum65mt[n2*(n2+1)/2+n2]/x;
                }
                else if(n2==-1) {
                        am+=blosum65mt[n1*(n1+1)/2+n1]/x;
                }
                else {
                        am+=blosum65mt[n1*(n1+1)/2+n1]/x;
                        bm+=blosum65mt[n2*(n2+1)/2+n2]/x;
                        if(n1<n2) a+=blosum65mt[n1*(n1+1)/2+n2]/x;
                        else      a+=blosum65mt[n2*(n2+1)/2+n1]/x;
                }
        }

	if(blosum65mt) delete [] blosum65mt;
        if(am*bm==0) return 0;

        float c=a/sqrt(am*bm);

        return c;
}

float ModAlgnFmt::getalgnscore() {

	int len=strlen(seqngap);

	ModAlgnFmt *t=getparentModAlgnFmt();
	ModAlgnFmt *s0=t->findsequencefmt();

	Algn algn;

	float *blosum65mt=algn.getscoretable("blosum65mt");

	float a=0;
	float am=0;
	float bm=0;
	for(int i=0;i<len;i++) 	{
		if(seqngap[i]=='-'&&s0->seqngap[i]=='-') continue;
		int n1=	TRES.getresid(seqngap[i]);
		int n2= TRES.getresid(s0->seqngap[i]);
		if(n1==-1&&n2==-1) continue;
		if(n1==-1) {
			bm+=blosum65mt[n2*(n2+1)/2+n2];
		}
		else if(n2==-1) {
			am+=blosum65mt[n1*(n1+1)/2+n1];
		}
		else {
			am+=blosum65mt[n1*(n1+1)/2+n1];
			bm+=blosum65mt[n2*(n2+1)/2+n2];
			if(n1<n2) a+=blosum65mt[n1*(n1+1)/2+n2];
			else 	  a+=blosum65mt[n2*(n2+1)/2+n1];
		}
	}

	if(blosum65mt) delete [] blosum65mt;

	if(am*bm==0) return 0;

	float c=a/sqrt(am*bm);

	return c;
}


void ModAlgnFmt::setconservedres() {

	//prepare score
	float *score=0;

	if(score) delete [] score; score=0;
	int n=strlen(seqngap);
	score=new float[n];

	int i=0;
	for(i=0;i<n;i++) { score[i]=0; }

	//calculate score of each alignment

	ModAlgnFmt *m,*s0;


	s0=findsequencefmt();

	if(s0==0) return;

	Algn algn;
	float *blosum65mt=algn.getscoretable("blosum65mt");

	for(i=0;i<n;i++){

		int c0=TRES.getresid(s0->seqngap[i]);
		if(c0==-1) continue;
		for(m=this;m;m=m->more) {
			if(m==s0) continue;
			if(m->seqngap[i]=='-'||s0->seqngap[i]=='-') continue;

			int ci=TRES.getresid(m->seqngap[i]);
			if(ci==-1) continue;

			if(ci<c0) score[i]=blosum65mt[ci*(ci+1)/2+c0];
			else      score[i]=blosum65mt[c0*(c0+1)/2+ci];
		}
	}
	if(blosum65mt)delete [] blosum65mt;
}

void ModAlgnFmt::setallconservedres() {

	ModAlgnFmt *t;
	for(t=this;t;t=t->next) {
		t->setconservedres();
	}
}

int ModAlgnFmt::isconservedpair(int n,int n0) {

	/*
	n=compare[n];
	n0=compare[n0];

	ModAlgnFmt *t,*t0;
	for(t=this;t;t=t->next) {
		for(t0=t;t0;t0=t0->more) {
			if(strcmp(t0->token,"structure")!=0) continue;
			if(t0->seqngap[n]=='-'||t0->seqngap[n0]=='-') continue;
			int i=match[n];
			int i0=match[n0];
			if(resn[i]->sec!='-'&&resn[i0]->sec!='-') {

			}
			//t0->pdb->chn->buildhbond();
			//t0->pdb->chn->setdsspstr();
		}
	}
	*/
	return 1;

}

int ModAlgnFmt::getcreditofhbond(Atm *a,Atm *a0) {
	
	int ii=0;

	if(a->tatm->id==0&&a0->tatm->id==3) ii=1;
	else if(a->tatm->id==3&&a0->tatm->id==0) ii=1;
	else if(strcmp(a->tatm->name," SG ")==0&&strcmp(a0->tatm->name," SG ")==0) ii=2;
	else ii=3;

	if(ii==3) return 1; //general salt bridges

	if(ii==2) return 10; //ssbond

	Res *r=a->res;
	Res *r0=a0->res;
	
	if(r->sec=='h'&&r0->sec=='h'&&fabs(r->id-r0->id)<=5&&fabs(r->id-r0->id)>=3) {
		return 100; //secondary bond
	}
	else if(r->sec=='e'&&r0->sec=='e') {
		return 50;
	}
	else return 5;
}

float *ModAlgnFmt::gethbondbounds(Res *r,Res *r0) {

	//take all hydrogens out including s-s bond and ionic bond

        ModAlgnFmt *mtt[1000];
	Atm	   *aaa[1000];
	float       aln[500];
	int         credit[500];

	ModAlgnFmt *t,*t0;

	int ntt=0;

        for(t=this;t;t=t->next) {
                ModAlgnFmt *s=t->findsequencefmt();
                int n=s->compare[r->id0];
                int n0=s->compare[r0->id0];
                for(t0=t;t0;t0=t0->more) {

			if(ntt>300) break;
                        if(strcmp(t0->token,"structure")!=0) continue;

                        if(t0->seqngap[n]=='-'||t0->seqngap[n0]=='-') continue;

                        int i=t0->match[n];

                        int i0=t0->match[n0];

                        Res *rr=t0->resn[i];
                        Res *rr0=t0->resn[i0];
                        if(rr==0||rr0==0) continue;

                        if(rr->ishbondedwithres(rr0)==0&&rr0->ishbondedwithres(rr)==0) continue;
			
			HBondList *h;

			for(h=rr->hbond;h;h=h->next) {
				if(h->donor->res==rr&&h->acceptor->res!=rr0) continue;
				if(h->donor->res==rr0&&h->acceptor->res!=rr) continue;
				
				if(h->donor->res==rr) {
					aaa[ntt*2]=h->donor;
					aaa[ntt*2+1]=h->acceptor;
				}
				else {
					aaa[ntt*2]=h->acceptor;
                                        aaa[ntt*2+1]=h->donor;
				}
				float d=t0->getalgnweight(n)*t0->getalgnweight(n0);
				d=sqrt(d);
				mtt[ntt]=t0;
				aln[ntt]=d;
				credit[ntt]=getcreditofhbond(aaa[ntt*2],aaa[ntt*2+1]);
				ntt++;
			}			
                }
        }

	if(ntt==0) return 0;

	int i,j;
	
	//calculate the distance constraint

	DistPopular *dst=TRES.popbin->gethbond();
	float *dp=new float[ntt*6+10];
	for(i=0;i<ntt*6+10;i++) dp[i]=0;
	int htt=0;
	//add
	//return dp;
	//add
	for(i=0;i<ntt;i++) {
		if(mtt[i]==0) continue;
		float ad=0,bd=0,cd=0,dd=0;
		float de=0;
		float dt=0;
		float dx=0;
		for(j=0;j<ntt;j++) {

			if(mtt[j]==0) continue;
			if(aaa[2*i]!=aaa[2*j]) continue;
			if(aaa[2*i+1]!=aaa[2*j+1]) continue;
			
			t=mtt[j];
			float d=aln[j];
			float x=TRES.distance(aaa[2*i],aaa[2*i+1]);
			ad+=d*x;
			x=TRES.distance(aaa[2*i]->bond[0],aaa[2*i+1]);
			bd+=d*x;
			x=TRES.distance(aaa[2*i],aaa[2*i+1]->bond[0]);
			cd+=d*x;
			x=TRES.distance(aaa[2*i]->bond[0],aaa[2*i+1]->bond[0]);
                        dd+=d*x;
			if(credit[j]==100) de+=d*2*(1-d);
			else if(credit[j]==50) de+=d*3*(1-d);
			else if(credit[j]==5) de+=d*5*(1-d);
			else if(credit[j]==10) de+=d*4*(1-d);
			else if(credit[j]==1) de+=d*7*(1-d);
			dt+=d;
			dx=max(dx,d);
		}
		for(j=0;j<ntt;j++) {
			if(aaa[2*i]!=aaa[2*j]) continue;
			if(aaa[2*i+1]!=aaa[2*j+1]) continue;
			mtt[i]=0;	
		}
		
		if(dt>0) {
			Atm *aa0=0,*bb0=0;

			if(aaa[2*i]->tatm->id==0) aa0=aaa[2*i]->bond[1];
			else aa0=aaa[2*i]->bond[0];

			if(aaa[2*i+1]->tatm->id==0) bb0=aaa[2*i+1]->bond[1];
                        else bb0=aaa[2*i+1]->bond[0];

			de=de/dt;
			ad=ad/dt;
			DistPopular *t=dst->getDistPopular(0,aaa[2*i]->name,aaa[2*i+1]->name,0);
			
			float low=max(t->dist[0],ad-de);
			float high=min(t->dist[1],ad+de);
			
			dp[6*htt]=aaa[2*i]->id0+0.1;
			dp[6*htt+1]=aaa[2*i+1]->id0+0.1;
			dp[6*htt+2]=ad;
			dp[6*htt+3]=low;
			dp[6*htt+4]=high;
			dp[6*htt+5]=5+dt;

			bd=bd/dt;
			
			t=dst->getDistPopular(0,aa0->tatm->name,aaa[2*i+1]->tatm->name,0);
                        low=max(t->dist[0],bd-de);
                        high=min(t->dist[1],bd+de);
                        dp[6*htt]=aa0->id0+0.1;
                        dp[6*htt+1]=aaa[2*i+1]->id0+0.1;
                        dp[6*htt+2]=bd;
                        dp[6*htt+3]=low;
                        dp[6*htt+4]=high;
                        dp[6*htt+5]=5+dt;

			cd=cd/dt;
			t=dst->getDistPopular(0,aaa[2*i]->tatm->name,bb0->tatm->name,0);
                        low=max(t->dist[0],cd-de);
                        high=min(t->dist[1],cd+de);
                        dp[6*htt]=aaa[2*i]->id0+0.1;
                        dp[6*htt+1]=bb0->id0+0.1;
                        dp[6*htt+2]=cd;
                        dp[6*htt+3]=low;
                        dp[6*htt+4]=high;
                        dp[6*htt+5]=5+dt;

			dd=dd/dt;
			t=dst->getDistPopular(0,aa0->tatm->name,bb0->tatm->name,0);
                        low=max(t->dist[0],dd-de);
                        high=min(t->dist[1],dd+de);
                        dp[6*htt]=aa0->id0+0.1;
                        dp[6*htt+1]=bb0->id0+0.1;
                        dp[6*htt+2]=dd;
                        dp[6*htt+3]=low;
                        dp[6*htt+4]=high;
                        dp[6*htt+5]=5+dt;
                        htt++;
		}
	}

	if(htt==0) {
		delete [] dp;
		dp=0;
	}
        return dp;
}

int ModAlgnFmt::ishbondexist(Res *r, Res *r0) {
	
        ModAlgnFmt *t,*t0;
        for(t=this;t;t=t->next) {
		ModAlgnFmt *s=t->findsequencefmt();
                int n=s->compare[r->id0];
                int n0=s->compare[r0->id0];
                for(t0=t;t0;t0=t0->more) {
                        if(strcmp(t0->token,"structure")!=0) continue;
 
                        if(t0->seqngap[n]=='-'||t0->seqngap[n0]=='-') continue;

                        int i=t0->match[n];

                        int i0=t0->match[n0];
			
			Res *rr=t0->resn[i];
			Res *rr0=t0->resn[i0];
			if(rr==0||rr0==0) continue;
			
			if(rr->ishbondedwithres(rr0)||rr0->ishbondedwithres(rr)) return 1;
                }
        }
        return 0;
}

void ModAlgnFmt::setsecstruct() {

	ModAlgnFmt *t,*t0;
        for(t=this;t;t=t->next) {
		for(t0=t;t0;t0=t0->more) {
			if(strcmp(t0->token,"structure")!=0) continue;
			t0->pdb->chn->buildhbond();
			t0->pdb->chn->setdsspstr();
		}
	}
}

int  ModAlgnFmt::getstructurenumber(){
	
	ModAlgnFmt *t,*t0;

	int n=0;
        for(t0=this;t0;t0=t0->next) 
        for(t=t0;t;t=t->more) {
		if(strcmp(t->token,"structure")) continue;
		n++;
	}	
	return n;
}
