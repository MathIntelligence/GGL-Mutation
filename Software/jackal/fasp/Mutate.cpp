#include"source.h"

Mutate::Mutate()
{
	initial();
}

void Mutate::initial() {
	
	aloop=1;
	xyzself=0;
        refine=0;
	mpdb=0;
	sqnto=0;
	owner=0;
	next=0;
	flag=0;
	match=0;
	resn=0;
	compare=0;
	dssp=0;
	rely=0;
	segen=0;
	boundscore=0;
	pdbcopy=0;
	pdbold=0;
	fapr=0;
	seqn=0;
	resid=0;
	autotmp=0;
	autoalgn=1;
	stepnum=1;
	window=20;
	mold=0;

	test=0;
	sharp=1;
	out=1;
	asloop=0;
	rmsd=2.0;
	code=0;
	seglen=100; 	 
	refineid=0;
	//
	restraint=1;	
}

Mutate::Mutate(Mutate *t)
{
	initial();
	aloop=t->aloop;
	xyzself=0;
        refine=t->refine;
	mpdb=0;
	if(t->mpdb) {
		mpdb=new Pdb(t->mpdb);
		mpdb->configure();
	}
	sqnto=0;
	if(t->sqnto) sqnto=strdup(t->sqnto);
	owner=0;	
	next=0;
   	flag=t->flag;
	
	match=0;
	resn=0;
	compare=0;
	dssp=0;
	if(t->match&&t->compare&&t->dssp) {
		int nlen;nlen=strlen(sqnto);
		match=new int[nlen];
		compare=new int[nlen];
		resn=new Res*[nlen];
		dssp=new char[nlen];
		int i;
		for(i=0;i<nlen;i++) {
			match[i]=t->match[i];
			compare[i]=t->compare[i];
			dssp[i]=t->dssp[i];
			if(t->resn[i]==0) resn[i]=0;
			else {
				resn[i]=mpdb->chn->isres(t->resn[i]->id0);
			}
		}
	}

	

	rely=0;	

	segen=0;
	boundscore=0;
	pdbcopy=0;
	pdbold=0;
	fapr=t->fapr;
	seqn=0;
	if(t->seqn) seqn=strdup(t->seqn);
	resid=0;
	autotmp=t->autotmp;
	autoalgn=t->autoalgn;
	stepnum=t->stepnum;
	window=t->window;
	mold=0;

	test=t->test;
	sharp=t->sharp;
	out=t->out;
	asloop=t->asloop;
	rmsd=t->rmsd;
	code=0;
	if(t->code) {
		code=strdup(t->code);
	}
	seglen=t->seglen; 	 
	refineid=0;
	//
	restraint=t->restraint;	
}
Mutate::~Mutate()
{
	if(code) delete [] code;code=0;
	if(xyzself) delete []  xyzself;xyzself=0;
	if(mpdb&&(flag&TRES.constant->pdbundeletable)==0){ delete mpdb;mpdb=0;}
	if(sqnto) delete [] sqnto;sqnto=0;
	if(next) delete next;	next=0;
	if(match)delete [] match;match=0;
	if(resn) delete [] resn;resn=0;
	if(compare) delete [] compare;compare=0;
	if(dssp) delete [] dssp;dssp=0;
	if(rely) delete [] rely;rely=0;
	if(segen) delete segen;segen=0;
	if(boundscore) delete boundscore;boundscore=0;
	if(pdbcopy) delete pdbcopy;pdbcopy=0;
	if(pdbold) delete pdbold;pdbold=0;
	if(seqn) delete [] seqn;seqn=0;	 
	if(resid) delete [] resid;resid=0;
	mold=0;
}



void Mutate::setboundscore() {
	
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	int nlen;nlen=strlen(sqnto);
	
	int i;

	if(boundscore) delete [] boundscore;  boundscore=0;

	boundscore=new float[nlen];

	for(i=0;i<nlen;i++) boundscore[i]=100;
		
	for(i=0;i<nlen;i++) {

		if(owner->seqngap[i]=='-'&&se->seqngap[i]=='-') goto re200;
 		if(parent->sitescore[i]==0) goto re200;
		int n1,n2;

		for(n1=i-1;n1>=0;n1--) {
			if(owner->seqngap[n1]=='-'&&se->seqngap[n1]=='-') continue;
			if(owner->seqngap[n1]=='-'||se->seqngap[n1]=='-') break;
		}	 	

		for(n2=i+1;n2<nlen;n2++) {
			if(owner->seqngap[n2]=='-'&&se->seqngap[n2]=='-') continue;
			if(owner->seqngap[n2]=='-'||se->seqngap[n2]=='-') break;
		}	 	
		float xn1;xn1=n1-0.5;
		float xn2;xn2=n2+0.5;
		float a,x;a=0;x=0;

		if(n1==0||n2==nlen) {
			x=1./(i-xn1+0.5);			
			a+=2*pow(x,0.75);
			x=1./(xn2-i+0.5);		
			a+=2*pow(x,0.75);				
		}
		else {
			x=1./(i-xn1+0.25);			
			a+=2*pow(x,0.75);
			x=1./(xn2-i+0.25);		
			a+=2*pow(x,0.75);			
		}
						
		
		boundscore[i]=pow(1/parent->sitescore[i],0.5)*a;
		re200:
		if(TRES.logg>3) cerr<<i<<" bound ..."<<boundscore[i]<<" "<<se->seqngap[i]<<" "<<sqnto[i]<<" "<<endl; 
	}

	
}




void Mutate::setchiangle(int *stem) {

	
}

void Mutate::printsegment(Res **stem) {

	int out[2];
	out[0]=compare[stem[0]->id0];
	out[1]=compare[stem[1]->id0];
	printsegment(out);
}

void Mutate::printsegment(int *stem){

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

	if(stem==0) return;
	int m1=stem[0];
	int m2=stem[1];
	//int mm=stem[2];

	int i=0;
	
	for(i=m1;i<=m2;i++) {

		cerr<<dssp[i];
	}
	cerr<<endl;

	for(i=m1;i<=m2;i++) {

		cerr<<sqnto[i];
	}
	cerr<<endl;

	for(i=m1;i<=m2;i++) {

                cerr<<owner->seqngap[i];
        }
        cerr<<endl;

	for(i=m1;i<=m2;i++) {

		cerr<<se->seqngap[i];
	}
	cerr<<endl;
}

void Mutate::printsegment(FILE *fp,int *stem){

        StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

        if(stem==0) return;
        int m1=stem[0];
        int m2=stem[1];
        //int mm=stem[2];

        int i=0;

        for(i=m1;i<=m2;i++) {

                fprintf(fp,"%c",dssp[i]);
        }
	fprintf(fp,"\n");

        for(i=m1;i<=m2;i++) {
                fprintf(fp,"%c",sqnto[i]);
        }
	fprintf(fp,"\n");

        for(i=m1;i<=m2;i++) {

                fprintf(fp,"%c",owner->seqngap[i]);
        }
	fprintf(fp,"\n");

        for(i=m1;i<=m2;i++) {

                fprintf(fp,"%c",se->seqngap[i]);
        }
	fprintf(fp,"\n");
}




Res ** Mutate::reordertmp(Res **tmp) {
	
	int n=0;
	while(tmp&&tmp[n]) n++;
	if(n==0) return tmp;

	int i,j,m;
	while(1) {
		m=0;
		for(i=0;i<n;i++)  
		for(j=i;j<n;j++) {
			if(tmp[i]->id>tmp[j]->id) {
				Res *a=tmp[i];
				tmp[i]=tmp[j];
				tmp[j]=a;
				m++;
			}
		}
		if(m) continue;
		else break;
	}
	return tmp;
}

int Mutate::calcinsert(Res **tmp) {

	int n=0;
	while(tmp&&tmp[n]) {
		n++;
	}
	return n;
}

int Mutate::calcinsert() {

	int n=strlen(sqnto);
	 
	int i=0;
	
	Res *r;
	for(r=mpdb->chn->res;r;r=r->next) {
		if(rely[r->id0]==0) i++;
	}
	return i;
}





Res *Mutate::findminscore(Res *rr,int nn) {
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
	if(se==0) return 0;
  
	Res *r,*t=0;
	float a=100000;
	for(r=rr;r;r=r->next) {
		if(r->id0>nn) break;
		int n=compare[r->id0];
		if(parent->sitescore[n]<a) {
			a=parent->sitescore[n];
			t=r;
		}
		 
	}

	return t;
}
 
Res *Mutate::findminscore(int *status,Res *rr,int nn) {
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
	if(se==0) return 0;
  
	Res *r,*t=0;
	float a=100000;
	for(r=rr;r;r=r->next) {
		if(r->id0>nn) break;
		if(status[r->id0]==0) continue;
		if(rely[r->id0]==2&&sharp==2) continue;		
		int n=compare[r->id0];		 
		if(parent->sitescore[n]<a) {
			a=parent->sitescore[n];
			t=r;
		}		 
	}

	return t;
}




int Mutate::calcinsert(Res *rr,int nn) {
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
	if(se==0) return 0;
 
	int i=0;
	
	Res *r;
	for(r=rr;r;r=r->next) {
		if(r->id0>nn) break;
		int n=compare[r->id0];
		if(se->seqngap[n]!='-'&&owner->seqngap[n]=='-') i++;
	}
	return i;
}

int Mutate::calcdelete(Res *rr,int nn) {
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
	if(se==0) return 0;
 
	int i=0;
	int n1=compare[rr->id0];
	int n2=compare[nn];
	int n; 
	for(n=n1;n<=n2;n++) {				
		if(se->seqngap[n]=='-'&&owner->seqngap[n]!='-') i++;
	}
	return i;
}


int Mutate::calcnewres(Res *s,int n) {
		 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	Res *r;
	int mm=0;
	for(r=s;r;r=r->next) {
		if(r->id0>n) break;
		int j=compare[r->id0];
		if(owner->seqngap[j]=='-'&&se->seqngap[j]!='-') mm++;		 
	}
	return mm;
}

void Mutate::updatepdbcopy() {
	if(pdbcopy) delete pdbcopy;pdbcopy=0;
	pdbcopy=new Pdb(mpdb);
}

void Mutate::updatepdbold() {
	if(pdbold) delete pdbold;pdbold=0;
	pdbold=new Pdb(mpdb);
}

void Mutate::setsegbed(SegBed *segbed,int *stem) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;
	
	segbed->num=0;
	int m1,m2;
	m1=match[stem[0]];
	m2=match[stem[1]];
			
	if(m1==-1||m2==-1) return;
	Res *r1=resn[match[stem[0]]];
	Res *r2=resn[match[stem[1]]];
	
	int n1,n2;
	
	int insert=calcinsert();
 	segbed->add((float *)0,r1->id0,r2->id0);
	if(r1->last==0||r2->next==0||fapr==5) {	 	
		return;
	}
	
	Res *r;
	if(insert) {
		n1=r1->id0;
		for(r=r1;r;r=r->next) {
			if(rely[r->id0]==0) break;
			n1=r->id0;
		}
		n2=r2->id0;
		for(r=r2;r;r=r->last) {
			if(rely[r->id0]==0) break;
			n2=r->id0;
		}		
	}
	else {
		n1=r1->id0;
		for(r=r1;r;r=r->next) {
			if(r->next==0) break;
			if(r->next&&mpdb->chn->islinked(r,r->next)==0) break;
			n1=r->next->id0;
		}
		n2=r2->id0;
		for(r=r2;r;r=r->last) {
			if(r->last==0) break;
			if(r->last&&mpdb->chn->islinked(r->last,r)==0) break;
			n2=r->last->id0;
		}	 
	}

	int low=mpdb->chn->res->id0;
	int high=mpdb->chn->lastres()->id0;

	int p1,p2;
	p1=compare[n1];
	p2=compare[n2];
	if(p1==-1||p2==-1) return;
	if(parent->sitescore[p1]>parent->sitescore[p2]+0.2) {
		p1=2;p2=3;
	}
	else if(parent->sitescore[p1]<parent->sitescore[p2]-0.2){
		p1=3;p2=2;
	}
	else {
		p1=3;p2=3;
	}

	int i,j;
	
	for(i=n1;i>=n1-p1&&i>=low;i--)	
	for(j=n2;j<=n2+p2&&j<=high;j++) {
		//if((j-i+1-insert<3||j-i+1-insert>3)&&insert) continue;	
		//if(j-i+1-insert>3) continue;
		int a1=n1-i;
		int a2=j-n2;
		if(abs(a1-a2)>1) continue;
		if(j-i+1-insert!=3) continue;
		//if(j-i+1-insert<3||j-i+1-insert>4) continue;
		//if(j-i+1-insert!=3&&j-i+1-insert!=4) continue;
		if(segbed->ifexist(i,j)) continue;
		int t1=compare[i];
		int t2=compare[j];
		if(t1==-1||t2==-1) continue;
		if(dssp[t1]=='h'||dssp[t1]=='e') continue;
		if(dssp[t2]=='h'||dssp[t2]=='e') continue;
		segbed->add((float *)0,i,j);
	}
}


int *Mutate::findsegbed(SegBed *s,int done) {

	if(done>=s->num) return 0;
	 
	int *stem=new int[2];
	int n1=s->start[done];
	int n2=s->end[done];
	stem[0]=compare[n1];
	stem[1]=compare[n2];
	return stem;
}

void Mutate::putsegbed(SegBed *s) {
	
	int i;
	if(TRES.logg)mpdb->chn->write("semp.pdb");
	

	int start,end;
	start=1000000;end=-1000000;	
	for(i=0;i<s->num;i++) {	
		int n1=s->start[i];
		int n2=s->end[i];
		start=min(start,n1);
		end=max(end,n2);
	}
	
	Res *r0=mpdb->chn->isres0(start);
	if(r0==0) r0=mpdb->chn->res;
	Res *r1=mpdb->chn->isres0(end);
	if(r1==0) r1=mpdb->chn->lastres();
	start=r0->id0;
	end=r1->id0;

	s->first=start;
	s->last=end;
	float *xyzt=mpdb->chn->gettransfer(r0,end);

	for(i=0;i<s->num;i++) {
		mpdb->chn->transfer(xyzt,r0,end,1);
		int n1=s->start[i];
		int n2=s->end[i];
		Res *r1=mpdb->chn->isres0(n1);
		if(r1==0) continue;
		Res *te1=mpdb->chn->isres0(n2);
		if(te1==0) continue;
		float *xyz=s->xyzout[i];
		if(xyz==0) continue;		
		mpdb->chn->transfer(xyz,r1,n2,1);			
		s->xyzout[i]=mpdb->chn->gettransfer(r0,end);
		delete [] xyz;xyz=0;		
	}
	s->xyzout[s->num]=0;
	mpdb->chn->transfer(xyzt,r0,end,1);
	//segen->start=start;
	//segen->end=end;
	segen->resort(s);
	delete [] xyzt;xyzt=0;
}
void Mutate::putsegbedmean(SegBed *s) {
	
	int i;
	if(TRES.logg)mpdb->chn->write("semp.pdb");
	

	int start,end;
	start=1000000;end=-1000000;	
	for(i=0;i<s->num;i++) {	
		int n1=s->start[i];
		int n2=s->end[i];
		start=min(start,n1);
		end=max(end,n2);
	}
	
	Res *r0=mpdb->chn->isres0(start);
	if(r0==0) r0=mpdb->chn->res;
	Res *r1=mpdb->chn->isres0(end);
	if(r1==0) r1=mpdb->chn->lastres();
	start=r0->id0;
	end=r1->id0;

	s->first=start;
	s->last=end;
	float *xyzt=mpdb->chn->gettransfer(r0,end);

	for(i=0;i<s->num;i++) {
		mpdb->chn->transfer(xyzt,r0,end,1);
		int n1=s->start[i];
		int n2=s->end[i];
		Res *r1=mpdb->chn->isres0(n1);
		if(r1==0) continue;
		Res *te1=mpdb->chn->isres0(n2);
		if(te1==0) continue;
		float *xyz=s->xyzout[i];
		if(xyz==0) continue;		
		mpdb->chn->transfer(xyz,r1,n2,1);			
		s->xyzout[i]=mpdb->chn->gettransfer(r0,end);
		delete [] xyz;xyz=0;		
	}
	s->xyzout[s->num]=0;
	mpdb->chn->transfer(xyzt,r0,end,1);
	//segen->start=start;
	//segen->end=end;
	segen->resortmean(s);
	delete [] xyzt;xyzt=0;
}

void Mutate::unitesegbed(SegBed *s) {
	
	int i;
	if(TRES.logg)mpdb->chn->write("semp.pdb");
	
	int start,end;
	start=1000000;end=-1000000;	
	for(i=0;i<s->num;i++) {	
		int n1=s->start[i];
		int n2=s->end[i];
		start=min(start,n1);
		end=max(end,n2);
	}
	
	Res *r0=mpdb->chn->isres0(start);
	if(r0==0) r0=mpdb->chn->res;
	Res *r1=mpdb->chn->isres0(end);
	if(r1==0) r1=mpdb->chn->lastres();
	start=r0->id0;
	end=r1->id0;

	s->first=start;
	s->last=end;
	float *xyzt=mpdb->chn->gettransfer(r0,end);

	for(i=0;i<s->num;i++) {
		mpdb->chn->transfer(xyzt,r0,end,1);
		int n1=s->start[i];
		int n2=s->end[i];
		Res *r1=mpdb->chn->isres0(n1);
		if(r1==0) {
			if(s->xyzout[i]) delete [] s->xyzout[i];
			s->xyzout[i]=0;
			continue;
		}
		Res *te1=mpdb->chn->isres0(n2);
		if(te1==0) {
			if(s->xyzout[i]) delete [] s->xyzout[i];
			s->xyzout[i]=0;
			continue;
		}
		float *xyz=s->xyzout[i];
		if(xyz==0) continue;		
		mpdb->chn->transfer(xyz,r1,n2,1);			
		s->xyzout[i]=mpdb->chn->gettransfer(r0,end);
		if(xyz) delete [] xyz;xyz=0;		
	}
	s->xyzout[s->num]=0;
	mpdb->chn->transfer(xyzt,r0,end,1);
	//delete [] xyzt;xyzt=0;

	for(i=0;i<s->num;i++) {	
		s->start[i]=s->first;
		s->end[i]=s->last;
	}

	int nn=0;
	for(i=0;i<s->num;i++) {
		if(s->xyzout[i]==0) continue;
		s->xyzout[nn++]=s->xyzout[i];
	}
	s->xyzout[nn]=0;
	s->num=nn;

	nn=calcinsert(r0,end)*10;

	if(s->num>nn*10) {
		segen->start=s->first;
		segen->end=s->last;
		segen->pdb=mpdb;
		segen->cid=mpdb->chn->id;

		if(segen->disc) delete segen->disc;segen->disc=0;
		segen->disc=new Disc;
		segen->disc->setupall(mpdb,segen->cutoff,1);
		strcpy(segen->dforce,"u");
		float *ent=segen->ensort(s->xyzout);
		if(ent) delete [] ent;ent=0;
 	}

	nn=calcinsert(r0,end)*10;
	for(i=nn;i<s->num;i++) {
		delete [] s->xyzout[i];
		s->xyzout[i]=0;
	}
	
	//float avgrmsd=segen->calcavgrmsd(s->xyzout,20);
	//segen->noclose(s->xyzout,avgrmsd/2);

	nn=0;
	for(i=0;i<s->num;i++) {
		if(s->xyzout[i]==0) continue;
		s->xyzout[nn++]=s->xyzout[i];
	}
	s->xyzout[nn]=0;
	s->num=nn;	
	
	//s->size=s->num+1;
	
	int nlen=strlen(owner->seqngap);
	if(s->post) delete s->post;s->post=0;
	s->post=new Res*[nlen];
	for(i=0;i<nlen;i++)s->post[i]=0;
	nn=0;Res *r;
	for(r=r0;r;r=r->next) {
		if(r->id0>end) break;
		s->post[nn++]=r;			
	}
	int nt=0;
	while(s->head&&s->head[nt])nt++;
	for(i=0;i<nt;i++) {
		if(s->head[i]) delete [] s->head[i];
		s->head[i]=0; 
	}
	for(i=0;i<s->num;i++) {
		char line[100];
		sprintf(line,"semp.%i.pdb",i);
		mpdb->chn->transfer(s->xyzout[i],r0,end,1);
		mpdb->chn->header(r0,end+1);
		s->head[i]=mpdb->chn->gettransfertemp(r0,end);
		if(TRES.logg) mpdb->chn->write(line);
	}
	mpdb->chn->transfer(s->xyzout[0],r0,end,1);
	mpdb->chn->transfertemp(s->head[0],r0,end,1);
	//mpdb->chn->transfer(xyzt,r0,end,1);
	if(xyzt) delete [] xyzt;xyzt=0;
	if(TRES.logg) mpdb->chn->write("temp.pdb");
	//mpdb->chn->header(r0,end+1);
	return;
}

void Mutate::setnoclose(float **xyzout,int a,int b) {
	
	int nn=0;
	if(xyzout==0) return;
	while(xyzout[nn])nn++;
 	if(nn<2) return;

	segen->start=a;
	segen->end=b;
	segen->pdb=mpdb;
	segen->cid=mpdb->chn->id;
 
	//float avgrmsd=segen->calcavgrmsd(xyzout,20);
	//float avgrmsd=0.1;
	segen->noclose(xyzout,0.2);
	segen->nolarge(xyzout,1.0);
	nn=0;while(xyzout[nn])nn++;

	if(nn<2) return;
 	/*
	if(segen->disc) delete segen->disc;segen->disc=0;
	segen->disc=new Disc;
	segen->disc->setupall(mpdb,segen->cutoff,1);
	*/ 
	strcpy(segen->dforce,"udDW");
	float *ent=segen->ensort(xyzout);
	int i;
	for(i=1;i<nn;i++) {
		if(ent[i]-ent[0]>5) {
			delete [] xyzout[i];
			xyzout[i]=0;
		}
	}
	if(ent) delete [] ent;ent=0;

	int j=0;
	for(i=0;i<nn;i++) {
		if(xyzout[i]==0) continue;
		xyzout[j++]=xyzout[i];
	}	
	xyzout[j]=0;
}
void Mutate::resetsegbed(SegBed *s) {
	
	int nn=0;
	while(s->post&&s->post[nn])nn++;
	if(nn==0) {
		s->num=0;return;
	} 

	//Res *r0=s->post[0];	
	//int end=s->post[nn-1]->id0;
	 
	Res *r;
	int i,j,k,ii;
	
	int end0=mpdb->chn->lastres()->id0;

	float *xyz=mpdb->chn->gettransfer(mpdb->chn->res,100000);
	float *xyztemp=mpdb->chn->gettransfertemp(mpdb->chn->res,100000);
	for(k=0;k<s->num;k++) {
		i=0;ii=0;
		mpdb->chn->transfer(xyz,mpdb->chn->res,100000,1);
		mpdb->chn->transfertemp(xyztemp,mpdb->chn->res,100000,1);
		for(j=0;j<nn;j++) {
			r=s->post[j];
 			r->transfer(s->xyzout[k]+i,1);
 			i+=r->tres->number*3;
			r->transfertemp(s->head[k]+ii,1);
			ii+=9;
		}
		if(s->xyzout[k])delete [] s->xyzout[k];s->xyzout[k]=0;
		if(s->head[k]) delete [] s->head[k];s->head[k]=0;
		s->xyzout[k]=mpdb->chn->gettransfer(mpdb->chn->res,100000);
		s->head[k]=mpdb->chn->gettransfertemp(mpdb->chn->res,100000);
		
		//mpdb->chn->transfer(s->xyzout[k],r0,end,0);
		char line[100];
		sprintf(line,"temp.%i.pdb",k);
		if(TRES.logg) mpdb->chn->write(line);
	}

	for( k=0;k<s->num;k++) {
		s->start[k]=mpdb->chn->res->id0;
		s->end[k]=end0;
	}	
	s->first=mpdb->chn->res->id0;
	s->last=end0;
	mpdb->chn->transfer(xyz,mpdb->chn->res,100000,1);
	mpdb->chn->transfertemp(xyztemp,mpdb->chn->res,100000,1);
	if(xyz) delete [] xyz;xyz=0;
	
	return;
}



int *Mutate::getnewstem(int *stem) {
	 
	int m1,m2;
	m1=match[stem[0]];
	m2=match[stem[1]];	
	if(m1==-1||m2==-1) return stem;
	Res *r1=resn[m1];
	Res *r2=resn[m2];
	
	if(r1==0||r2==0) return stem;
 
	for(r1=r1;r1->last;r1=r1->last) {
		if(rely[r1->id0]==2) break;
	}
	 
	for(r2=r2;r2->next;r2=r2->next) {
		if(rely[r2->id0]==2) break;
	} 

	m1=compare[r1->id0];
	m2=compare[r2->id0];

	stem[0]=m1;
	stem[1]=m2;
	return stem;
}

void Mutate::checkburiedpolars() {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	cerr<<"calculate solvent accessible area for each residue of the template..."<<endl;
	LatSurfv latsurfv;
	latsurfv.ready(owner->pdb,2,1.4);
	latsurfv.calcarea(owner->pdb);
	latsurfv.assignarea(owner->pdb);	

	Chn *c;
	Res *r;
	Atm *a;
	int nbt=0;
	for(c=owner->pdb->chn;c;c=c->next)
	for(r=c->res;r;r=r->next) {
		float e=r->bury(0,100);
		if(e<0) e=0;
		e=(1-e)*100;
		float ee=r->bury(4,100);
                if(ee<0) ee=0;
                ee=(1-ee)*100;
		if(r->name!='G') 
		fprintf(stderr,"residue: %c%i total exposed area(A^2): %8.3f  percent exposed of sidechain(all):%i%c(%i%c)\n",
		r->tres->name,r->oid,r->area,int(ee),'%',int(e),'%');
		else 
		fprintf(stderr,"residue: %c%i total exposed area(A^2): %8.3f  percent exposed of sidechain(all):*(%i%c)\n",
                r->tres->name,r->oid,r->area,int(e),'%');
		/*
		int bry=1;
		for(a=r->atm;a;a=a->next) {
			if(a->tatm->id<=4) continue;
			if(a->tatm->name[1]!='O'&&a->tatm->name[1]!='N') continue;
			if(a->area>a->tatm->area*0.1) {bry=0;break;}
		}
		*/
		if(ee<10&&strchr("KRDE",r->name))  {
			nbt++;
		}
	
	}
	cerr<<"the total buried charges: "<<nbt<<endl;
	cerr<<"the total solvent exposed area of the template is: "<<owner->pdb->area<<" A^2"<<endl;
	cerr<<"================================="<<endl;
	cerr<<endl;
	cerr<<"check salt bridges and disulfide bond in the template..."<<endl;
	for(c=owner->pdb->chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {
		if(strchr("KRDEC",r->name)==0) continue;
		char u=getmodelresname(r);
		Res *rr;
		for(rr=r->next;rr;rr=rr->next) {
			if(strchr("KRDEC",rr->name)==0) continue;	
			if(r->name=='C'&&rr->name!='C') continue;
			if(r->name!='C'&&rr->name=='C') continue;
			char uu=getmodelresname(rr);
			int cs=0;
			if(u==r->name&&uu==rr->name) cs=1;
			Atm *a,*aa;
			for(a=r->atm;a;a=a->next) {
				if(a->tatm->id==0||a->tatm->id==3) continue;
				if(a->tatm->name[1]!='N'&&a->tatm->name[1]!='O'&&a->tatm->name[1]!='S') continue;
				for(aa=rr->atm;aa;aa=aa->next) {
					if(aa->tatm->id==0||aa->tatm->id==3) continue;
					if(aa->tatm->name[1]!='N'&&aa->tatm->name[1]!='O'&&aa->tatm->name[1]!='S') continue;
					float d=TRES.distance(a,aa);
					if(d>3.5) continue;
					if(r->name=='C'&&d>3.0) continue;
					if(r->name=='C') { 
					cerr<<"disulfide bond found in the template: " <<r->name<<r->oid<<"("<<a->name<<")"<<"--";
					}
					else { 
					cerr<<"salt bridge found in the template: " <<r->name<<r->oid<<"("<<a->name<<")"<<"--";
					}
					if(cs==1) 
						cerr<<rr->name<<rr->oid<<"("<<aa->name<<")"<<" : "<<d<<"A. conserved!"<<endl;
					else 
						cerr<<rr->name<<rr->oid<<"("<<aa->name<<")"<<" : "<<d<<"A."<<endl;
					goto re200;
				}
			}
			re200:
			continue;
		}
	}
}


Res **Mutate::checkmodelburiedpolars(int ot,int srf) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	cerr<<"calculate solvent accessible area for each residue of the model..."<<endl;
	LatSurfv latsurfv;
	if(srf==1) {
		latsurfv.ready(mpdb,2,1.4);
		latsurfv.calcarea(mpdb);
		latsurfv.assignarea(mpdb);	
	}
	//
	int maxn=mpdb->chn->lastres()->id0*2+1000;
	Res **badres=new Res*[maxn];
	int i;
	for(i=0;i<maxn;i++) badres[i]=0;
	int nbad=0;
	//

	Chn *c;
	Res *r;
	Atm *a;
	 
	for(c=mpdb->chn;c;c=c->next)
	for(r=c->res;r;r=r->next) {
		float e=r->bury(0,100);
		if(e<0) e=0;
		e=(1-e)*100;
		float ee=r->bury(4,100);
                if(ee<0) ee=0;
                ee=(1-ee)*100;
		if(r->name!='G')  {
			if(ot) fprintf(stderr,"residue: %c%i total exposed area(A^2): %8.3f percent exposed of sidechain(all):%i%c(%i%c)\n", r->tres->name,r->oid,r->area,int(ee),'%',int(e),'%');
		}
		else {
			if(ot) fprintf(stderr,"residue: %c%i total exposed area(A^2): %8.3f  percent exposed of sidechain(all):*(%i%c)\n", r->tres->name,r->oid,r->area,int(e),'%');
		}
		/* 
		int bry=1;
		for(a=r->atm;a;a=a->next) {
			if(a->tatm->id<=4) continue;
			if(a->tatm->name[1]!='O'&&a->tatm->name[1]!='N') continue;
			if(a->area>a->tatm->area*0.1) {
				bry=0;
				break;
			}
		}
		*/
		if(ee<10&&strchr("KRDE",r->name)) {
			if(TRES.logg>3) cerr<<r->name<<r->oid<<" "<<ee<<endl;
			badres[nbad++]=r;	
		}
	}
	int nbury=nbad;
	
	cerr<<"the total solvent exposed area of the model is: "<<mpdb->area<<" A^2"<<endl;
	cerr<<"================================="<<endl;
	cerr<<endl;
	cerr<<"check salt bridges and disulfide bond in the model..."<<endl;
	int nonc=0;
	for(c=owner->pdb->chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {
		if(strchr("KRDEC",r->name)==0) continue;
		Res *u=getmodelres(r);
		Res *rr;
		for(rr=r->next;rr;rr=rr->next) {
			if(strchr("KRDEC",rr->name)==0) continue;	
			if(r->name=='C'&&rr->name!='C') continue;
			if(r->name!='C'&&rr->name=='C') continue;
			Res *uu=getmodelres(rr);
			int cs=0;
			if(uu&&u&&u->name==r->name&&rr->name==uu->name) cs=1;
			if(cs==0) continue;
			if(isionbridge(r,rr)==0) continue;
			if(isionbridge(u,uu)==0) {
				if(r->name=='C') { 
					if(ot) cerr<<"disulfide bond in the template not conserved: " <<u->name<<r->oid<<"--";
				}
				else { 
					if(ot) cerr<<"salt bridge in the template not conserved: " <<u->name<<r->oid<<"--";
				}
					if(ot) cerr<<uu->name<<uu->oid<<endl;
				nonc++;
                                badres[nbad++]=u;
                                badres[nbad++]=uu;
			} 
			else if(isionbridge(u,uu)) {
				if(r->name=='C') {
                                        if(ot) cerr<<"disulfide bond in the model conserved: " <<u->name<<u->oid<<"--";
                                }
                                else {
                                        if(ot) cerr<<"salt bridge in the model conserved: " <<u->name<<u->oid<<"--";
                                }
                                if(ot) cerr<<uu->name<<uu->oid<<endl;
			}
		}
	}
	cerr<<"the number of buried charges in the model: "<<nbury<<endl;
	cerr<<"the number of lost conserved salt bridges or disulfide bond:"<<nonc<<endl;
	if(nbad==0) {
		delete [] badres;
		badres=0;
	}
	return badres;
}

int Mutate::isionbridge(Res *r,Res *rr) {

	Atm *a,*aa;
	for(a=r->atm;a;a=a->next) {
		if(a->tatm->id==0||a->tatm->id==3) continue;
		if(a->tatm->name[1]!='N'&&a->tatm->name[1]!='O'&&a->tatm->name[1]!='S') continue;
		for(aa=rr->atm;aa;aa=aa->next) {
			if(aa->tatm->id==0||aa->tatm->id==3) continue;
			if(aa->tatm->name[1]!='N'&&aa->tatm->name[1]!='O'&&aa->tatm->name[1]!='S') continue;
			float d=TRES.distance(a,aa);
			if(d>3.5) continue;
			if(r->name=='C'&&d>3.0) continue;
			return 1;
			
		}
	}
	return 0;
}


void Mutate::actmutate0() {

	checkburiedpolars();

	if(segen) delete segen;segen=0;
        segen=new Segen();
	segen->flex=1;
	segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
	segen->randcoil=0;
	//TRES.smoothclash=1;
	segen->smoothclash=1; 
	segen->part=0;
 

	//segbed define 	
	SegBed segbed;
	segbed.next=new SegBed();
	segbed.setsize(1000);
	segbed.next->setsize(1000);
		
	segbed.next->next=new SegBed();
	segbed.next->next->setsize(1000);
    	
	//refine
	refine=0;

	StrFmt *root=owner->getrootStrFmt();
	if(root->onlyrefine==0) {
		cerr<<endl;
		cerr<<"model building for the "<<owner->getblockids()<<"th alignment blocks"<<endl;
		cerr<<"this step may take minutes to hours depending on your alignment and size.."<<endl;
		cerr<<"please be patient..."<<endl;
		cerr<<endl; 
	}

   	Strhandler cc;  
	
	FILE *fp=0;
	if(TRES.logg)fp=stderr;
	if(TRES.logg)mpdb->write("start.pdb");
        StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;
		
	int dloop=0;
        int nlen=strlen(sqnto);
        mpdb->setlastresorder();
	
 	int m=0;
	int *stem=0;
	if(fapr<3)setrotamerorder(0); 
	else      setrotamerorder(1);
	int n=actsidechainmutation(100);
	setrotamerorder(1); 
	//int n=actsidechainmutation(-99999);
	//minimize(mpdb->chn->res,100000,-100000); 
	//mpdb->chn->setbackbonetorsion(mpdb->chn->res,10000);
	
	int insert=0;//int newres=0;
	int isbrk=0; int midr=0;
	int rstem[2];
	rstem[0]=rstem[1]=-1;
	int dstem[2];
	dstem[0]=dstem[1]=-1;
	//int mloop=0;
	int asloop0=asloop;
	int mtt=0;
	int toted=0;
	while(1) {
		Res **tmp;	
		segen->secd='-';

		//m=findbadsite();
		//dloop
		if(asloop==0) dstem[0]=dstem[1]=-1;
		if(dstem[0]==-1||dstem[1]==-1)  m=findbadsite(); //find the best insertion and deletion region,most conserved
		else m=findbadsite(dstem);
		if(calctotalinsert(m)<3) asloop=0;
		else asloop=asloop0;
		// 
		if(m==-1) break;	
		midr=m;	
		//int tlen=getlength(m); tlen=2;
		int tlen=seglen;	
		//if(fapr==5) tlen=1000;
		isbrk=0; 	 
		if(isbreak(m)) { tlen=1000; isbrk=1;}
		//
		if(isbrk) findrealstem(m,rstem);
		//
		//mpdb->chn->header();
		tmp=setinitialstructure(m,tlen);	
		tmp=reordertmp(tmp);			 
        	mpdb->setlastresorder();
		linkstructureviaboth(tmp);
		//updatepdbcopy();
		//mpdb->chn->dihedral();		
		setreliable(tmp);				
		insert=calcinsert(); 
		//if(insert>=3) segen->randcoil=0;
		//else	      segen->randcoil=0;
		char lite[100];
		sprintf(lite,"thead.%i.pdb",mtt++);
		if(TRES.logg>3)mpdb->chn->write(lite);  
 		if(tmp==0)   stem=getsegmentworkon(m);		
		else 	     stem=getsegmentworkon(tmp,isbrk);
		if(TRES.logg>3)mpdb->chn->write("s");  
		//missing residues
		int ntmp=0;while(tmp&&tmp[ntmp])ntmp++;
		int it=0;
		if(rstem[0]==-1||rstem[1]==-1) ntmp=0;
		for(it=0;it<ntmp;it++) {
			int iu=compare[tmp[it]->id0];
			if(iu<rstem[0]||iu>rstem[1]) {
				rstem[0]=rstem[1]=-1;
			}
		}
		if(ntmp==0) rstem[0]=rstem[1]=-1;
		 

		//dloop
		if(tmp==0) dloop=0;
		else if(asloop)  {
			dloop=checkdloop(m);
			if(dloop==1&&segbed.num==0) dloop=1;
			else if(dloop==1) dloop=2;
			else if(dloop==0) dloop=3;
		}			 
		else dloop=0;
		if(dloop==0||dloop==3) dstem[0]=dstem[1]=-1;
		else	     setdloopstem(m,dstem); 
		if(asloop&&tmp) resetsegbed(segbed.next->next);		
		//		
		if(TRES.logg>3)mpdb->chn->write("s");  
		segbed.next->num=0;
		segbed.num=0;
		setsegbed(segbed.next,stem);	
		//if(asloop&&tmp) resetsegbed(segbed.next,segbed.next->next);	
		 
		if(tmp) delete [] tmp;tmp=0;
		int done=1;
 		if(TRES.logg>3)mpdb->chn->write("s");  
		float *xyzorg=mpdb->chn->gettransfer(mpdb->chn->res,1000000); 
		//dloop
		if(segbed.next->next->num==0&&dloop!=0) {				
				int i=0;
				int first=10000;
				int last=-1;
				for(i=0;i<segbed.next->num;i++) {
					first=min(segbed.next->start[i],first);
					last=max(segbed.next->end[i],last);
				}
				Res *rr1=mpdb->chn->isres0(first);
				float *xyzout=mpdb->chn->gettransfer(rr1,last);
				float *head=mpdb->chn->gettransfertemp(rr1,last);
				segbed.next->next->add(xyzout,first,last);
				segbed.next->next->head[0]=head;
				segbed.next->next->head[1]=0;
				segbed.next->next->first=first;
				segbed.next->next->last=last;
				//unitesegbed(segbed.next->next);				
		}
		//
		if(TRES.logg)mpdb->chn->write("s");  
		while(1) {		
			//stem=getnewstem(stem);
			int m1,m2;
			m1=match[stem[0]];
			m2=match[stem[1]];
			if(segbed.ifexist(m1,m2)) {			
				//done++;
				if(done==segbed.next->num) goto redone; 
				stem=findsegbed(segbed.next,done);
				done++;
				continue;
			}
			if(m1==-1||m2==-1) break;
			Res *r1;r1=resn[match[stem[0]]];
			Res *r2;r2=resn[match[stem[1]]];
		 	if(r1&&r2&&insert==0&&mpdb->chn->isalllinked(r1,r2->id0)==1) break;  		
			
			//newres=calcnewres(r1,r2->id0); 			

			if(segen->chiangle) {
                        	delete segen->chiangle;
                        	segen->chiangle=0;
                	}
                	else if(segen->bound) {
                        	delete segen->bound;
                        	segen->bound=0;
                	}						
			//setchiangle(stem);
			
			if(TRES.logg) printalign(fp,stem);
			int nn;nn=getreslength(stem); 
			if(TRES.logg) printsegment(stem);
			delsidechain(stem);
			
 	 	 	if(fapr==5) { 		
				if(insert>=3) segen->arbt=20/4;
                        	else if(insert>=2) segen->arbt=10/4;
                        	else if(insert>=1) segen->arbt=5/4;			
				else segen->arbt=5/4;
			}			
			else if(fapr==4) { 		
				if(insert>=3) segen->arbt=20/2;
                        	else if(insert>=2) segen->arbt=10/2;
                        	else if(insert>=1) segen->arbt=5/2;			
				else segen->arbt=5/2;
			}
			else if(fapr==3) { 		
				if(insert>=3) segen->arbt=20;
                        	else if(insert>=2) segen->arbt=10;
                        	else if(insert>=1) segen->arbt=5;			
				else segen->arbt=5;
			}
			else if(fapr==2) {
				if(insert>=3) segen->arbt=20*2;
                        	else if(insert>=2) segen->arbt=10*2;
                        	else if(insert>=1) segen->arbt=5*2;			
				else segen->arbt=5*2;				
			}
			else if(fapr==1) {
				if(insert>=3) segen->arbt=20*4;
                        	else if(insert>=2) segen->arbt=10*4;
                        	else if(insert>=1) segen->arbt=5*4;			
				else segen->arbt=5*4;				
			}
			else if(fapr==0) {
				if(insert>=3) segen->arbt=20*8;
                        	else if(insert>=2) segen->arbt=10*8;
                        	else if(insert>=1) segen->arbt=5*8;			
				else segen->arbt=5*8;				
			}
			
			//
 			if(segbed.num&&rstem[0]!=-1&&rstem[1]!=-1) {				
				if(stem[0]<rstem[0]||stem[1]>rstem[1]) {
					if(done==segbed.next->num) goto redone; 
					stem=findsegbed(segbed.next,done);
					done++;
					continue;		
				}		
			}
			// 
			toted++;
			cerr<<".";
			if(toted%50==0) cerr<<endl;
			
			if(asloop==0||dloop==0) {
				if(TRES.logg>3)mpdb->write("s1");				 
				float *xyzout=segen->myfixsegment(stem);
				if(TRES.logg>3)mpdb->write("s");
				if(xyzout&&TRES.logg>3) mpdb->chn->transfer(xyzout,r1,segen->end,1);
			 	if(TRES.logg>3)mpdb->write("ss1");
				mpdb->chn->transfer(xyzorg,mpdb->chn->res,1000000,1);
			 	if(TRES.logg>3)mpdb->write("ss");
				if(stem) {delete [] stem; stem=0;}
				if(xyzout) {
					segbed.add(xyzout,segen->start,segen->end);
					xyzout=0;
					//done++;
					if(done==segbed.next->num) goto redone; 
					stem=findsegbed(segbed.next,done);done++;
					continue;
				}
				else {
					stem=findloosesegment(midr,r1,r2);
					continue;
				}
			}
			else {				
				int it=0;
				if(TRES.logg>3)mpdb->write("se1");	
				Res *rr0=mpdb->chn->isres0(segbed.next->next->first);
				int gett=0;
				for(it=0;it<segbed.next->next->num;it++) {
					mpdb->chn->transfer(segbed.next->next->xyzout[it],rr0,segbed.next->next->last,1);
					mpdb->chn->transfertemp(segbed.next->next->head[it],rr0,segbed.next->next->last,1);
					if(TRES.logg)mpdb->write("s2");
					float **xyzout=segen->newfixsegment(stem);
					//float *xyzout=segen->myfixsegment(stem);
					if(TRES.logg>3)mpdb->write("s");
					if(xyzout&&xyzout[0]) mpdb->chn->transfer(xyzout[0],r1,segen->end,1);
			 		if(TRES.logg>3)mpdb->write("ss1");
					mpdb->chn->transfer(xyzorg,mpdb->chn->res,1000000,1);
			 		if(TRES.logg>3)mpdb->write("ss");			
					if(xyzout) {
						gett++;
						int a=segen->start;
						int b=segen->end;
						setnoclose(xyzout,segen->start,segen->end);
						segen->start=a;
						segen->end=b;
						segbed.add(xyzout,segen->start,segen->end);
						xyzout=0;
						//done++;						
					}
					
				}
				if(stem) {delete [] stem; stem=0;}
				if(gett==0) {
					stem=findloosesegment(midr,r1,r2);
					continue;
				}
				if(segbed.num) {										 
					if(done==segbed.next->num) goto redone; 
					stem=findsegbed(segbed.next,done);done++;
					continue;
				}
				else {					
					stem=findloosesegment(midr,r1,r2);
					continue;
				}				
			}
			redone:
			if(done==segbed.next->num) {	
				//modified on 12/07/2002 in case of no side-chain assembled. 
				hooksidechain();
				//end
			 	if(dloop==0||dloop==3) {
					putsegbed(&segbed);
					segbed.next->next->clear();
				}
				else  {
					segbed.next->next->clear();
					segbed.next->next->transfer(&segbed);
					unitesegbed(segbed.next->next);
				}	

				if(fapr<3){
					setrotamerorder(0); 
					actsidechainmutation(100);
					setrotamerorder(1);	
				}				 				 
				//minimize(mpdb->chn->res,100000,-100000); 				
				break;	
			}						 			 
		}
		if(xyzorg) delete [] xyzorg;xyzorg=0; 
	}		
	cerr<<endl;
	segen->mutate=0;
	
	asloop=asloop0;
	
	if(fapr<3) setrotamerorder(0); 
	else 	   setrotamerorder(1);
        actsidechainmutation(-99999);
	setrotamerorder(1);
	
	setlabel(mpdb);

	sethbonddssp();
	updatepdbold();
	reindex(pdbold);

	char line[100];

	reindex();
	PdbFix pdbfix;
	int nu=pdbfix.checkpdb(mpdb);
	if(nu) {
                        cerr<<"the initial model has missing atoms..."<<endl;
                        cerr<<"try profix program to fix it..."<<endl;
                        pdbfix.pdb=mpdb;
                        pdbfix.onlybackbone=0;
                        pdbfix.pdb->configure();
                        pdbfix.ready();
                        pdbfix.myfixnow();
                        mpdb=pdbfix.pdb;
                        pdbfix.pdb=0;
                        if(TRES.logg>3) mpdb->write("fix.pdb");
         }

	if(root->onlyrefine==0) {
	cerr<<endl<<"finish initial model building..."<<endl;
	if(code&&out) {		 		
		sprintf(line,"%s_ini_model.%i.pdb",code,refineid++);
		cerr<<endl<<"initial model written to: "<<line<<endl;
		cerr<<"model before any minimization..."<<endl<<endl;
		FILE *fpp=fopen(line,"w");
		fprintf(fpp,"!initial built model before minimizaton\n");
		StrFmt *parent=owner->getparentStrFmt();
        	StrFmt *se = parent->getsequenceStrFmt();
		char c=mpdb->chn->id;
		mpdb->chn->id=se->cid;
		if(mpdb->chn->id=='0') mpdb->chn->id=' ';
		mpdb->write(fpp);
		fclose(fpp);		
		mpdb->chn->id=c;
 	} 
	}
}

void Mutate::hooksidechain() {

	Res *r;
	int n=0;
	for(r=mpdb->chn->res;r;r=r->next) {
		if(r->issidechaincomplete()==1) continue;
		r->addsidechain();
		n++;
	}
	if(n) mpdb->configure();
}
void Mutate::setlabel(Pdb *s) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();

	Res *r;
	float t=0,tt=1;
	
	for(r=s->chn->res;r;r=r->next) {
		int i=compare[r->id0];
		int n=owner->match[i];
		if(n==-1) {
			t=99;
		} 
		else {
			Res *rr=owner->resn[n];
			float d=rr->directrmsdanyway(r,0,3);
			t=d*5;			
			if(t>50) t=50;			
		}	
		Atm *a;
		for(a=r->atm;a;a=a->next){
			a->bfact=t;
			a->occup=tt;
		} 
	}
}

void Mutate::optimize() {

	Res *refstart=mpdb->chn->res;
	int  refend=100000; 
 
	StrFmt *parent=owner->getparentStrFmt();
	StrFmt *root=owner->getrootStrFmt();

	if(root->refstart>-90000&&root->refend<90000)  {	
		refstart=mpdb->chn->isresoid(root->refstart);
		if(refstart==0) {
			cerr<<endl<<"no residue found with id: "<<root->refstart<<endl<<endl;
			exit(0);
		}
		refend=min(root->refend,mpdb->chn->lastres()->oid);
		Res *r=mpdb->chn->findsmallres(refend);
		if(r) refend=r->id0;
		else {
			cerr<<endl<<"no residue found with id: "<<root->refend<<endl<<endl;
			exit(0);
		} 
	}
	

	if(segen) {
		segen->mutate=0;
		delete segen;segen=0;
	}
        segen=new Segen();
	segen->flex=1;
	segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
	segen->randcoil=0;
	//TRES.smoothclash=1;
	segen->smoothclash=0; 
	segen->part=0;
  	segen->onlyenergy=0;
 	
	cerr<<"minimizing only unconserved residues..."<<endl;	
	
	//int ny=fapr;fapr=2;
	if(fapr<=3) minbadsidechain(refstart,refend,0);
	else minbadsidechain(refstart,refend,1);
	//fapr=ny;

	int nt;
	
 	for(nt=0;nt<1;nt++) {
		//minimize(mpdb->chn->res,100000,0); 		
		cerr<<endl;
		cerr<<"minimizing clashes..."<<endl;
		if(root->onlyrefine==0) minimizeold(refstart,refend,0); 
		else  			minimizeold(refstart,refend,-10);  	    			 			 		
	}
	segen->mutate=0; 
	if(out&&code) {
		char line[1000];
		sprintf(line,"%s_ini_min.%i.pdb",code,refineid++);
		cerr<<endl<<"minimized initial model written to: "<<line<<endl;
		FILE *fpp=fopen(line,"w");
                fprintf(fpp,"!minimized initial model\n");
		char c=mpdb->chn->id;
        	StrFmt *parent=owner->getparentStrFmt();
        	StrFmt *se = parent->getsequenceStrFmt();
        	mpdb->chn->id=se->cid;
        	if(mpdb->chn->id=='0') mpdb->chn->id=' ';

                mpdb->write(fpp);
                fclose(fpp);
		mpdb->chn->id=c;
	}	
}

void Mutate::setconsistentdssp(){

	Res *r;
	for(r=mpdb->chn->res;r;r=r->next) {
		int n=compare[r->id0];
		int m=owner->match[n];
		if(m==-1) {//insert
			r->sec='-';
		}
		if(dssp[n]!=r->sec) r->sec='-';
		if(TRES.logg>3) cerr<<"consistent: "<<r->id0<<" "<<r->sec<<endl;
	}
}

void Mutate::setownerdssp(){

	Res *r;
	for(r=mpdb->chn->res;r;r=r->next) {
		int n=compare[r->id0];
		r->sec=dssp[n];
		if(TRES.logg>3) cerr<<"consistent: "<<r->id0<<" "<<r->sec<<endl;
	}
}

void Mutate::refinement() {

	Res *refstart=mpdb->chn->res;
	int  refend=100000; 
 
	StrFmt *parent=owner->getparentStrFmt();
	StrFmt *root=owner->getrootStrFmt();

	cerr<<endl;
	cerr<<"begin structure refinement..."<<endl;
	cerr<<"the secondary structure regions are defined as consensus between template and the built model..."<<endl;
	cerr<<endl;
	if(root->refstart>-90000&&root->refend<90000)  {	
		refstart=mpdb->chn->isresoid(root->refstart);
		if(refstart==0) {
			cerr<<endl<<"no residue found with id: "<<root->refstart<<endl<<endl;
			exit(0);
		}
		refend=min(root->refend,mpdb->chn->lastres()->oid);
		Res *r=mpdb->chn->findsmallres(refend);
		if(r) refend=r->id0;
		else {
			cerr<<endl<<"no residue found with id: "<<root->refend<<endl<<endl;
			exit(0);
		} 
	}

	//TRES.smoothclash=0;
	if(parent->iscom==0) {
		restraint=1;
	}
	else if(parent->iscom==1) {
		restraint=2;
	}
	else {
		restraint=1;
	}
	char line[100];
 	
	//save the original
	updatepdbold();
	
	//set dssp
	if(parent->iscom==0) {
		sethbonddssp();
	}
	else {
		sethbonddssp(0,4);
	}
 	setconsistentdssp();
	//refine loop regions, 2, or 3;
	if(sharp>=2) {
		cerr<<"minimizing unconserved sidechains..."<<endl;
		minbadsidechain(refstart,refend,0);
	}
	if(sharp>=2) {
		if(sharp==2) {
			cerr<<endl<<"refinement loop regions with insertion and deletions.."<<endl<<endl;
		}
		else {
			cerr<<endl<<"refinement all loop regions.."<<endl<<endl;
		}
		//refineloop(mpdb->chn->res,100000);
		refineloop(refstart,refend);
		segen->mutate=0; 
		//refineminloop(mpdb->chn->res,100000);
		refineminloop(refstart,refend);
		segen->mutate=0; 
		if(code&&out) {
			char line[1000];
			sprintf(line,"%s_loop_refine.%i.pdb",code,refineid++);
			cerr<<endl<<"output structure after loop refinement..."<<line<<endl<<endl;
			FILE *fpp=fopen(line,"w");
			fprintf(fpp,"!model after loop refinement\n");
			char c=mpdb->chn->id;
        		StrFmt *parent=owner->getparentStrFmt();
        		StrFmt *se = parent->getsequenceStrFmt();
        		mpdb->chn->id=se->cid;
        		if(mpdb->chn->id=='0') mpdb->chn->id=' ';
			mpdb->write(fpp);
                	fclose(fpp);
			mpdb->chn->id=c;
		}
	}
	
	if(parent->iscom==0) {
		sethbonddssp();
	}
	else {
		sethbonddssp(0,4);
	} 
	//sethbonddssp();
	setconsistentdssp(); 
	//refine secondary regions
	if(sharp>=4) {
		cerr<<endl<<"refinement secondary structure regions.."<<endl<<endl;
		//refinesec(mpdb->chn->res,100000);
		refinesec(refstart,refend);
		segen->mutate=0;
		//refineminsec(mpdb->chn->res,100000);
		refineminsec(refstart,refend);
		segen->mutate=0; 
		if(code&&out) {
                        char line[1000];
                        sprintf(line,"%s_sec_refine.%i.pdb",code,refineid++);
                        cerr<<endl<<"output structure after secondary structure refinement..."<<line<<endl<<endl;
                        FILE *fpp=fopen(line,"w");
                        fprintf(fpp,"!model after secondary structure refinement\n");
			char c=mpdb->chn->id;
        		StrFmt *parent=owner->getparentStrFmt();
        		StrFmt *se = parent->getsequenceStrFmt();
        		mpdb->chn->id=se->cid;
        		if(mpdb->chn->id=='0') mpdb->chn->id=' ';
                        mpdb->write(fpp);
                        fclose(fpp);
			mpdb->chn->id=c;
                }
	}
 
	cerr<<endl<<"refinement sidechains with scap..."<<endl;
	cerr<<"with original conformation set as initial conformations for sidechain refinement..."<<endl<<endl;
	
	//minsidechain(mpdb->chn->res,100000);
	minsidechain(refstart,refend);
	if(code&&out) {	
		sprintf(line,"%s_side_refine.%i.pdb",code,refineid++);
		cerr<<endl<<"output structure after sidechain refinement..."<<line<<endl<<endl;
                FILE *fpp=fopen(line,"w");
                fprintf(fpp,"!model after sidechain structure refinement\n");
		char c=mpdb->chn->id;
        	StrFmt *parent=owner->getparentStrFmt();
        	StrFmt *se = parent->getsequenceStrFmt();
        	mpdb->chn->id=se->cid;
        	if(mpdb->chn->id=='0') mpdb->chn->id=' ';
                mpdb->write(fpp);
                fclose(fpp);
		mpdb->chn->id=c;
	} 
}

void Mutate::refineminment() {

	Res *refstart=mpdb->chn->res;
	int  refend=100000; 
 
	StrFmt *parent=owner->getparentStrFmt();
	StrFmt *root=owner->getrootStrFmt();

	cerr<<endl;
	cerr<<"begin structure minimization..."<<endl;
	cerr<<"the secondary structure regions are defined as consensus between template and the built model..."<<endl;
	cerr<<endl;
	if(root->refstart>-90000&&root->refend<90000)  {	
		refstart=mpdb->chn->isresoid(root->refstart);
		if(refstart==0) {
			cerr<<endl<<"no residue found with id: "<<root->refstart<<endl<<endl;
			exit(0);
		}
		refend=min(root->refend,mpdb->chn->lastres()->oid);
		Res *r=mpdb->chn->findsmallres(refend);
		if(r) refend=r->id0;
		else {
			cerr<<endl<<"no residue found with id: "<<root->refend<<endl<<endl;
			exit(0);
		} 
	}

	//TRES.smoothclash=0;
	if(parent->iscom==0) {
		restraint=1;
	}
	else if(parent->iscom==1) {
		restraint=2;
	}
	else {
		restraint=1;
	}
	char line[100];
 	
	//save the original
	updatepdbold();
	
	//set dssp
	if(parent->iscom==0) {
		sethbonddssp();
	}
	else {
		sethbonddssp(0,4);
	}
 	setconsistentdssp();
	//refine loop regions, 2, or 3;
	
	if(sharp>=2) {
		if(sharp==2) {
			cerr<<endl<<"minimize loop regions with insertion and deletions.."<<endl<<endl;
		}
		else {
			cerr<<endl<<"minimize all loop regions.."<<endl<<endl;
		}
		//refineloop(mpdb->chn->res,100000);
		//refineloop(refstart,refend);
		segen->mutate=0; 
		//refineminloop(mpdb->chn->res,100000);
		refineminloop(refstart,refend);
		segen->mutate=0; 
		if(code&&out) {
			char line[1000];
			sprintf(line,"%s_loop_min.%i.pdb",code,refineid++);
			cerr<<endl<<"output structure after loop minimization..."<<line<<endl<<endl;
			FILE *fpp=fopen(line,"w");
			fprintf(fpp,"!model after loop minimization\n");
			char c=mpdb->chn->id;
        		StrFmt *parent=owner->getparentStrFmt();
        		StrFmt *se = parent->getsequenceStrFmt();
        		mpdb->chn->id=se->cid;
        		if(mpdb->chn->id=='0') mpdb->chn->id=' ';
			mpdb->write(fpp);
                	fclose(fpp);
			mpdb->chn->id=c;
		}
	}
	
	if(parent->iscom==0) {
		sethbonddssp();
	}
	else {
		sethbonddssp(0,4);
	} 
	//sethbonddssp();
	setconsistentdssp(); 
	//refine secondary regions
	if(sharp>=4) {
		cerr<<endl<<"minimize secondary structure regions.."<<endl<<endl;
		//refinesec(mpdb->chn->res,100000);
		//refinesec(refstart,refend);
		segen->mutate=0;
		//refineminsec(mpdb->chn->res,100000);
		refineminsec(refstart,refend);
		segen->mutate=0; 
		if(code&&out) {
                        char line[1000];
                        sprintf(line,"%s_sec_min.%i.pdb",code,refineid++);
                        cerr<<endl<<"output structure after secondary structure minimization..."<<line<<endl<<endl;
                        FILE *fpp=fopen(line,"w");
                        fprintf(fpp,"!model after secondary structure minimization\n");
			char c=mpdb->chn->id;
        		StrFmt *parent=owner->getparentStrFmt();
        		StrFmt *se = parent->getsequenceStrFmt();
        		mpdb->chn->id=se->cid;
        		if(mpdb->chn->id=='0') mpdb->chn->id=' ';
                        mpdb->write(fpp);
                        fclose(fpp);
			mpdb->chn->id=c;
                }
	}
 	/*
	cerr<<endl<<"refinement sidechains with scap..."<<endl;
	cerr<<"with original conformation set as initial conformations for sidechain refinement..."<<endl<<endl;
	*/
	//minsidechain(mpdb->chn->res,100000);
	//minsidechain(refstart,refend);

	if(code&&out&&0) {	
		sprintf(line,"%s_side_refine.%i.pdb",code,refineid++);
		cerr<<endl<<"output structure after sidechain minimization..."<<line<<endl<<endl;
                FILE *fpp=fopen(line,"w");
                fprintf(fpp,"!model after sidechain structure minimization\n");
		char c=mpdb->chn->id;
        	StrFmt *parent=owner->getparentStrFmt();
        	StrFmt *se = parent->getsequenceStrFmt();
        	mpdb->chn->id=se->cid;
        	if(mpdb->chn->id=='0') mpdb->chn->id=' ';
                mpdb->write(fpp);
                fclose(fpp);
		mpdb->chn->id=c;
	} 
}


void Mutate::smoothmin() {

	Res *refstart=mpdb->chn->res;
	int  refend=100000; 
 
	StrFmt *parent=owner->getparentStrFmt();
	StrFmt *root=owner->getrootStrFmt();

	if(root->refstart>-90000&&root->refend<90000)  {	
		refstart=mpdb->chn->isresoid(root->refstart);
		if(refstart==0) {
			cerr<<endl<<"no residue found with id: "<<root->refstart<<endl<<endl;
			exit(0);
		}
		refend=min(root->refend,mpdb->chn->lastres()->oid);
		Res *r=mpdb->chn->findsmallres(refend);
		if(r) refend=r->id0;
		else {
			cerr<<endl<<"no residue found with id: "<<root->refend<<endl<<endl;
			exit(0);
		} 
	}


	//StrFmt *parent=owner->getparentStrFmt();
	if(parent->iscom==0) {
		restraint=1;
	}
	else if(parent->iscom==1) {
		restraint=2;
	}
	else {
		restraint=1;
	}
 	
	//save the original
	updatepdbold();
	
	//set dssp
	
	//refine loop regions, 2, or 3;
	int smt=0;
	int smtot=10;
	 
	if(fapr==0) smtot=10;
	else if(fapr==1) smtot=5;
	else if(fapr==2) smtot=0;
		 
	for(smt=smtot;smt>=0;smt--) {
		 
		TRES.smoothclash=1;
		 
		//set coefficient
		float a,b;
		if(smtot==10) {
			a=smt*0.1;
			b=(smtot-smt)*0.9/10+0.1; 
		}
		else if(smtot==5) {
			a=smt*0.2;
			b=(smtot-smt)*0.8/5+0.2;
		}
		else {
			TRES.smoothclash=0;
		}
		TRES.setnewengcoeff(a,b); 
		 

		//set hbond list
		sethbonddssp();
 		setconsistentdssp();
		if(sharp>=2) {
			if(sharp==2) {
				cerr<<endl<<"refinement loop regions with insertion and deletions.."<<endl<<endl;
			}
			else {
				cerr<<endl<<"refinement all loop regions.."<<endl<<endl;
			}
			//refineloop(mpdb->chn->res,100000);
			//segen->mutate=0; 
			//refinesmoothloop(mpdb->chn->res,100000);
			refinesmoothloop(refstart,refend);
			segen->mutate=0; 
			if(code&&out) {
				char line[1000];
				sprintf(line,"%s_loop_refine.%i.pdb",code,refineid++);
				cerr<<endl<<"output structure after loop refinement..."<<line<<endl<<endl;
				FILE *fpp=fopen(line,"w");
				fprintf(fpp,"!model after loop refinement\n");
				char c=mpdb->chn->id;
        			StrFmt *parent=owner->getparentStrFmt();
        			StrFmt *se = parent->getsequenceStrFmt();
        			mpdb->chn->id=se->cid;
        			if(mpdb->chn->id=='0') mpdb->chn->id=' ';
				mpdb->write(fpp);
                		fclose(fpp);
				mpdb->chn->id=c;
			}
		}
	
	 	//set hbond list
		sethbonddssp();
		setconsistentdssp(); 
		//refine secondary regions
 
		if(sharp>=4) {
			cerr<<endl<<"refinement secondary structure regions.."<<endl<<endl;
			//refinesec(mpdb->chn->res,100000);
			//segen->mutate=0;
			//refinesmoothsec(mpdb->chn->res,100000);
			refinesmoothsec(refstart,refend);
			segen->mutate=0; 
			if(code&&out) {
                        	char line[1000];
                        	sprintf(line,"%s_sec_refine.%i.pdb",code,refineid++);
                        	cerr<<endl<<"output structure after secondary structure refinement..."<<line<<endl<<endl;
                        	FILE *fpp=fopen(line,"w");
                        	fprintf(fpp,"!model after secondary structure refinement\n");
				char c=mpdb->chn->id;
        			StrFmt *parent=owner->getparentStrFmt();
        			StrFmt *se = parent->getsequenceStrFmt();
        			mpdb->chn->id=se->cid;
        			if(mpdb->chn->id=='0') mpdb->chn->id=' ';
                       	 	mpdb->write(fpp);
                        	fclose(fpp);
				mpdb->chn->id=c;
               	 	}
		}
	}
	
 	/*
	cerr<<endl<<"refinement sidechains with scap..."<<endl;
	cerr<<"with original conformation set as initial conformations for sidechain refinement..."<<endl<<endl;

	minsidechain(mpdb->chn->res,100000);
	
	if(code&&out) {	
		sprintf(line,"%s_side_refine.%i.pdb",code,refineid++);
		cerr<<endl<<"output structure after sidechain refinement..."<<line<<endl<<endl;
                FILE *fpp=fopen(line,"w");
                fprintf(fpp,"!model after sidechain structure refinement\n");
		
                mpdb->write(fpp);
                fclose(fpp);
	} 
	*/
}

void Mutate::writefinalout() {

	char line[100];
	sprintf(line,"%s_final.pdb",code);
	cerr<<endl<<"output the final model: "<<line<<endl<<endl;
	char c=mpdb->chn->id;
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
	mpdb->chn->id=se->cid;
	if(mpdb->chn->id=='0') mpdb->chn->id=' ';
	int ne=mpdb->chn->start;
	mpdb->chn->start=se->start;
	mpdb->write(line);
	mpdb->chn->start=ne;
	mpdb->chn->id=c;
}

void Mutate::beforewritefinalout() {
	cerr<<"======================"<<endl;
	char line[100];
	sprintf(line,"%s_final.pdb",code);
	 
	cerr<<"check out the final output..."<<line<<endl;	
	cerr<<endl;
	Res **tt=checkmodelburiedpolars(1,1);	
	cerr<<endl;
	
	if(tt) {
		cerr<<"using scap on buried or unconserved salt-bridges..."<<endl;
		Scap scprd;
		
		strcpy(scprd.force,"124");
		TRES.setdonaldcharge();
		scprd.dielectric=10;
		scprd.includeself=1;
        	scprd.singletorsion=1;
        	scprd.colonyline=1;
       	 	scprd.colony=2;
        	scprd.ncolony=1;
        	scprd.nncolony=1;
        	scprd.nummore=100;
        	scprd.bmax=2;
        	scprd.tormax=2;
        	scprd.ring=1;
        	scprd.pdb=mpdb;
		scprd.resultout=0;
		setrotamerorder(0);
		mpdb->setflgr(-99999);
		int n=0;
		while(tt[n]) {
			tt[n]->flag=10000;
			n++;
		} 
		scprd.scpred(1);
		TRES.switchcharge(1);
		setrotamerorder(1);
		delete [] tt;tt=0;
	}
	
}

void Mutate::actmutateloop(int dloop) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
	segen->flex=1;
	segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
	segen->randcoil=0;
        //TRES.smoothclash=1;
	segen->smoothclash=1; 
	segen->part=0;
	//segen->charge=1;
 
	 
	SegBed segbed;
	segbed.next=new SegBed();
	segbed.setsize(1000);
	segbed.next->setsize(1000);
        //refine: 0
   	 
	refine=0;
	FILE *fp=0;
	if(TRES.logg)fp=fopen("malgn","w");
	if(TRES.logg)mpdb->write("start.pdb");
        StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

        int nlen=strlen(sqnto);
        mpdb->setlastresorder();
	
 	int m=0;
	int *stem=0;
	int n=actsidechainmutation(100);
	//minimize(mpdb->chn->res,100000,-100000); 
	//mpdb->chn->setbackbonetorsion(mpdb->chn->res,10000);
	
	int insert=0;//int newres=0;
	int isbrk=0; int midr=0;
	int rstem[2];
	rstem[0]=rstem[1]=-1;
	while(1) {
		Res **tmp;	
		segen->secd='-';			
		m=findbadsite(); //find the best insertion and deletion region,most conserved
		if(m==-1) break;	
		midr=m;	
		int tlen=getlength(m);	tlen=2;	
		isbrk=0; 	 
		if(isbreak(m)) { tlen=1000; isbrk=1;}
		//
		if(isbrk) findrealstem(m,rstem);
		//
		tmp=setinitialstructure(m,tlen);	
		tmp=reordertmp(tmp);			 
        	mpdb->setlastresorder();
		linkstructureviaboth(tmp);
		//updatepdbcopy();		
		setreliable(tmp);				
		insert=calcinsert(); 
		if(TRES.logg)mpdb->chn->write("s");  
 		if(tmp==0)   stem=getsegmentworkon(m);		
		else 	     stem=getsegmentworkon(tmp,isbrk);

		//missing residues
		int ntmp=0;while(tmp&&tmp[ntmp])ntmp++;
		int it=0;
		if(rstem[0]==-1||rstem[1]==-1) ntmp=0;
		for(it=0;it<ntmp;it++) {
			int iu=compare[tmp[it]->id0];
			if(iu<rstem[0]||iu>rstem[1]) {
				rstem[0]=rstem[1]=-1;
			}
		}
		if(ntmp==0) rstem[0]=rstem[1]=-1;
		//

		segbed.next->num=0;
		segbed.num=0;
		setsegbed(segbed.next,stem);	 
		if(tmp) delete [] tmp;tmp=0;
		int done=1;
 		
		float *xyzorg=mpdb->chn->gettransfer(mpdb->chn->res,1000000); 
		while(1) {		
			//stem=getnewstem(stem);
			int m1,m2;
			m1=match[stem[0]];
			m2=match[stem[1]];
			if(segbed.ifexist(m1,m2)) {			
				//done++;
				if(done==segbed.next->num) goto redone; 
				stem=findsegbed(segbed.next,done);
				done++;
				continue;
			}
			if(m1==-1||m2==-1) break;
			Res *r1;r1=resn[match[stem[0]]];
			Res *r2;r2=resn[match[stem[1]]];
		 	if(r1&&r2&&insert==0&&mpdb->chn->isalllinked(r1,r2->id0)==1) break;  		
			
			//newres=calcnewres(r1,r2->id0); 			

			if(segen->chiangle) {
                        	delete segen->chiangle;
                        	segen->chiangle=0;
                	}
                	else if(segen->bound) {
                        	delete segen->bound;
                        	segen->bound=0;
                	}						
			//setchiangle(stem);
			
			if(TRES.logg) printalign(fp,stem);
			int nn;nn=getreslength(stem); 
			if(TRES.logg) printsegment(stem);
			delsidechain(stem);
 	 	 	 		
			if(insert>=3) segen->arbt=20;
                        else if(insert>=2) segen->arbt=10;
                        else if(insert>=1) segen->arbt=5;			
			else segen->arbt=5;
			
                        if(fapr==1) {
				if(insert>=3) segen->arbt=5;
                        	else if(insert>=2) segen->arbt=1;
                        	else if(insert>=1) segen->arbt=1;			 
				else segen->arbt=1; 
			}
			//
 			if(rstem[0]!=-1&&rstem[1]!=-1) {				
				if(stem[0]<rstem[0]||stem[1]>rstem[1]) {
					if(done==segbed.next->num) goto redone; 
					stem=findsegbed(segbed.next,done);
					done++;
					continue;		
				}		
			}
			// 
			float *xyzout;xyzout=segen->myfixsegment(stem);
			if(TRES.logg)mpdb->write("s");
			mpdb->chn->transfer(xyzorg,mpdb->chn->res,1000000,1);
			 
			if(stem) {delete [] stem; stem=0;}
			if(xyzout) {
				segbed.add(xyzout,segen->start,segen->end);
				xyzout=0;
				//done++;
				if(done==segbed.next->num) goto redone; 
				stem=findsegbed(segbed.next,done);done++;
				continue;
			}
			else {
				stem=findloosesegment(midr,r1,r2);
				continue;
			}
			redone:
			if(done==segbed.next->num) {	
			 	putsegbed(&segbed);						 
				//actsidechainmutation(100);
				//minimize(mpdb->chn->res,100000,-100000); 				
				break;	
			}						 			 
		}
		if(xyzorg) delete [] xyzorg;xyzorg=0; 
	}	
	
        actsidechainmutation(-99999);
	sethbonddssp();
	updatepdbold();
	reindex(pdbold);
	if(se->source!='T'&&fapr!=1) {
		minimize(mpdb->chn->res,100000,-100000);       
		//char line[100];
		//sprintf(line,"%s.pdb",code);			 
		reindex();	
		//mpdb->write(line);
	}
	else {
		reindex();
		//char line[100];
		//sprintf(line,"%s.pdb",code);	
		//mpdb->write(line);
	}
	
	sethbonddssp();

	StrFmt *composite=owner->getStrFmt("composite");

	if(sharp>0&&code&&composite==0) {
		//minloop(mpdb->chn->res,100000);		
		minwindowloop(mpdb->chn->res,100000);
		char line[100];
		if(code&&out) {
			sprintf(line,"%s_loop.pdb",code);
			char c=mpdb->chn->id;
        		StrFmt *parent=owner->getparentStrFmt();
        		StrFmt *se = parent->getsequenceStrFmt();
        		mpdb->chn->id=se->cid;
        		if(mpdb->chn->id=='0') mpdb->chn->id=' ';
			mpdb->write(line);
			mpdb->chn->id=c;
		}
	}

	if(sharp&&code&&composite==0) {		 
		minwinrefine0(mpdb->chn->res,100000);
		char line[100];
		if(code&&out) {
			sprintf(line,"%s_refine.pdb",code);
			char c=mpdb->chn->id;
        		StrFmt *parent=owner->getparentStrFmt();
        		StrFmt *se = parent->getsequenceStrFmt();
        		mpdb->chn->id=se->cid;
        		if(mpdb->chn->id=='0') mpdb->chn->id=' ';
			mpdb->write(line);
			mpdb->chn->id=c;
		}
	}
	
	if(sharp==4&&code&&composite==0) {
		minsidechain(mpdb->chn->res,100000);
		char line[100];
		if(code&&out) {
			sprintf(line,"%s_side.pdb",code);
			char c=mpdb->chn->id;
        		StrFmt *parent=owner->getparentStrFmt();
        		StrFmt *se = parent->getsequenceStrFmt();
        		mpdb->chn->id=se->cid;
        		if(mpdb->chn->id=='0') mpdb->chn->id=' ';
			mpdb->write(line);
			mpdb->chn->id=c;
		}
	}
 
	if(code&&composite==0) {
		char line[100];
		sprintf(line,"%s.pdb",code);
		char c=mpdb->chn->id;
        	StrFmt *parent=owner->getparentStrFmt();
        	StrFmt *se = parent->getsequenceStrFmt();
        	mpdb->chn->id=se->cid;
        	if(mpdb->chn->id=='0') mpdb->chn->id=' ';
		mpdb->write(line);
		 
		mpdb->chn->id=c;
	}	

	else if(code&&composite&&out) {
		char line[100];
		sprintf(line,"%s_org.pdb",code);
		char c=mpdb->chn->id;
        	StrFmt *parent=owner->getparentStrFmt();
        	StrFmt *se = parent->getsequenceStrFmt();
        	mpdb->chn->id=se->cid;
        	if(mpdb->chn->id=='0') mpdb->chn->id=' ';
		mpdb->write(line);
		mpdb->chn->id=c;
	}
}
void Mutate::minshift(Res *rr,int nn) {

		
	minwinrefine0(mpdb->chn->res,100000);
}

void Mutate::sethbonddssp() {
	mpdb->chn->clearhbond();
	mpdb->chn->header();
	mpdb->chn->buildhbond();
	mpdb->chn->setdsspstr();
	mpdb->chn->buildssbond();
	mpdb->chn->setthreestatesec();
	if(TRES.logg>3) mpdb->chn->writesecondary(stderr);
	//setdssp(1);
}

void Mutate::sethbonddssp(int from0,int end0) {
	mpdb->chn->clearhbond();
	mpdb->chn->header();
	mpdb->chn->buildhbond(from0,end0);
	mpdb->chn->setdsspstr();
	mpdb->chn->buildssbond(from0,end0);
	mpdb->chn->setthreestatesec();
	if(TRES.logg>3) mpdb->chn->writesecondary(stderr);
	//setdssp(1);
}

void Mutate::buildhbond() {
	mpdb->chn->clearhbond();
	mpdb->chn->header();
	mpdb->chn->buildhbond();
	mpdb->chn->buildssbond();
}

void Mutate::sethbonddssp(int n) {
	mpdb->chn->clearhbond();
	mpdb->chn->header();
	mpdb->chn->buildhbond();
	mpdb->chn->setdsspstr(n);
	mpdb->chn->buildssbond();
	mpdb->chn->setthreestatesec();
	mpdb->chn->writesecondary(stderr);
	//setdssp(1);
}
void Mutate::combine() {

	cerr<<endl<<"working on composite model building :"<<code<<endl<<endl;

	if(segen) delete segen;segen=0;
        segen=new Segen();
	segen->flex=1;
	segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
	segen->randcoil=0;
        //TRES.smoothclash=0;
	segen->smoothclash=0; 
	segen->part=0;
 	segen->cutoff=6.;
	//
	restraint=1;
	//	 
	StrFmt *root=owner->getrootStrFmt();
	StrFmt *parent=owner->getparentStrFmt();
	StrFmt *composite=owner->getStrFmt("composite");
	StrFmt *se=owner->getsequenceStrFmt();	
	if(composite==0||se==0) {
		cerr<<"there is no composite or sequence pir in the alignment set:"<<endl;
		parent->printout(stderr);
		return;
	}

	refine=0;
	 
	int n=strlen(composite->seqngap);
	int i;
	setrightids();	
 
	i=0;
	while(1) {
		i++;
		int *stem=composite->getresalncomposite(i);
		
		if(stem==0) break;
		stem=resetstemcomposite(stem);
		composite->start=composite->resaln[i*2-2];
		composite->end=composite->resaln[i*2-1];
		if(stem[0]>stem[1]) continue;
		int i1=stem[0];
		int i2=stem[1];
		if(match[i1]==-1||match[i2]==-1) continue;
		int ne=handlecomposite(stem);			
		delete [] stem;stem=0;	
	} 
	
	if(out>0&&code) {
		char line[100];
		sprintf(line,"%s_composite.%i.pdb",code,refineid++);
		cerr<<endl<<"output the composite structure: "<<line<<endl<<endl;
		char c=mpdb->chn->id;
        	StrFmt *parent=owner->getparentStrFmt();
        	StrFmt *se = parent->getsequenceStrFmt();
        	mpdb->chn->id=se->cid;
        	if(mpdb->chn->id=='0') mpdb->chn->id=' ';
		mpdb->write(line);
		mpdb->chn->id=c;
	}

	segen->mutate=0;
}
void Mutate::combineold() {
	if(segen) delete segen;segen=0;
        segen=new Segen();
	segen->flex=1;
	segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
	segen->randcoil=0;
	//TRES.smoothclash=0;
	segen->smoothclash=0; 
	segen->part=0;
 	segen->cutoff=6.;
	//
	restraint=1;
	//	 
	StrFmt *root=owner->getrootStrFmt();
	StrFmt *parent=owner->getparentStrFmt();
	StrFmt *composite=owner->getStrFmt("composite");
	StrFmt *se=owner->getsequenceStrFmt();	
	if(composite==0||se==0) {
		cerr<<"there is no composite or sequence pir in the alignment set:"<<endl;
		parent->printout(stderr);
		return;
	}

	refine=0;
	 
	int n=strlen(composite->seqngap);
	int i;
	setrightids();	
 
	 
	char c=se->zipcode;
	for(i=0;i<n;i++) {
		if(composite->seqngap[i]==c) continue;
		if(se->seqngap[i]=='-') continue;
		if(composite->seqngap[i]!=se->zipcode) {
			int *stem=getcompositestem(i);
			if(stem==0) goto re200;
			int ne=handlecomposite(stem);	
			delete [] stem;stem=0;	
		}	
		re200:
		c=composite->seqngap[i];
	}
	 
	
	if(out>0&&code) {
		char line[100];
		sprintf(line,"%s_refine_ini.pdb",code);
		char c=mpdb->chn->id;
        	StrFmt *parent=owner->getparentStrFmt();
        	StrFmt *se = parent->getsequenceStrFmt();
        	mpdb->chn->id=se->cid;
        	if(mpdb->chn->id=='0') mpdb->chn->id=' ';
		mpdb->write(line);
		mpdb->chn->id=c;
	}

	segen->mutate=0;
}
Res *Mutate::findresid(int n) {

	int nlen=strlen(seqn);
	
	int i;
 
	for(i=0;i<nlen;i++) {
		if(resid[i]==-1) continue;
		if(resid[i]!=n) continue;
		return resn[i]; 
	}
	return 0;
}

int Mutate::findresid(Res *t) {

	int nlen=strlen(seqn);
	
	int i;
 
	for(i=0;i<nlen;i++) {
		if(resn[i]==t) {
			return resid[i];
		}		 
	}
	return -1;
}



void Mutate::setrightids() {
	Strhandler cc;
	StrFmt *root=owner->getrootStrFmt();
	StrFmt *t;
	
	int nlen=strlen(seqn);
	
	int i;
	if(resid) delete [] resid;resid=0;
	resid=new int[nlen];
 	for(i=0;i<nlen;i++) resid[i]=-1;
	for(i=0;i<nlen;i++) {
		Res *r=resn[i];
		if(r) resid[i]=r->id0;
	}

	StrFmt *sqA=owner->getStrFmt("sequence");
	
	for(t=root;t;t=t->next) {

		if(t->mutate==0) continue;
		if(t->mutate==this) continue;
		if(t->iscom==1)continue;
 		if(isfather(t->mutate)==1) continue;

		int tlen=strlen(t->mutate->seqn);
		if(t->mutate->resid) delete [] t->mutate->resid;t->mutate->resid=0;
		t->mutate->resid=new int[tlen];

		for(i=0;i<tlen;i++) {		
			t->mutate->resid[i]=-1;
		}
		StrFmt *sqB=t->getStrFmt("sequence");
		cerr<<endl;
		cerr<<"try to find corresponding residues between two models using needle-wunsch alignment..."<<endl;
		cerr<<"the sequence of the composite segment in one model should be aligned with the corresponding segment in the other model"<<endl;
		cerr<<"otherwise, there could be problems"<<endl;
		if(sqA&&sqB&&sqA->zipcode!=' '&&sqA->zipcode!='\0'&&sqB->zipcode!=' '&&sqB->zipcode!='\0') {
			cerr<<"aligning sequence from model: "<<sqA->zipcode<<"----"<<sqB->zipcode<<endl;  
		}
		cerr<<endl;
		Algn algn;
		algn.setsequence(seqn,t->mutate->seqn);		
		algn.defalgn();
		int slen=algn.getroutelength();
		char **result=algn.output(stderr);
		int n1=-1;
		int n2=-1;
		int nn=0;
		for(i=0;i<slen;i++) {
			if(result[0][i]!='-'&&result[0][i]!='X')  n1++;
			if(result[1][i]!='-'&&result[1][i]!='X')  n2++;
			//if(result[1][i]!='-'&&result[1][i]!='X'&&result[0][i]!='-'&&result[0][i]!='X') n2++;
			//if(result[1][i]!='-'&&result[1][i]!='X'&&result[1][i]==result[0][i]) n2++;
			if(result[0][i]=='-'||result[0][i]=='X') continue;
			if(result[1][i]=='-'||result[1][i]=='X') continue;
			if(result[1][i]!=result[0][i]) continue;			
			t->mutate->resid[n2]=resid[n1];			 						 
			nn++;
		}
		if(nn<strlen(seqn)*0.5||nn<strlen(t->mutate->seqn)*0.5) {
			cerr<<endl;
			cerr<<"*******Warning!*******"<<endl;
			cerr<<"the two models do not align well."<<endl;
			cerr<<endl; 
			cerr<<endl;
		}
		cc.strdel(result);
	}
}



int *Mutate::getseqnstem(Mutate *c,Res *r1,Res *r2) {
		
	//from c->seqn, find the sequence start and end number from r1,r2;

	Res *t1=r1;
	Res *t2=r2;
	int n=0;
	while(t2) {	
		char *s=mpdb->chn->getseqn(t1,t2->id0);
		if(s==0)  {
			t2=t2->next;
			continue;
		}
		n=uniquestrstr(c->seqn,s);
		if(s) delete [] s;s=0;
		if(n==1) break;
		t2=t2->next;
	}
	
	if(n==1) goto re200;

	t1=r1;
	t2=r2;
	n=0;
	while(t1) {	
		char *s=mpdb->chn->getseqn(t1,t2->id0);
		if(s==0)  {
			t2=t2->next;
			continue;
		}
		n=uniquestrstr(c->seqn,s);
		if(s) delete [] s;s=0;
		if(n==1) break;
		t1=t1->last;
	}

	if(n==1) goto re200;

	t1=r1;
	t2=r2;
	n=0;
	while(t1) {	
		char *s=mpdb->chn->getseqn(t1,t2->id0);
		if(s==0)  {
			t2=t2->next;
			continue;
		}
		n=uniquestrstr(c->seqn,s);
		if(s) delete [] s;s=0;
		if(n==1) break;
		t1=t1->last;
	}
	
 	if(n!=1)  {
		char *s=mpdb->chn->getseqn(r1,r2->id0);
		cerr<<"\n\n\nthe first sequence could be uniquely defined in the second sequnece: "<<endl;
		cerr<<"\nthe first sequence:"<<endl;
		cerr<<s<<endl;
		cerr<<"\nthe second sequence:"<<endl;
		cerr<<c->seqn<<endl;
		exit(0);
	}
	re200:
	char *s=mpdb->chn->getseqn(t1,t2->id0);
	char *ss=strstr(c->seqn,s);
	int ii=0;
	int nlen=strlen(c->seqn);

	for(ii=0;ii<nlen;ii++) {
		if(c->seqn+ii==ss) {
			break;
		}
	}
	
	int n1=ii;

	Res *e1=t1;
	while(e1) {
		if(e1==r1) break;
		n1++;
		e1=e1->next;
	}
		
	int n2=n1;

	e1=r1;

	while(e1) {
		if(e1==r2) break;
		n2++;
		e1=e1->next;
	}
	 
	int i1=c->compare[n1];
	int i2=c->compare[n2];
	int *stem=new int[2];
	stem[0]=i1;
	stem[1]=i2;
	return stem;
}
int Mutate::uniquestrstr(char *a,char *b) {
	//if a unique exist in b;
	char *s=strstr(b,a);
	if(s==0) return -1;
	s++;
	char *s1=strstr(s,a);
	if(s==0) return 1;
	else return 0;
}


int *Mutate::findresidstems(int n1,int n2) {

	int nlen=strlen(seqn);

 	int m1=n2+1000000;
	int m2=-1;
	int i;
	int ii=0;
	for(i=0;i<nlen;i++) {
		if(resid[i]==-1) continue;
		if(resid[i]>=n1&&resid[i]<=n2) {
			m1=min(i,m1);
			m2=max(i,m2);
			ii++;
		}
	}

	if(ii==0) {
		return 0;
	}
	//if(m2==-1) return 0;

	int *stem=new int[2];

	stem[0]=compare[m1];
	stem[1]=compare[m2];

	return stem;
}
 
int Mutate::handlecomposite(int *stem){
 
	StrFmt *root=owner->getrootStrFmt();
	StrFmt *composite=owner->getStrFmt("composite");
	StrFmt *se=owner->getStrFmt("sequence");
	StrFmt *c,*t;

	int n1=stem[0];
	int n2=stem[1];
	int n=-1;

	int i;
	for(i=n1;i<=n2;i++) {	
		if(composite->seqngap[i]!=' '&&composite->seqngap[i]!=se->zipcode) {
			n=i;break;	
		}
	}

	if(n==-1) return 0;

	if(composite->seqngap[n]!='-'&&composite->seqngap[n]!='+') {
		if(composite->seqngap[n]==se->zipcode) c=se;
		else c=root->getStrFmtWithId(composite->seqngap[n]);
		if(c==0) {
			cerr<<"\n\n\ncan not find sequence pir with id:"<<composite->seqngap[n]<<endl;
			exit(0);
		}
		 
		t=c->getparentStrFmt();
		if(t->mutate==0||t->mutate->mpdb==0||t->mutate->mpdb->chn==0) {
			cerr<<"\n\n\nthe sequence pir does not have predicted strucutre:"<<endl;
			t->printout(stderr);
			exit(0);
		}
		int i1=match[n1];
		int i2=match[n2];
		Res *r1=resn[i1];		
		Res *r2=resn[i2];
 		cerr<<endl<<"trying to replace segments: "<<r1->name<<r1->id0<<"---"<<r2->name<<r2->id0<<endl;
		cerr<<"with segment from model with id:"<<composite->seqngap[n]<<endl<<endl;
		assemble(t->mutate,r1,r2);		
		root->setupatmid0();
		if(t->mutate->mpdb) t->mutate->mpdb->configureatmid0();
		if(mpdb) mpdb->configureatmid0();
		if(code&&out>1) {
			char nn[100];
                	sprintf(nn,"%s_composite_tmp.%i.pdb",code,refineid++);
                	cerr<<endl<<"output intermediate composite structures: "<<nn<<endl<<endl;
			char c=mpdb->chn->id;
        		StrFmt *parent=owner->getparentStrFmt();
        		StrFmt *se = parent->getsequenceStrFmt();
        		mpdb->chn->id=se->cid;
        		if(mpdb->chn->id=='0') mpdb->chn->id=' ';
                	mpdb->write(nn);    
			mpdb->chn->id=c;
    
		}
	}
	else if(composite->seqngap[n]=='-'){
		int i1=match[n1];
		int i2=match[n2];
		Res *r1=resn[i1];		
		Res *r2=resn[i2];
		cerr<<endl<<"trying to replace segments: "<<r1->name<<r1->id0<<"---"<<r2->name<<r2->id0<<endl;
                cerr<<"with segment from all other models..."<<endl<<endl;
		assembleall(r1,r2);	
		root->setupatmid0();
		if(code&&out>1) {
                        char nn[100];
                        sprintf(nn,"%s_composite_tmp.%i.pdb",code,refineid++);
                        cerr<<endl<<"output intermediate structures of the lowest energy from multiple templates: "<<nn<<endl<<endl;
  			char c=mpdb->chn->id;
        		StrFmt *parent=owner->getparentStrFmt();
        		StrFmt *se = parent->getsequenceStrFmt();
        		mpdb->chn->id=se->cid;
        		if(mpdb->chn->id=='0') mpdb->chn->id=' ';
                        mpdb->write(nn);
			mpdb->chn->id=c;
                }
	}
	else if(composite->seqngap[n]=='+'){
		int i1=match[n1];
		int i2=match[n2];
		Res *r1=resn[i1];		
		Res *r2=resn[i2];
		cerr<<endl<<"trying to average segments: "<<r1->name<<r1->id0<<"---"<<r2->name<<r2->id0<<endl;
                cerr<<"with segment from all other models..."<<endl<<endl;
		assembleallmean(r1,r2);	
		root->setupatmid0();
		if(code&&out>1) {
                        char nn[100];
                        sprintf(nn,"%s_composite_tmp.%i.pdb",code,refineid++);
                        cerr<<endl<<"output intermediate structures of the lowest energy from multiple templates: "<<nn<<endl<<endl;
			char c=mpdb->chn->id;
        		StrFmt *parent=owner->getparentStrFmt();
        		StrFmt *se = parent->getsequenceStrFmt();
        		mpdb->chn->id=se->cid;
        		if(mpdb->chn->id=='0') mpdb->chn->id=' ';
                        mpdb->write(nn);
			mpdb->chn->id=c;
                }
	} 
	return 1;
}
 
void Mutate::algnsegment(int *stem) {

	int nr1,nr2;
        nr1=stem[0];nr2=stem[1];
	int i1=match[nr1];
	int i2=match[nr2];
	Res *r1=resn[i1];		
	Res *r2=resn[i2];
	StrFmt *root=owner->getrootStrFmt();
	StrFmt *composite=owner->getStrFmt("composite");
	StrFmt *t;
	
	for(t=root;t;t=t->next) {
		if(t->iscom==1) continue;		 
		if(t->mutate==0||t->mutate->mpdb==0||t->mutate->mpdb->chn==0) continue;
		if(t->mutate==this) continue;
		if(isfather(t->mutate)==1) continue;
		
		int *stem0=t->mutate->findresidstems(r1->id0,r2->id0);
		int i1=stem0[0];
		int i2=stem0[1];
		delete [] stem0;stem0=0;
		i1=t->mutate->match[i1];
		i2=t->mutate->match[i2];
		Res *t1=t->mutate->resn[i1];
		Res *t2=t->mutate->resn[i2];
		int n1=t->mutate->resid[i1];
		int n2=t->mutate->resid[i2]; 
		Res *e1=mpdb->chn->isres0(n1);
		Res *e2=mpdb->chn->isres0(n2);
		char *s1=e1->chn->getseqn(e1,e2->id0);
		char *s2=t1->chn->getseqn(t1,t2->id0);
		if(s1==0||s2==0||strcmp(s1,s2)!=0) {
			cerr<<"\n\n\ncomposite replacement could not be performed :"<<endl;			
			cerr<<composite->seqngap<<endl;
			cerr<<"\nthe reason is: "<<endl;
			int i;
			for(i=n1;i<=n2;i++) cerr<<composite->seqngap[i];			
			cerr<<endl;
			cerr<<"does not have identical residues in the two subsequencs between";
			cerr<<" original and the replacement: "<<endl;
			cerr<<"original: \n"<<s1<<endl;
			cerr<<"replace: \n"<<s2<<endl;
			exit(0);
		}
		if(s1) delete [] s1;s1=0;
		if(s2) delete [] s2;s2=0;
		//Stralg algn;algn.flag=1;
		//algn.superimposesimple(t1,t2->id0,e1,e2->id0,1);
		if(composite->start||composite->end) {
			segexcludesuperimpose(t1,t2,e1,e2);
			if(out>1&&code) {
                       	 	char nn[100];
                       	 	sprintf(nn,"%s_align.%i.pdb",code,refineid++);
                        	FILE *fp;
                        	fp=fopen(nn,"w");
                        	char a;
                       	 	a=mpdb->chn->id;
                        	mpdb->chn->id='A';
                        	mpdb->chn->write(fp);
                        	mpdb->chn->id=a;
                        	//
                        	a=t->mutate->mpdb->chn->id;
                        	t->mutate->mpdb->chn->id='B';
                        	t->mutate->mpdb->chn->write(fp);
                        	t->mutate->mpdb->chn->id=a;
                        	//
                        	fclose(fp);
                	}
		} 
		if(TRES.logg>3) t1->chn->write("sa");
		if(TRES.logg>3) e1->chn->write("sb");
	}
}

int *Mutate::findrightids(Mutate *t, int *stem) {

	//stem is the gap position

	if(t==0||stem==0) return 0;
	
	int nr1,nr2;
        nr1=stem[0];
	nr2=stem[1];

	//get sequence position
	int i1=match[nr1];
	int i2=match[nr2];

	//get sequecne position
	Res *r1=resn[i1];		
	Res *r2=resn[i2];
	if(stem) delete [] stem;stem=0;

	//find the segment between r1 and r2;

	//get t's gap position
	int *stem0=t->findresidstems(r1->id0,r2->id0);
	if(stem0==0) return 0;
	i1=stem0[0];
	i2=stem0[1];
	
	//go t's sequence position
	i1=t->match[i1];
	i2=t->match[i2];
	 
	//go t's original id
	int n1=t->resid[i1];
	int n2=t->resid[i2];
	if(n1==-1||n2==-1) {
		delete [] stem0;
		stem0=0;
		return 0;
	}

	//get orignal gap position

	stem0[0]=compare[n1];
	stem0[1]=compare[n2];
	 
	return stem0;
}
void Mutate::segexcludesuperimpose(Res *r1,Res *r2,Res *t1,Res *t2) {
	Res **alga,**algb;
	Chn *chna,*chnb; 
  	Res *r,*t;
  	int n,na,nb,window;
 
  	chna=r1->chn;
  	chnb=t1->chn;
  	na=0; for(r=chna->res;r;r=r->next) na++;
  	nb=0; for(r=chnb->res;r;r=r->next) nb++;
  	 
  	n=max(na,nb)+100;
 	alga=new Res*[2*n];
  	algb=new Res*[2*n];
	for(na=0;na<2*n;na++) {alga[na]=0;algb[na]=0;}
	
	na=0; 
 	int start,end;
	StrFmt *sc=owner->getStrFmt("composite");
	start=sc->start;
	end=sc->end;
	if(start>0) {	
		window=start;
		nb=0;
		r=r1->last;t=t1->last;
		while(r&&t&&r->name==t->name) {				 	
			alga[na]=r;
			algb[na]=t;
			r=r->last;t=t->last;
			if(nb>window) break; 
			nb++;na++;
		}
	}
	else if(start<0) {
		window=-start;
		nb=0;
		r=r1->next;t=t1->next;
		while(r&&t&&r->name==t->name) {		
			alga[na]=r;
			algb[na]=t;
			r=r->next;t=t->next;
			if(na>window) break; 
			na++;nb++;			
		}
	}

	if(end>0) {
		window=end;
		nb=0;
		r=r2->next;t=t2->next;
		while(r&&t&&r->name==t->name) {		
			alga[na]=r;
			algb[na]=t;
			r=r->next;t=t->next;
			if(nb>window) break; 
			na++;nb++;
		}
	}
	else if(end<0) {
		window=-end;
		nb=0;
		r=r2->last;t=t2->last;
		while(r&&t&&r->name==t->name) {		
			alga[na]=r;
			algb[na]=t;
			r=r->next;t=t->next;
			if(nb>window) break; 
			na++;nb++;
		}
	}

	int i,j;
	for(i=0;i<na;i++) {
		for(j=i+1;j<na;j++) {
			if(alga[i]==alga[j]&&algb[i]==algb[j]) {
				alga[j]=0;algb[j]=0;
			}
		}
	}

	nb=0;
	for(i=0;i<na;i++) {
		if(alga[i]==0||algb[i]==0) continue;
		alga[nb]=alga[i];
		algb[nb]=algb[i];
		nb++;
	}
	alga[nb]=0;
	algb[nb]=0;	 
	if(nb==0) {
		cerr<<endl<<"there is no corresponding residues specified for aligning  two structures"<<endl<<endl;
		if(alga) delete [] alga;alga=0;
		if(algb) delete [] algb;algb=0;	
		return;
	}
	//write out
	cerr<<endl<<"superimpose the two segments with corresponding residues..."<<endl;
	/*
	for(i=0;i<nb;i++) {
		cerr<<"corresponding residues: "<<alga[i]->name<<alga[i]->id0<<"---"<<algb[i]->name<<algb[i]->id0<<endl;
	}
	*/
	//
	Stralg algn;
	algn.flag=1;
	algn.alga=alga;
	algn.algb=algb;
	float dd=algn.superimpose(1);
	for(i=0;i<nb;i++) {
		float d=alga[i]->directrmsdanyway(algb[i],0,3);
                cerr<<"corresponding residues: "<<alga[i]->name<<alga[i]->id0<<"---"<<algb[i]->name<<algb[i]->id0<<" with rmsd:"<<d<<endl;
        }

	cerr<<endl<<"the above rmsd :"<<dd<<" A"<<endl<<endl;
	if(TRES.logg>3) r1->chn->write("r1");
	if(TRES.logg>3) t1->chn->write("t1");	
	if(alga) delete [] alga;alga=0;
	if(algb) delete [] algb;algb=0;	
	algn.alga=alga;algn.algb=algb;
}


int *Mutate::findnewseg(Res *rr,Res *rr1, Res *rr2,int f) {
 
	int n1=rr->id0-rr1->id0;
	int n2=rr2->id0-rr->id0;
	
	Res *first=mpdb->chn->res;
	Res *last=mpdb->chn->lastres();
	
	 
	if(rr1->last==0&&rr2->next==0) {
		return 0;
	}
	else if(rr1->last==0&&rr2->next) {
		rr2=rr2->next;
	}
	else if(rr1->last&&rr2->next==0) {
		rr1=rr1->last;
	} 
	else if(rr1->id0-first->id0<2&&last->id0-rr2->id0>2) {
		rr2=rr2->next; 
	}
	else if(rr1->id0-first->id0>2&&last->id0-rr2->id0<2) {
		rr1=rr1->last; 
	}
	else if(n2-n1>2&&f==0) {
		rr1=rr1->last; 
	}
	else if(n1-n2>2&&f==1) {
		rr2=rr2->next;
	}
	else if(f==0) {
		rr2=rr2->next;
	}
	else if(f==1) {
		rr1=rr1->last;
	}	 

	int *stem=new int[2];
	stem[0]=compare[rr1->id0];
	stem[1]=compare[rr2->id0];
	return stem;
}

int *Mutate::findnewseg(Res *rr,Res *tt,Res *rr1, Res *rr2,int f) {
 
	//int n1=rr->id0-rr1->id0;
	//int n2=rr2->id0-rr->id0;
		
	Res *first=mpdb->chn->res;
	Res *last=mpdb->chn->lastres();
 
	if(rr1->last==0&&rr2->next==0) {
		return 0;
	}
	else if(rr1->last==0&&rr2->next) {
		rr2=rr2->next;
	}
	else if(rr1->last&&rr2->next==0) {
		rr1=rr1->last;
	}
	else if(rr1->id0-first->id0<2&&last->id0-rr2->id0>2) {
		rr2=rr2->next; 
	}
	else if(rr1->id0-first->id0>2&&last->id0-rr2->id0<2) {
		rr1=rr1->last; 
	}	
	else if(f==0&&rr2->id0<tt->last->id0-1) {
		rr2=rr2->next;
	}
	else if(f==1&&rr1->id0>rr->next->id0+1) {
		rr1=rr1->last;
	}
	else if(f==0) {
		rr1=rr1->last;
	}
	else if(f==1) {
		rr2=rr2->next;
	}	 

	int *stem=new int[2];
	stem[0]=compare[rr1->id0];
	stem[1]=compare[rr2->id0];
	return stem;
}


int Mutate::assemble(Mutate *t,Res *r1,Res *r2) {
	 	
        //set model
	mold=t;
        //
	
	mpdb->setlastresorder();
	mold->mpdb->setlastresorder();
	
	refine=10; //when refine>=10: composite refinement

	//find composite
 	StrFmt *composite=owner->getStrFmt("composite");
 	
	//find ids
	int n1=compare[r1->id0];
	int n2=compare[r2->id0];

	//find corresponding ids in model
	int *stem0=t->findresidstems(r1->id0,r2->id0);  
	
	//
	if(stem0==0) {
		cerr<<"\n\ncomposite replacement of the sequence could not be performed:";
		char *s1=mpdb->chn->getseqn(r1,r2->id0);
		if(s1) {
			cerr<<s1<<endl;
			delete [] s1;s1=0;
		}
		cerr<<endl;
		cerr<<"the reason is no corresponding segment found in the other templates"<<endl;			
		return 0;
	}

	
	int i1=stem0[0];
	int i2=stem0[1];
	delete [] stem0;stem0=0;
	
	//	
	i1=t->match[i1];
	i2=t->match[i2]; 
	Res *t1=t->resn[i1];
	Res *t2=t->resn[i2];	
	char *s1=mpdb->chn->getseqn(r1,r2->id0);
	char *s2=t->mpdb->chn->getseqn(t1,t2->id0);

	//the sequence of the replacement must be equal
	if(s1==0||s2==0||strcmp(s1,s2)!=0) {		
		cerr<<"\n\n\ncomposite replacement could not be performed :"<<endl;			
		cerr<<composite->seqngap<<endl;
		cerr<<"\nthe reason is: "<<endl;
		int i;
		for(i=n1;i<=n2;i++) cerr<<composite->seqngap[i];			
		cerr<<endl;
		cerr<<"does not have identical residues in the two subsequencs between";
		cerr<<" original and the replacement: "<<endl;
		cerr<<"original: \n"<<s1<<endl;
		cerr<<"replace: \n"<<s2<<endl;
		if(s1) delete [] s1;s1=0;	
		if(s2) delete [] s2;s2=0;
		return 0;
		//exit(0);
	}
 	 
	if((strlen(s1)<3&&r1->last&&r2->next)||strlen(s1)==1) {
		cerr<<"the subsequence to be replaced is too few. ignored"<<endl;
		cerr<<"original: \n"<<s1<<endl;
		cerr<<"replace: \n"<<s2<<endl;
		if(s1) delete [] s1;s1=0;	
		if(s2) delete [] s2;s2=0;
		return 0;
	}
	 	
	if(s1) delete [] s1;s1=0;	
	if(s2) delete [] s2;s2=0;
	
	//superimpose the segments or ends
	if(composite->start||composite->end) {
		segexcludesuperimpose(t1,t2,r1,r2);
		if(out>1&&code) {
			char nn[100];
			sprintf(nn,"%s_align.%i.pdb",code,refineid++);
			FILE *fp;
			fp=fopen(nn,"w");
			char a;
			a=mpdb->chn->id;
			mpdb->chn->id='A';
			mpdb->chn->write(fp);
			mpdb->chn->id=a;
			//
			a=t->mpdb->chn->id;
			t->mpdb->chn->id='B';
			t->mpdb->chn->write(fp);
			t->mpdb->chn->id=a;
			//
			cerr<<"write out strucutre after segment superimposition during composite replacement process:"<<nn<<endl;
			fclose(fp);
		}
	}

	//find how far away the ends
	float rmd1=r1->directrmsdanyway(t1,0,3);
	if(r1->last) {
		cerr<<endl;
		cerr<<"the distance between the starting residues of the segment which are to"<<endl;
		cerr<<" be replaced after superimposition is:"<<rmd1<<endl;
		cerr<<"the name of the two residues: "<<r1->name<<r1->id0<<"--"<<t1->name<<t1->id0<<endl;
		if(rmd1>10) {
			cerr<<endl;	
			cerr<<"Warning! ....."<<endl;
			cerr<<"the distance is large. the composite replacement may not be successful"<<endl;
			cerr<<"you may need to superimpose the segments manually or"<<endl;
			cerr<<"change the start and end values in composite token line"<<endl;
			cerr<<endl;
		}
	}
	float rmd2=r2->directrmsdanyway(t2,0,3);
	if(r2->next) {
		cerr<<endl;
                cerr<<"the distance between the starting residues of the segment which are to"<<endl;
                cerr<<" be replaced after superimposition is:"<<rmd2<<endl;
                cerr<<"the name of the two residues: "<<r2->name<<r2->id0<<"--"<<t2->name<<t2->id0<<endl;
                if(rmd2>10) {
                        cerr<<endl;             
                        cerr<<"Warning! ....."<<endl;
                        cerr<<"the distance is large. the composite replacement may not be successful"<<endl;
			cerr<<"you may need to superimpose the segments manually or"<<endl;
                        cerr<<"change the start and end values in composite token line"<<endl;
                        cerr<<endl;
                }
        }


	//header
	mpdb->header();
	t->mpdb->header();

	//
	float *xyz0,*xyzt,*xyz1,*xyz2; 
	int nlen=strlen(sqnto);

	//set rely
	int i;	
	if(rely) delete [] rely;rely=0;
	if(rely==0) rely=new int[nlen];
	for(i=0;i<nlen;i++) rely[i]=1;

	//change the coordinates
	xyz0=0;
	xyzt=0;
	xyz1=0;
	xyz2=0;
 
	xyz0=mpdb->chn->gettransfer(mpdb->chn->res,1000000);
	xyzt=mpdb->chn->gettransfertemp(mpdb->chn->res,1000000);

	if(r1->last&&r2->next) {		 
		xyz1=t1->chn->gettransfer(t1->next,t2->last->id0);
		xyz2=t1->chn->gettransfertemp(t1->next,t2->last->id0);
		r1->chn->transfer(xyz1,r1->next,r2->last->id0,1);		
		r1->chn->transfertemp(xyz2,r1->next,r2->last->id0,1); 	
		for(i=r1->id0+2;i<=r2->id0-2;i++) rely[i]=2;
	}
	else if(r1->last==0) { 		 
		xyz1=t1->chn->gettransfer(t1,t2->last->id0);
		xyz2=t1->chn->gettransfertemp(t1,t2->last->id0);
		r1->chn->transfer(xyz1,r1,r2->last->id0,1);		
		r1->chn->transfertemp(xyz2,r1,r2->last->id0,1); 	
		for(i=r1->id0;i<=r2->id0-2;i++) rely[i]=2;
	}
	else if(r2->next==0) {
		xyz1=t1->chn->gettransfer(t1->next,t2->id0);
		xyz2=t1->chn->gettransfertemp(t1->next,t2->id0);
		r1->chn->transfer(xyz1,r1->next,r2->id0,1);		
		r1->chn->transfertemp(xyz2,r1->next,r2->id0,1); 	
		for(i=r1->id0+2;i<=r2->id0;i++) rely[i]=2;		  		
	}
	
 	if(TRES.logg>3) mpdb->write("temp.pdb");

	int nend=0;
	if(rmd1<3&&rmd2<3) nend=1;
	nend=1;
	if(r1->last&&r2->next&&nend==0) {
		refine=10;	
		int nlt=2;
		 
		re300:
		int nonerot=0;	
		if(1) {			
			segen->rotatm=new Atm*[1000];
			Res *r;
			int nn=0;
			for(r=r1;r;r=r->next) { 
				if(r->id0>r2->id0) break;
				 
				if(abs(r->id0-r1->id0)<nlt||abs(r->id0-r2->id0)<nlt) {
					segen->rotatm[nn++]=r->isatmid(1);
					segen->rotatm[nn++]=r->isatmid(2);
				}
				else if(rely[r->id0]==1) {
					segen->rotatm[nn++]=r->isatmid(1);
					segen->rotatm[nn++]=r->isatmid(2);
				}	
				else {
					nonerot++;
				}			 						
			}
			segen->rotatm[nn++]=0;					 
		}
		int *stem=new int[2];
		stem[0]=compare[r1->id0];
		stem[1]=compare[r2->id0];
		delsidechain(stem);
		// 
		float coeff=1;

         	if(fapr==5)      coeff=0.25;
         	else if(fapr==4) coeff=0.5;
      		else if(fapr==3) coeff=1;
       		else if(fapr==2) coeff=2;
        	else if(fapr==1) coeff=4;
           	else if(fapr==0) coeff=8;

             	//
        	coeff=2*coeff;
            	//
         	if(r2->id0-r1->id0+1>=10) {
              		segen->arbt=(int)(200*coeff);
            	}
          	else if(r2->id0-r1->id0+1>=8) {
                	segen->arbt=(int)(100*coeff);
            	}
      		else if(r2->id0-r1->id0+1>=6) {
           		segen->arbt=(int)(50*coeff);
            	}
       		else if(r2->id0-r1->id0+1>=4) {
                 	segen->arbt=(int)(30*coeff);
               	}
             	else{
                    	segen->arbt=(int)(20*coeff);
             	}
 
             	if(test) segen->arbt=5;

		//
		float *xyzout=segen->myfixsegment(stem);
		if(stem) delete [] stem;stem=0;
		if(segen->rotatm) {
			delete [] segen->rotatm;
			segen->rotatm=0;
			segen->onlybound=0;
		}	
		if(xyzout) {
			r1->chn->transfer(xyzout,r1,r2->id0,1);
			goto re200;
		}
		else {
			if(nonerot==0) {
				mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
				mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
				cerr<<"the composite replacement is not successful at:"<<endl;
				cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
				goto re200;	
			}	 
			nlt++;
			goto re300;
		}
	} 	

	//fuse the connector of the first one
	if((r1->last&&r2->next==0)||(r1->last&&r2->next)) {

		int *stem=new int[2];
		stem[0]=compare[r1->id0];
		stem[1]=compare[r1->id0];
		Res *r;
		for(r=r1;r;r=r->next) {
			if(r->id0>r2->id0) break;
			if(r->id0-r1->id0>2) break;
			stem[1]=compare[r->id0];
		}
			 	
		while(1) {
			int m1,m2;
                        m1=match[stem[0]];
                        m2=match[stem[1]];
			if(m1==-1||m2==-1) {
				mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
				mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);				 
				//r1->chn->transfer(xyz0,r1,r2->id0,1); 
				//r1->chn->transfertemp(xyzt,r1,r2->id0,1); 
				cerr<<"the composite replacement is not successful at:"<<endl;
				cerr<<"between: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
				goto re200;				 				 	
			}
			Res *rr1=resn[m1];
                        Res *rr2=resn[m2];
			
			if(rr2->id0-rr1->id0>8) {
				mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
				mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
				//r1->chn->transfer(xyz0,r1,r2->id0,1); 
				//r1->chn->transfertemp(xyzt,r1,r2->id0,1); 
				cerr<<"the composite replacement is not successful at:"<<endl;
				cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
				goto re200;
			}
			// 
			float coeff=1;

         		if(fapr==5)      coeff=0.25;
         		else if(fapr==4) coeff=0.5;
      			else if(fapr==3) coeff=1;
       			else if(fapr==2) coeff=2;
        		else if(fapr==1) coeff=4;
           		else if(fapr==0) coeff=8;

             		//
        		coeff=2*coeff;
            		//
         		if(rr2->id0-rr1->id0+1>=10) {
              			segen->arbt=(int)(50*coeff);
            		}
          		else if(rr2->id0-rr1->id0+1>=8) {
                		segen->arbt=(int)(25*coeff);
            		}
      			else if(rr2->id0-rr1->id0+1>=6) {
           			segen->arbt=(int)(15*coeff);
            		}
       			else if(rr2->id0-rr1->id0+1>=4) {
                 		segen->arbt=(int)(10*coeff);
               		}
             		else{
                    		segen->arbt=(int)(5*coeff);
             		}
 
             		if(test) segen->arbt=5;

			//
			float *xyzorg=mpdb->chn->gettransfer(rr1,rr2->id0);			
			segen->updatexyzsave(stem);
			updatepdbcopy();
			delsidechain(stem);
			if(TRES.logg>3) printsegment(stem);
			float *xyzout=segen->myfixsegment(stem);
			if(stem) delete [] stem;stem=0;
			if(xyzout==0) {
				//set sidechain				 
				mpdb->chn->transfer(xyzorg,rr1,rr2->id0,1);				 
				segen->hooksidechain();
				mpdb->chn->transfer(xyzorg,rr1,rr2->id0,1);				 
				mpdb->configure();				
				segen->setallnear();
				if(xyzorg) delete [] xyzorg;xyzorg=0;
				
				//try more
				stem = findnewseg(r1,r2,rr1,rr2,0);
				if(stem==0) {
					mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
					mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
					cerr<<"the composite replacement is not successful at:"<<endl;
					cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
					goto re200;
				}
				continue;
			}
			else {
				if(xyzorg) delete [] xyzorg;xyzorg=0;				
				mpdb->chn->transfer(xyzout,rr1,rr2->id0,1); 
				mpdb->chn->header(rr1->id0-1,rr2->id0+1);
				break;
			}
		}
	}
	
	//fuse the other connector
	if((r2->next&&r1->last==0)||(r1->last&&r2->next)) {
		int *stem=new int[2];
		stem[0]=compare[r2->id0];
		stem[1]=compare[r2->id0];
		Res *r;
		for(r=r2;r;r=r->last) {
			if(r->id0<r1->id0) break;
			if(r2->id0-r->id0>2) break;
			stem[0]=compare[r->id0];
		}
			 	
		while(1) {
			int m1,m2;
                        m1=match[stem[0]];
                        m2=match[stem[1]];
			if(m1==-1||m2==-1) {
				mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
				mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
				//r1->chn->transfer(xyz0,r1,r2->id0,1); 
				//r1->chn->transfertemp(xyzt,r1,r2->id0,1);
				cerr<<"the composite replacement is not successful at:"<<endl;
				cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
				goto re200;
			}
			Res *rr1=resn[m1];
                        Res *rr2=resn[m2];
			if(rr2->id0-rr1->id0>8) {
				mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
				mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
				//r1->chn->transfer(xyz0,r1,r2->id0,1); 
				//r1->chn->transfertemp(xyzt,r1,r2->id0,1); 
				cerr<<"the composite replacement is not successful at:"<<endl;
				cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
				goto re200;
			}
			// 
			float coeff=1;

         		if(fapr==5)      coeff=0.25;
         		else if(fapr==4) coeff=0.5;
      			else if(fapr==3) coeff=1;
       			else if(fapr==2) coeff=2;
        		else if(fapr==1) coeff=4;
           		else if(fapr==0) coeff=8;

             		//
        		coeff=2*coeff;
            		//
         		if(rr2->id0-rr1->id0+1>=10) {
              			segen->arbt=(int)(50*coeff);
            		}
          		else if(rr2->id0-rr1->id0+1>=8) {
                		segen->arbt=(int)(25*coeff);
            		}
      			else if(rr2->id0-rr1->id0+1>=6) {
           			segen->arbt=(int)(15*coeff);
            		}
       			else if(rr2->id0-rr1->id0+1>=4) {
                 		segen->arbt=(int)(10*coeff);
               		}
             		else{
                    		segen->arbt=(int)(5*coeff);
             		}
 
             		if(test) segen->arbt=5;

			//
			float *xyzorg=mpdb->chn->gettransfer(rr1,rr2->id0);			 
			segen->updatexyzsave(stem);
			updatepdbcopy();
			delsidechain(stem);
			if(TRES.logg>3) printsegment(stem);
			float *xyzout=segen->myfixsegment(stem);
			if(stem) delete [] stem;stem=0;
			if(xyzout==0) {
				//set sidechain
				mpdb->chn->transfer(xyzorg,rr1,rr2->id0,1);				 
				segen->hooksidechain();
				mpdb->chn->transfer(xyzorg,rr1,rr2->id0,1);				 
				mpdb->configure();				
				segen->setallnear();
				if(xyzorg) delete [] xyzorg;xyzorg=0;
				 
				//try more
				stem = findnewseg(r1,r2,rr1,rr2,1);
				if(stem==0) {					
					mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
					mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
					cerr<<"the composite replacement is not successful at:"<<endl;
					cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
					goto re200;
				}
				continue;
			}
			else {
				if(xyzorg) delete [] xyzorg;xyzorg=0;				 
				mpdb->chn->transfer(xyzout,rr1,rr2->id0,1);
				mpdb->chn->header(rr1->id0-1,rr2->id0+1);
				break;
			}
		}
	}
	 
        if(TRES.logg>3) r1->chn->write("se0");
	if(segen->bound) delete segen->bound;segen->bound=0;
	//segen->clear();
	cerr<<"minimizing the composite segment..."<<r1->name<<r1->id0<<"--"<<r2->name<<r2->id0<<endl;
	segen->onlysidechain=1;
	segen->readyminfix(r1,r2->id0);
	float d;d=segen->predtmin();
	segen->onlysidechain=0;
	//strcpy(segen->boundtag,"P");	
	//segen->segment->next=segen->mycreate(r1,r2->id0);
	segen->readyminfix(r1,r2->id0);
	d=segen->predtmin();
        //segen->segment->chn->transfer(xyz0,0); 
	//mpdb->chn->transfer(xyz0,r1,r2->id0,1);
	 
	re200:
 	if(TRES.logg>3) mpdb->chn->write("se");
	if(xyz0) delete [] xyz0;
	if(xyzt) delete [] xyzt;
	if(xyz1) delete [] xyz1;
	if(xyz2) delete [] xyz2;
	return 1;
}

int Mutate::assemblemean(Mutate *t,Res *r1,Res *r2) {
	 	
        //set model
	mold=t;
        //
	
	mpdb->setlastresorder();
	mold->mpdb->setlastresorder();
	
	refine=10; //when refine>=10: composite refinement

	//find composite
 	StrFmt *composite=owner->getStrFmt("composite");
 	
	//find ids
	int n1=compare[r1->id0];
	int n2=compare[r2->id0];

	//find corresponding ids in model
	int *stem0=t->findresidstems(r1->id0,r2->id0);  
	
	//
	if(stem0==0) {
		cerr<<"\n\ncomposite replacement of the sequence could not be performed:";
		char *s1=mpdb->chn->getseqn(r1,r2->id0);
		if(s1) {
			cerr<<s1<<endl;
			delete [] s1;s1=0;
		}
		cerr<<endl;
		cerr<<"the reason is no corresponding segment found in the other templates"<<endl;			
		return 0;
	}

	
	int i1=stem0[0];
	int i2=stem0[1];
	delete [] stem0;stem0=0;
	
	//	
	i1=t->match[i1];
	i2=t->match[i2]; 
	Res *t1=t->resn[i1];
	Res *t2=t->resn[i2];	
	char *s1=mpdb->chn->getseqn(r1,r2->id0);
	char *s2=t->mpdb->chn->getseqn(t1,t2->id0);

	//the sequence of the replacement must be equal
	if(s1==0||s2==0||strcmp(s1,s2)!=0) {		
		cerr<<"\n\n\ncomposite replacement could not be performed :"<<endl;			
		cerr<<composite->seqngap<<endl;
		cerr<<"\nthe reason is: "<<endl;
		int i;
		for(i=n1;i<=n2;i++) cerr<<composite->seqngap[i];			
		cerr<<endl;
		cerr<<"does not have identical residues in the two subsequencs between";
		cerr<<" original and the replacement: "<<endl;
		cerr<<"original: \n"<<s1<<endl;
		cerr<<"replace: \n"<<s2<<endl;
		if(s1) delete [] s1;s1=0;	
		if(s2) delete [] s2;s2=0;
		return 0;
		//exit(0);
	}
 	 
	if((strlen(s1)<3&&r1->last&&r2->next)||strlen(s1)==1) {
		cerr<<"the subsequence to be replaced is too few. ignored"<<endl;
		cerr<<"original: \n"<<s1<<endl;
		cerr<<"replace: \n"<<s2<<endl;
		if(s1) delete [] s1;s1=0;	
		if(s2) delete [] s2;s2=0;
		return 0;
	}
	 	
	if(s1) delete [] s1;s1=0;	
	if(s2) delete [] s2;s2=0;
	
	//superimpose the segments or ends
	if(composite->start||composite->end) {
		segexcludesuperimpose(t1,t2,r1,r2);
		if(out>1&&code) {
			char nn[100];
			sprintf(nn,"%s_align.%i.pdb",code,refineid++);
			FILE *fp;
			fp=fopen(nn,"w");
			char a;
			a=mpdb->chn->id;
			mpdb->chn->id='A';
			mpdb->chn->write(fp);
			mpdb->chn->id=a;
			//
			a=t->mpdb->chn->id;
			t->mpdb->chn->id='B';
			t->mpdb->chn->write(fp);
			t->mpdb->chn->id=a;
			//
			cerr<<"write out strucutre after segment superimposition during composite replacement process:"<<nn<<endl;
			fclose(fp);
		}
	}

	//find how far away the ends
	float rmd1=r1->directrmsdanyway(t1,0,3);
	if(r1->last) {
		cerr<<endl;
		cerr<<"the distance between the starting residues of the segment which are to"<<endl;
		cerr<<" be replaced after superimposition is:"<<rmd1<<endl;
		cerr<<"the name of the two residues: "<<r1->name<<r1->id0<<"--"<<t1->name<<t1->id0<<endl;
		if(rmd1>10) {
			cerr<<endl;	
			cerr<<"Warning! ....."<<endl;
			cerr<<"the distance is large. the composite replacement may not be successful"<<endl;
			cerr<<"you may need to superimpose the segments manually or"<<endl;
			cerr<<"change the start and end values in composite token line"<<endl;
			cerr<<endl;
		}
	}
	float rmd2=r2->directrmsdanyway(t2,0,3);
	if(r2->next) {
		cerr<<endl;
                cerr<<"the distance between the starting residues of the segment which are to"<<endl;
                cerr<<" be replaced after superimposition is:"<<rmd2<<endl;
                cerr<<"the name of the two residues: "<<r2->name<<r2->id0<<"--"<<t2->name<<t2->id0<<endl;
                if(rmd2>10) {
                        cerr<<endl;             
                        cerr<<"Warning! ....."<<endl;
                        cerr<<"the distance is large. the composite replacement may not be successful"<<endl;
			cerr<<"you may need to superimpose the segments manually or"<<endl;
                        cerr<<"change the start and end values in composite token line"<<endl;
                        cerr<<endl;
                }
        }


	//header
	mpdb->header();
	t->mpdb->header();

	//
	float *xyz0,*xyzt,*xyz1,*xyz2; 
	int nlen=strlen(sqnto);

	//set rely
	int i;	
	if(rely) delete [] rely;rely=0;
	if(rely==0) rely=new int[nlen];
	for(i=0;i<nlen;i++) rely[i]=1;

	//change the coordinates
	xyz0=0;
	xyzt=0;
	xyz1=0;
	xyz2=0;
 
	xyz0=mpdb->chn->gettransfer(mpdb->chn->res,1000000);
	xyzt=mpdb->chn->gettransfertemp(mpdb->chn->res,1000000);

	if(r1->last&&r2->next) {		 
		xyz1=t1->chn->gettransfer(t1->next,t2->last->id0);
		xyz2=t1->chn->gettransfertemp(t1->next,t2->last->id0);
		r1->chn->transfer(xyz1,r1->next,r2->last->id0,1);		
		r1->chn->transfertemp(xyz2,r1->next,r2->last->id0,1); 	
		for(i=r1->id0+2;i<=r2->id0-2;i++) rely[i]=2;
	}
	else if(r1->last==0) { 		 
		xyz1=t1->chn->gettransfer(t1,t2->last->id0);
		xyz2=t1->chn->gettransfertemp(t1,t2->last->id0);
		r1->chn->transfer(xyz1,r1,r2->last->id0,1);		
		r1->chn->transfertemp(xyz2,r1,r2->last->id0,1); 	
		for(i=r1->id0;i<=r2->id0-2;i++) rely[i]=2;
	}
	else if(r2->next==0) {
		xyz1=t1->chn->gettransfer(t1->next,t2->id0);
		xyz2=t1->chn->gettransfertemp(t1->next,t2->id0);
		r1->chn->transfer(xyz1,r1->next,r2->id0,1);		
		r1->chn->transfertemp(xyz2,r1->next,r2->id0,1); 	
		for(i=r1->id0+2;i<=r2->id0;i++) rely[i]=2;		  		
	}
	
 	if(TRES.logg>3) mpdb->write("temp.pdb");

	int nend=0;
	if(rmd1<3&&rmd2<3) nend=1;
	nend=1;
	if(r1->last&&r2->next&&nend==0) {
		refine=10;	
		int nlt=2;
		 
		re300:
		int nonerot=0;	
		if(1) {			
			segen->rotatm=new Atm*[1000];
			Res *r;
			int nn=0;
			for(r=r1;r;r=r->next) { 
				if(r->id0>r2->id0) break;
				 
				if(abs(r->id0-r1->id0)<nlt||abs(r->id0-r2->id0)<nlt) {
					segen->rotatm[nn++]=r->isatmid(1);
					segen->rotatm[nn++]=r->isatmid(2);
				}
				else if(rely[r->id0]==1) {
					segen->rotatm[nn++]=r->isatmid(1);
					segen->rotatm[nn++]=r->isatmid(2);
				}	
				else {
					nonerot++;
				}			 						
			}
			segen->rotatm[nn++]=0;					 
		}
		int *stem=new int[2];
		stem[0]=compare[r1->id0];
		stem[1]=compare[r2->id0];
		delsidechain(stem);
		// 
		float coeff=1;

         	if(fapr==5)      coeff=0.25;
         	else if(fapr==4) coeff=0.5;
      		else if(fapr==3) coeff=1;
       		else if(fapr==2) coeff=2;
        	else if(fapr==1) coeff=4;
           	else if(fapr==0) coeff=8;

             	//
        	coeff=2*coeff;
            	//
         	if(r2->id0-r1->id0+1>=10) {
              		segen->arbt=(int)(200*coeff);
            	}
          	else if(r2->id0-r1->id0+1>=8) {
                	segen->arbt=(int)(100*coeff);
            	}
      		else if(r2->id0-r1->id0+1>=6) {
           		segen->arbt=(int)(50*coeff);
            	}
       		else if(r2->id0-r1->id0+1>=4) {
                 	segen->arbt=(int)(30*coeff);
               	}
             	else{
                    	segen->arbt=(int)(20*coeff);
             	}
 
             	if(test) segen->arbt=5;

		//
		float *xyzout=segen->myfixsegment(stem);
		if(stem) delete [] stem;stem=0;
		if(segen->rotatm) {
			delete [] segen->rotatm;
			segen->rotatm=0;
			segen->onlybound=0;
		}	
		if(xyzout) {
			r1->chn->transfer(xyzout,r1,r2->id0,1);
			goto re200;
		}
		else {
			if(nonerot==0) {
				mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
				mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
				cerr<<"the composite replacement is not successful at:"<<endl;
				cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
				goto re200;	
			}	 
			nlt++;
			goto re300;
		}
	} 	

	//fuse the connector of the first one
	if((r1->last&&r2->next==0)||(r1->last&&r2->next)) {

		int *stem=new int[2];
		stem[0]=compare[r1->id0];
		stem[1]=compare[r1->id0];
		Res *r;
		for(r=r1;r;r=r->next) {
			if(r->id0>r2->id0) break;
			if(r->id0-r1->id0>2) break;
			stem[1]=compare[r->id0];
		}
			 	
		while(1) {
			int m1,m2;
                        m1=match[stem[0]];
                        m2=match[stem[1]];
			if(m1==-1||m2==-1) {
				mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
				mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);				 
				//r1->chn->transfer(xyz0,r1,r2->id0,1); 
				//r1->chn->transfertemp(xyzt,r1,r2->id0,1); 
				cerr<<"the composite replacement is not successful at:"<<endl;
				cerr<<"between: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
				goto re200;				 				 	
			}
			Res *rr1=resn[m1];
                        Res *rr2=resn[m2];
			
			if(rr2->id0-rr1->id0>8) {
				mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
				mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
				//r1->chn->transfer(xyz0,r1,r2->id0,1); 
				//r1->chn->transfertemp(xyzt,r1,r2->id0,1); 
				cerr<<"the composite replacement is not successful at:"<<endl;
				cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
				goto re200;
			}
			// 
			float coeff=1;

         		if(fapr==5)      coeff=0.25;
         		else if(fapr==4) coeff=0.5;
      			else if(fapr==3) coeff=1;
       			else if(fapr==2) coeff=2;
        		else if(fapr==1) coeff=4;
           		else if(fapr==0) coeff=8;

             		//
        		coeff=2*coeff;
            		//
         		if(rr2->id0-rr1->id0+1>=10) {
              			segen->arbt=(int)(50*coeff);
            		}
          		else if(rr2->id0-rr1->id0+1>=8) {
                		segen->arbt=(int)(25*coeff);
            		}
      			else if(rr2->id0-rr1->id0+1>=6) {
           			segen->arbt=(int)(15*coeff);
            		}
       			else if(rr2->id0-rr1->id0+1>=4) {
                 		segen->arbt=(int)(10*coeff);
               		}
             		else{
                    		segen->arbt=(int)(5*coeff);
             		}
 
             		if(test) segen->arbt=5;

			//
			float *xyzorg=mpdb->chn->gettransfer(rr1,rr2->id0);			
			segen->updatexyzsave(stem);
			updatepdbcopy();
			delsidechain(stem);
			if(TRES.logg>3) printsegment(stem);
			float *xyzout=segen->myfixsegment(stem);
			if(stem) delete [] stem;stem=0;
			if(xyzout==0) {
				//set sidechain				 
				mpdb->chn->transfer(xyzorg,rr1,rr2->id0,1);				 
				segen->hooksidechain();
				mpdb->chn->transfer(xyzorg,rr1,rr2->id0,1);				 
				mpdb->configure();				
				segen->setallnear();
				if(xyzorg) delete [] xyzorg;xyzorg=0;
				
				//try more
				stem = findnewseg(r1,r2,rr1,rr2,0);
				if(stem==0) {
					mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
					mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
					cerr<<"the composite replacement is not successful at:"<<endl;
					cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
					goto re200;
				}
				continue;
			}
			else {
				if(xyzorg) delete [] xyzorg;xyzorg=0;				
				mpdb->chn->transfer(xyzout,rr1,rr2->id0,1); 
				mpdb->chn->header(rr1->id0-1,rr2->id0+1);
				break;
			}
		}
	}
	
	//fuse the other connector
	if((r2->next&&r1->last==0)||(r1->last&&r2->next)) {
		int *stem=new int[2];
		stem[0]=compare[r2->id0];
		stem[1]=compare[r2->id0];
		Res *r;
		for(r=r2;r;r=r->last) {
			if(r->id0<r1->id0) break;
			if(r2->id0-r->id0>2) break;
			stem[0]=compare[r->id0];
		}
			 	
		while(1) {
			int m1,m2;
                        m1=match[stem[0]];
                        m2=match[stem[1]];
			if(m1==-1||m2==-1) {
				mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
				mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
				//r1->chn->transfer(xyz0,r1,r2->id0,1); 
				//r1->chn->transfertemp(xyzt,r1,r2->id0,1);
				cerr<<"the composite replacement is not successful at:"<<endl;
				cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
				goto re200;
			}
			Res *rr1=resn[m1];
                        Res *rr2=resn[m2];
			if(rr2->id0-rr1->id0>8) {
				mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
				mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
				//r1->chn->transfer(xyz0,r1,r2->id0,1); 
				//r1->chn->transfertemp(xyzt,r1,r2->id0,1); 
				cerr<<"the composite replacement is not successful at:"<<endl;
				cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
				goto re200;
			}
			// 
			float coeff=1;

         		if(fapr==5)      coeff=0.25;
         		else if(fapr==4) coeff=0.5;
      			else if(fapr==3) coeff=1;
       			else if(fapr==2) coeff=2;
        		else if(fapr==1) coeff=4;
           		else if(fapr==0) coeff=8;

             		//
        		coeff=2*coeff;
            		//
         		if(rr2->id0-rr1->id0+1>=10) {
              			segen->arbt=(int)(50*coeff);
            		}
          		else if(rr2->id0-rr1->id0+1>=8) {
                		segen->arbt=(int)(25*coeff);
            		}
      			else if(rr2->id0-rr1->id0+1>=6) {
           			segen->arbt=(int)(15*coeff);
            		}
       			else if(rr2->id0-rr1->id0+1>=4) {
                 		segen->arbt=(int)(10*coeff);
               		}
             		else{
                    		segen->arbt=(int)(5*coeff);
             		}
 
             		if(test) segen->arbt=5;

			//
			float *xyzorg=mpdb->chn->gettransfer(rr1,rr2->id0);			 
			segen->updatexyzsave(stem);
			updatepdbcopy();
			delsidechain(stem);
			if(TRES.logg>3) printsegment(stem);
			float *xyzout=segen->myfixsegment(stem);
			if(stem) delete [] stem;stem=0;
			if(xyzout==0) {
				//set sidechain
				mpdb->chn->transfer(xyzorg,rr1,rr2->id0,1);				 
				segen->hooksidechain();
				mpdb->chn->transfer(xyzorg,rr1,rr2->id0,1);				 
				mpdb->configure();				
				segen->setallnear();
				if(xyzorg) delete [] xyzorg;xyzorg=0;
				 
				//try more
				stem = findnewseg(r1,r2,rr1,rr2,1);
				if(stem==0) {					
					mpdb->chn->transfer(xyz0,mpdb->chn->res,10000000,1);
					mpdb->chn->transfertemp(xyzt,mpdb->chn->res,10000000,1);
					cerr<<"the composite replacement is not successful at:"<<endl;
					cerr<<"from: "<<r1->id0<<r1->name<<" "<<r2->id0<<r2->name<<endl;
					goto re200;
				}
				continue;
			}
			else {
				if(xyzorg) delete [] xyzorg;xyzorg=0;				 
				mpdb->chn->transfer(xyzout,rr1,rr2->id0,1);
				mpdb->chn->header(rr1->id0-1,rr2->id0+1);
				break;
			}
		}
	}
	 
        if(TRES.logg>3) r1->chn->write("se0");
	if(segen->bound) delete segen->bound;segen->bound=0;
	//segen->clear();
	cerr<<"minimizing the composite segment..."<<r1->name<<r1->id0<<"--"<<r2->name<<r2->id0<<endl;
	segen->onlysidechain=1;
	segen->readyminfix(r1,r2->id0);
	float d;d=segen->predtmin();
	segen->onlysidechain=0;
	//strcpy(segen->boundtag,"P");	
	//segen->segment->next=segen->mycreate(r1,r2->id0);
	segen->readyminfix(r1,r2->id0);
	d=segen->predtmin();
        //segen->segment->chn->transfer(xyz0,0); 
	//mpdb->chn->transfer(xyz0,r1,r2->id0,1);
	 
	re200:
 	if(TRES.logg>3) mpdb->chn->write("se");
	if(xyz0) delete [] xyz0;
	if(xyzt) delete [] xyzt;
	if(xyz1) delete [] xyz1;
	if(xyz2) delete [] xyz2;
	return 1;
}

void Mutate::assemble0(Mutate *t,Res *r1,Res *r2) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
        segen->flex=1;
        segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
        segen->randcoil=0;
	//TRES.smoothclash=1;
        segen->smoothclash=1;
        segen->part=0;
	segen->flexrot=10;
	segen->arbt=20;
 
	 
	SegBed segbed;
	segbed.next=new SegBed();
	segbed.setsize(1000);
	segbed.next->setsize(1000);
     
	//Res *r;
	refine=10; //when refine>=10: composite refinement

 	StrFmt *composite=owner->getStrFmt("composite");
 	
	int n1=compare[r1->id0];
	int n2=compare[r2->id0];

	int *stem0=t->findresidstems(r1->id0,r2->id0);  
	if(stem0==0) return;
	int i1=stem0[0];
	int i2=stem0[1];
	delete [] stem0;stem0=0;
	
	i1=t->match[i1];
	i2=t->match[i2]; 
	Res *t1=t->resn[i1];
	Res *t2=t->resn[i2];	
	char *s1=r1->chn->getseqn(r1,r2->id0);
	char *s2=t1->chn->getseqn(t1,t2->id0);
	if(s1==0||s2==0||strcmp(s1,s2)!=0) {		
		cerr<<"\n\n\ncomposite replacement could not be performed :"<<endl;			
		cerr<<composite->seqngap<<endl;
		cerr<<"\nthe reason is: "<<endl;
		int i;
		for(i=n1;i<=n2;i++) cerr<<composite->seqngap[i];			
		cerr<<endl;
		cerr<<"does not have identical residues in the two subsequencs between";
		cerr<<" original and the replacement: "<<endl;
		cerr<<"original: \n"<<s1<<endl;
		cerr<<"replace: \n"<<s2<<endl;
		exit(0);
	}
 	 
	if((strlen(s1)<3&&r1->last&&r2->next)||strlen(s1)==1) {
		if(s1) delete [] s1;s1=0;	
		if(s2) delete [] s2;s2=0;
		return;
	}
	 	
	if(s1) delete [] s1;s1=0;	
	if(s2) delete [] s2;s2=0;
	if(autoalgn) {
		Stralg algn;
		algn.flag=1;
		algn.superimposesimple(t1,t2->id0,r1,r2->id0,1);
	}
	t1->chn->header(t1->id0-5,t2->id0+5);
	r1->chn->header(r1->id0-5,r2->id0+5);
	float *xyz0,*xyz1,*xyz2;
	int *stem=new int[2];	 
	int nlen=strlen(sqnto);
	int i;
	for(i=0;i<nlen;i++) rely[i]=1;
	xyz0=0;
	xyz1=0;
	xyz2=0;
 
	xyz0=r1->chn->gettransfer(r1,r2->id0);
	
	if(r1->last&&r2->next) {		 
		xyz1=t1->chn->gettransfer(t1->next,t2->last->id0);
		xyz2=t1->chn->gettransfertemp(t1->next,t2->last->id0);
		r1->chn->transfer(xyz1,r1->next,r2->last->id0,1);		
		r1->chn->transfertemp(xyz2,r1->next,r2->last->id0,1); 	
		stem[0]=compare[t1->id0];
		stem[1]=compare[t2->id0];
		for(i=r1->id0+2;i<=r2->id0-2;i++) rely[i]=2;
	}
	else if(r1->last==0) { 		 
		xyz1=t1->chn->gettransfer(t1,t2->last->id0);
		xyz2=t1->chn->gettransfertemp(t1,t2->last->id0);
		r1->chn->transfer(xyz1,r1,r2->last->id0,1);		
		r1->chn->transfertemp(xyz2,r1,r2->last->id0,1); 
		stem[0]=compare[t1->id0];
		stem[1]=compare[t2->id0];	
		for(i=r1->id0;i<=r2->id0-2;i++) rely[i]=2;
		//for(i=r1->id0;i<=r2->last->id0;i++) rely[i]=2; 	
	}
	else if(r2->next==0) {
		xyz1=t1->chn->gettransfer(t1->next,t2->id0);
		xyz2=t1->chn->gettransfertemp(t1->next,t2->id0);
		r1->chn->transfer(xyz1,r1->next,r2->id0,1);		
		r1->chn->transfertemp(xyz2,r1->next,r2->id0,1); 
		stem[0]=compare[t1->id0];
		stem[1]=compare[t2->id0];
		//test
		//stem[0]=compare[t1->next->id0];
		//end
		for(i=r1->id0+2;i<=r2->id0;i++) rely[i]=2;
		//for(i=r1->next->id0;i<=r2->id0;i++) rely[i]=2;	 		
	}
 
	if(1) {
		segen->rotatm=new Atm*[1000];
		Res *r;
		int nn=0;
		for(r=r1;r;r=r->next) { 
			if(r->id0>r2->id0) break;
			if(r1->last&&r2->next) {
				if(abs(r->id0-r1->id0)<2||abs(r->id0-r2->id0)<3) {
					segen->rotatm[nn++]=r->isatmid(1);
					segen->rotatm[nn++]=r->isatmid(2);
				}
				else if(rely[r->id0]==1) {
					segen->rotatm[nn++]=r->isatmid(1);
					segen->rotatm[nn++]=r->isatmid(2);
				}
			}
			else if(r1->last==0) {
				if(abs(r->id0-r2->id0)<3||rely[r->id0]==1) {
					segen->rotatm[nn++]=r->isatmid(1);
					segen->rotatm[nn++]=r->isatmid(2);
				}
			}
			else if(r2->next==0) {
				if(abs(r->id0-r1->id0)<3||rely[r->id0]==1) {
					segen->rotatm[nn++]=r->isatmid(1);
					segen->rotatm[nn++]=r->isatmid(2);
				}
			}
			
		}
		segen->rotatm[nn++]=0;		
		//segen->onlybound=1;
	}
	if(TRES.logg) mpdb->write("temp.pdb");
 	//segen->onlybound=1;
	float *xyzout=segen->myfixsegment(stem);	 
 	if(segen->rotatm) {
		delete [] segen->rotatm;
		segen->rotatm=0;
		segen->onlybound=0;
	}
	r1->chn->transfer(xyz0,r1,r2->id0,1);

	if(xyz0) delete [] xyz0;xyz0=0;
	if(xyz1) delete [] xyz1;xyz1=0;
	if(xyz2) delete [] xyz2;xyz2=0;

	if(xyzout) {
        		segbed.add(xyzout,segen->start,segen->end);
              		xyzout=0;
			putsegbed(&segbed);                      	             		
        } 
}

int Mutate::isfather(Mutate *t) {
	if(t==0||t->owner==0||owner==0) return 0;

	StrFmt *se=owner->getStrFmt("sequence");
	StrFmt *se0=t->owner->getStrFmt("sequence");
	 
	if(se==0||se0==0) return 0;
 	
	if(strcmp(t->owner->code,owner->code)!=0) return 0;
	if(strcmp(t->owner->seqngap,owner->seqngap)!=0) return 0; 
	if(strcmp(se->seqngap,se0->seqngap)!=0) return 0;
	
	if(t->owner->iscodeexist(code)==0) return 0;	
	
	return 1;		   
}

void Mutate::assembleall(Res *r1,Res *r2) {

	SegBed segbed;	 
	segbed.setsize(1000);
      
	StrFmt *root=owner->getrootStrFmt();
	
	float *xyzorg=mpdb->chn->gettransfer(r1,r2->id0);
		

	StrFmt *t;	 
	int ids=0;
	for(t=root;t;t=t->next) {
		if(t->iscom==1) continue;
		if(t->mutate==0) continue;
		if(t->mutate->mpdb==0) continue;
		if(t->mutate->mpdb->chn==0) continue;
		if(t->mutate==this) continue;	
		if(isfather(t->mutate)==1) continue;

		
		mpdb->chn->transfer(xyzorg,r1,r2->id0,1);	 
		mpdb->chn->header();
		cerr<<endl<<"replace the segment with corresponding segments from model in "<<ids<<"th alignment block..."<<endl;
		int succ=assemble(t->mutate,r1,r2);
		if(t->mutate->mpdb) t->mutate->mpdb->configureatmid0();
		if(mpdb) mpdb->configureatmid0();
		if(succ) {
			float *xyzout=mpdb->chn->gettransfer(r1,r2->id0);
			if(xyzout) {
        			segbed.add(xyzout,r1->id0,r2->id0);
              			xyzout=0;			                    	             		
        		} 
			if(out>1&&code) {
				char nn[100];
				sprintf(nn,"%s_composite_tmp.%i.pdb",code,refineid++);
                		cerr<<endl<<"output intermediate composite structures: "<<nn<<endl;
				cerr<<"segment replaced :"<<r1->name<<r1->id0<<"---"<<r2->name<<r2->id0<<endl<<endl;
				char c=mpdb->chn->id;
        			StrFmt *parent=owner->getparentStrFmt();
        			StrFmt *se = parent->getsequenceStrFmt();
        			mpdb->chn->id=se->cid;
        			if(mpdb->chn->id=='0') mpdb->chn->id=' ';
                		mpdb->write(nn); 
				mpdb->chn->id=c;
			}
		}
		else {
			mpdb->chn->transfer(xyzorg,r1,r2->id0,1);
		}
		ids++;
	}	 
 	mpdb->chn->transfer(xyzorg,r1,r2->id0,1);
 	mpdb->chn->header(r1->id0-1,r2->id0+1);
 
	segen->onlysidechain=1;
	segen->readyminfix(r1,r2->id0);
	float d=segen->predtmin();
	segen->onlysidechain=0;
	//strcpy(segen->boundtag,"P");
	//segen->segment->next=segen->mycreate(r1,r2->id0);
	segen->readyminfix(r1,r2->id0);
	d=segen->predtmin();
	//08/19/2002
	root->setupatmid0();
	//
	mpdb->chn->transfer(xyzorg,r1,r2->id0,0);
 
	segbed.add(xyzorg,r1->id0,r2->id0); 
	xyzorg=0;
	cerr<<endl<<"find the minimal segments among multiple choices: "<<endl<<endl;
	putsegbed(&segbed);   
	if(xyzorg) delete [] xyzorg;xyzorg=0;
}

void Mutate::assembleallmean(Res *r1,Res *r2) {

	SegBed segbed;	 
	segbed.setsize(1000);
      
	StrFmt *root=owner->getrootStrFmt();
	
	float *xyzorg=mpdb->chn->gettransfer(r1,r2->id0);
		

	StrFmt *t;	 
	int ids=0;
	for(t=root;t;t=t->next) {
		if(t->iscom==1) continue;
		if(t->mutate==0) continue;
		if(t->mutate->mpdb==0) continue;
		if(t->mutate->mpdb->chn==0) continue;
		if(t->mutate==this) continue;	
		if(isfather(t->mutate)==1) continue;

		
		mpdb->chn->transfer(xyzorg,r1,r2->id0,1);	 
		mpdb->chn->header();
		cerr<<endl<<"replace the segment with corresponding segments from model in "<<ids<<"th alignment block..."<<endl;
		int succ=assemblemean(t->mutate,r1,r2);
		if(t->mutate->mpdb) t->mutate->mpdb->configureatmid0();	
		if(mpdb) mpdb->configureatmid0();
		if(succ) {
			float *xyzout=mpdb->chn->gettransfer(r1,r2->id0);
			if(xyzout) {
        			segbed.add(xyzout,r1->id0,r2->id0);
              			xyzout=0;			                    	             		
        		} 
			if(out>1&&code) {
				char nn[100];
				sprintf(nn,"%s_composite_tmp.%i.pdb",code,refineid++);
                		cerr<<endl<<"output intermediate composite structures: "<<nn<<endl;
				cerr<<"segment replaced :"<<r1->name<<r1->id0<<"---"<<r2->name<<r2->id0<<endl<<endl;
				char c=mpdb->chn->id;
        			StrFmt *parent=owner->getparentStrFmt();
        			StrFmt *se = parent->getsequenceStrFmt();
        			mpdb->chn->id=se->cid;
        			if(mpdb->chn->id=='0') mpdb->chn->id=' ';
                		mpdb->write(nn); 
				mpdb->chn->id=c;
			}
		}
		else {
			mpdb->chn->transfer(xyzorg,r1,r2->id0,1);
		}
		ids++;
	}	 
 	mpdb->chn->transfer(xyzorg,r1,r2->id0,1);
 	mpdb->chn->header(r1->id0-1,r2->id0+1);
 
	segen->onlysidechain=1;
	segen->readyminfix(r1,r2->id0);
	float d=segen->predtmin();
	segen->onlysidechain=0;
	//strcpy(segen->boundtag,"P");
	//segen->segment->next=segen->mycreate(r1,r2->id0);
	segen->readyminfix(r1,r2->id0);
	d=segen->predtmin();
	//08/19/2002
	root->setupatmid0();
	//
	mpdb->chn->transfer(xyzorg,r1,r2->id0,0);
 
	segbed.add(xyzorg,r1->id0,r2->id0); 
	xyzorg=0;
	cerr<<endl<<"find the minimal segments among multiple choices: "<<endl<<endl;
	putsegbedmean(&segbed);   
	if(xyzorg) delete [] xyzorg;xyzorg=0;
}

float *Mutate::getassemble(Mutate *t,Res *r1,Res *r2) {

	segen->arbt=20;
 
	refine=10; //when refine>=10: composite refinement

 	StrFmt *composite=owner->getStrFmt("composite");
 	
	int n1=compare[r1->id0];
	int n2=compare[r2->id0];

	int *stem0=t->findresidstems(r1->id0,r2->id0);  
	if(stem0==0) return 0;
	int i1=stem0[0];
	int i2=stem0[1];
	delete [] stem0;stem0=0;
	
	i1=t->match[i1];
	i2=t->match[i2]; 
	Res *t1=t->resn[i1];
	Res *t2=t->resn[i2];
	i1=t->resid[i1];
	i2=t->resid[i2];
	Chn *chn=r1->chn;
	r1=0; r2=0;
	r1=(*chn)[i1];
	r2=(*chn)[i2];
	char *s1,*s2; 	
	s1=0; s2=0;

	if(r1&&r2) s1=r1->chn->getseqn(r1,r2->id0);
	if(t1&&t2) s2=t1->chn->getseqn(t1,t2->id0);
	if(s1==0||s2==0||strcmp(s1,s2)!=0) {		
		cerr<<"\n\n\ncomposite replacement could not be performed :"<<endl;			
		cerr<<composite->seqngap<<endl;
		cerr<<"\nthe reason is: "<<endl;
		int i;
		for(i=n1;i<=n2;i++) cerr<<composite->seqngap[i];			
		cerr<<endl;
		cerr<<"does not have identical residues in the two subsequencs between";
		cerr<<" original and the replacement: "<<endl;
		cerr<<"original: \n"<<s1<<endl;
		cerr<<"replace: \n"<<s2<<endl;
		exit(0);
	}	

	float *xyz0=r1->chn->gettransfer(r1,r2->id0);
	float *xyzn=r1->chn->gettransfertemp(r1,r2->id0);

	float *xyz1=t1->chn->gettransfer(t1,t2->id0);
	float *xyz2=t1->chn->gettransfertemp(t1,t2->id0);
	 
	r1->chn->transfer(xyz1,r1,r2->id0,1);		
	r1->chn->transfertemp(xyz2,r1,r2->id0,1); 	

	if(r1->last) t1=r1->last;
	else	     t1=r1;

	if(r2->next) t2=r2->next;
	else	t2=r2;
	
	int *stem=new int[2];
	stem[0]=compare[t1->id0];
	stem[1]=compare[t2->id0];
	int nlen=strlen(sqnto);
	int i;
	for(i=0;i<nlen;i++) rely[i]=1;
	for(i=r1->next->id0;i<=r2->last->id0;i++) rely[i]=2;
	float *xyzout=segen->myfixsegment(stem);	 
 	
	r1->chn->transfer(xyz0,r1,r2->id0,1);
	r1->chn->transfertemp(xyzn,r1,r2->id0,1); 
	if(xyzn) delete [] xyzn;xyzn=0;
	if(xyz0) delete [] xyz0;xyz0=0;
	if(xyz1) delete [] xyz1;xyz1=0;
	if(xyz2) delete [] xyz2;xyz2=0;

	return xyzout;
}

int *Mutate::getmaxrmsd(float *rms,int n1,int n2){

	int i;

	if(n2-n1+1<=window) {
		int *stem=new int[2];
		stem[0]=n1;
		stem[1]=n2;
		return stem;
	}

	int t1,t2;
	float a=0;
	for(i=n1;i<=n2-window+1;i++) {
			
		int j;
		float b=0;
		
		for(j=i;j<=i+window-1;j++) {
			 b+=rms[j];
		}
		if(b>a) {
			a=b;
			t1=i;
			t2=i+window-1;
		}
	}
	
	int *stem=new int[2];
	stem[0]=t1;
	stem[1]=t2;
	return stem;
}

int Mutate::getonlymaxrmsd(float *rms,int n1,int n2){

	int i;

	int   nd=n1;
	float a=rms[n1];
	for(i=n1;i<=n2;i++) {
		if(rms[i]>a) {
			nd=i;
			a=rms[i];
		}	
	}
	return nd;
}

void Mutate::assemble(Res *r1,Res *r2) {

        int i;

        float prev=-1000000;
        int ste=0;
        while(1) {
                ste++;


                float *rms=getavgrmsd(r1,r2);
                Res *r;
                for(r=r1;r&&r->id0<=r2->id0;r=r->next) {
                        cerr<<r->id0<<r->name<<" "<<rms[r->id0]<<endl;
                }
                while(1) {
			int *stem=new int[2];
			stem[0]=compare[r1->id0];
			stem[1]=compare[r2->id0];
			if(autoalgn) algnsegment(stem);
			if(stem) delete [] stem;stem=0;
                        stem=getmaxrmsd(rms,r1->id0,r2->id0);

                        int n1=stem[0];
                        int n2=stem[1];
			if(stem) delete [] stem;stem=0;
                        float d=0;
                        for(i=n1;i<=n2;i++)  d+=rms[i]*rms[i];
                        d=sqrt(d/(n2-n1+1));

                        if(d==0) {                                
                                break;
                        }
                        else if(d<0.1) {                                 
                                for(i=n1;i<=n2;i++) rms[i]=0;
                                continue;
                        }
                        for(i=n1;i<=n2;i++) rms[i]=0;
                        if(stem) delete [] stem;stem=0;
                        Res *t1=mpdb->chn->isres0(n1);
                        Res *t2=mpdb->chn->isres0(n2);
                        assembleall(t1,t2);
                }
                if(rms) delete [] rms;rms=0;
                if(ste>=stepnum) break;
                setcontact(r1,r2->id0);
                float ee=mpdb->chn->getresenergy(r1,r2->id0,5);
                if(fabs(ee-prev)<0.01*prev) break;
                prev=ee;
        }
}
 
void Mutate::assembleflex(Res *r1,Res *r2) {

 

	SegBed segbed;
        //segbed.next=new SegBed();
        segbed.setsize(1000);
       	//segbed.next->setsize(1000);

	StrFmt *t;
	
	StrFmt *root=owner->getrootStrFmt();
	
	float *xyzorg=mpdb->chn->gettransfer(mpdb->chn->res,100000);

	for(t=root;t;t=t->next) {
		if(t->iscom==1) continue;
		if(t->mutate==0||t->mutate->mpdb==0||t->mutate->mpdb->chn==0) continue;
		if(t->mutate==this) continue;
		if(isfather(t->mutate)==1) continue;
		
		int *stem=new int[2];
		stem[0]=compare[r1->id0];
		stem[1]=compare[r2->id0];
		stem=findrightids(t->mutate,stem);
		if(stem==0) continue;
		mpdb->chn->transfer(xyzorg,mpdb->chn->res,100000,1);
		mpdb->header();
		
		assemble(&segbed,t->mutate,stem);
		if(stem) delete [] stem;stem=0;
	}
	//segbed.add(xyzorg,mpdb->chn->res,100000,1);
	putsegbed(&segbed);
	if(xyzorg) delete [] xyzorg;xyzorg=0;
}
 
 
void Mutate::assemble (SegBed *segbed,Mutate *try0,int *stem) 
{
 
	if(segbed==0||try0==0||stem==0) return;

	int n1=stem[0];
	int n2=stem[1];

	int m1=match[n1];
	int m2=match[n2];

	Res *r1=resn[m1];
	Res *r2=resn[m2];

	float *rms=getrmsd(try0,r1,r2);
}
 
float *Mutate::getrmsd(Mutate *try0,Res *r1,Res *r2) {
	float *rms=new float[r2->id0-r1->id0+1000];
	int j=0;
	for(j=0;j<r2->id0-r1->id0+1000;j++) rms[j]=0;

	Res *r;
	for(r=r1;r;r=r->next) {
		if(r->id0>r2->id0) break;
		int n=strlen(try0->seqn);
		int i;
		for(i=0;i<n;i++) {
			if(try0->resid[i]!=r->id0) continue;
			Res *t=try0->resn[i];
			if(t==0) continue;
			rms[r->id0]=r->directrmsdanyway(t,0,3);			 
		}
	}

	return rms;
}

void Mutate::assembleflexold(Res *r1,Res *r2) {

        int i;

        float prev=-1000000;
        int ste=0;
        while(1) {
                ste++;


                float *rms=getavgrmsd(r1,r2);
                Res *r;
                for(r=r1;r&&r->id0<=r2->id0;r=r->next) {
                        cerr<<r->id0<<r->name<<" "<<rms[r->id0]<<endl;
                }
                while(1) {
			int *stem=new int[2];
			stem[0]=compare[r1->id0];
			stem[1]=compare[r2->id0];
			if(autoalgn) algnsegment(stem);
			if(stem) delete [] stem;stem=0;
                        stem=getmaxrmsd(rms,r1->id0,r2->id0);

                        int n1=stem[0];
                        int n2=stem[1];
			if(stem) delete [] stem;stem=0;
                        float d=0;
                        for(i=n1;i<=n2;i++)  d+=rms[i]*rms[i];
                        d=sqrt(d/(n2-n1+1));

                        if(d==0) {                                
                                break;
                        }
                        else if(d<0.1) {                                 
                                for(i=n1;i<=n2;i++) rms[i]=0;
                                continue;
                        }
                        for(i=n1;i<=n2;i++) rms[i]=0;
                        if(stem) delete [] stem;stem=0;
                        Res *t1=mpdb->chn->isres0(n1);
                        Res *t2=mpdb->chn->isres0(n2);
                        assembleall(t1,t2);
                }
                if(rms) delete [] rms;rms=0;
                if(ste>=stepnum) break;
                setcontact(r1,r2->id0);
                float ee=mpdb->chn->getresenergy(r1,r2->id0,5);
                if(fabs(ee-prev)<0.01*prev) break;
                prev=ee;
        }
}

int *Mutate::getresstem(int nsite){
	
	return 0;
}

void Mutate::assemble0(Res *r1,Res *r2) {

	int i;
	
	float prev=-1000000;
	int ste=0;
	while(1) {	
		ste++;	 
		 	 	
		float *rms=getavgrmsd(r1,r2);
		Res *r;
		for(r=r1;r&&r->id0<=r2->id0;r=r->next) {
			cerr<<r->id0<<r->name<<" "<<rms[r->id0]<<endl;
		}
 	 	while(1) {					
			int nsite=getonlymaxrmsd(rms,r1->id0,r2->id0); 
			int *stem=getresstem(nsite);	
			int n1=stem[0];
			int n2=stem[1];
	
			if(stem) delete [] stem;stem=0;
			Res *t1=mpdb->chn->isres0(n1);
			Res *t2=mpdb->chn->isres0(n2);
 
			float d=0;
			for(i=n1;i<=n2;i++)  d+=rms[i]*rms[i];
			d=sqrt(d/(n2-n1+1));
			
			if(d==0) {				 
				break;
			}
			else if(d<0.1) {
				for(i=n1;i<=n2;i++) rms[i]=0;
				continue;
			}
			for(i=n1;i<=n2;i++) rms[i]=0;					 
 
			assembleall(t1,t2);

			stem=new int[2];
			stem[0]=compare[r1->id0];
			stem[1]=compare[r2->id0];
			if(autoalgn) algnsegment(stem);
			delete [] stem;stem=0;
		}
		if(rms) delete [] rms;rms=0;
		if(ste>=stepnum) break;
		setcontact(r1,r2->id0);
		float ee=mpdb->chn->getresenergy(r1,r2->id0,5);
		if(fabs(ee-prev)<0.01*prev) break;
		prev=ee;
	}
}

float *Mutate::getavgrmsd(Res *r1,Res *r2) {
 
	StrFmt *root=owner->getrootStrFmt();
	StrFmt *composite=owner->getStrFmt("composite");
	StrFmt *t;
	float *rms=0;	
	int nlen=strlen(sqnto);
	int i;
	for(i=0;i<nlen;i++) rms=new float[nlen];
	for(i=0;i<nlen;i++) rms[i]=0;
	
	Res *r;
	for(r=r1;r&&r->id0<=r2->id0;r=r->next) {
		float a=0;
		int nn=0;
		for(t=root;t;t=t->next) {		 
			if(t->mutate==0||t->mutate->mpdb==0||t->mutate->mpdb->chn==0) continue;
			if(t->mutate==this) continue;
			int *stem0=t->mutate->findresidstems(r->id0,r->id0);
			if(stem0==0) continue;			
			int i=stem0[0];			 
			delete [] stem0;stem0=0;
			i=t->mutate->match[i];			 
			Res *rt=t->mutate->resn[i];
			float aa=r->directrmsd(rt,1,3);
			if(rt->id0==1) {
				rt->chn->write("s1");
				r->chn->write("s2");
			}
			a+=aa*aa;			
 			nn++;
		}	
		if(nn==0) nn=1;
		a=sqrt(a/nn);
		rms[r->id0]=a;	
	}
	return rms;
}

int *Mutate::getcompositestem(int n){
	
	StrFmt *composite=owner->getStrFmt("composite");
	StrFmt *se=owner->getsequenceStrFmt();
	
	int n1=0;
	for(n1=n;n1>=0;n1--) {
		if(se->seqngap[n1]=='-') continue;
		if(composite->seqngap[n1]==' ')continue;
		if(composite->seqngap[n1]!=composite->seqngap[n]) break;
		
	}
	n1++;

	int n2=0;
	int nlen=strlen(composite->seqngap);
	for(n2=n;n2<nlen;n2++) {
		if(se->seqngap[n2]=='-') continue;
		if(composite->seqngap[n2]==' ')continue;
		if(composite->seqngap[n2]!=composite->seqngap[n]) break;
	}
	n2--;

	if(n1==-1) n1=0;
	if(n2==nlen) n2=nlen-1;
	int i;
	for(i=n1;i<=n2;i++) {
		if(se->seqngap[i]=='-') continue;
		break;
	}
	n1=i;

	for(i=n2;i>=n1;i--) {
		if(se->seqngap[i]=='-') continue;
		break;
	}
	n2=i;	
			
	if(n2<n1) return 0;

	int *stem=new int[2];
	stem[0]=n1;
	stem[1]=n2;
	return stem;
}
 
void Mutate::reindex(Pdb *ss) {

	StrFmt *root=owner->getrootStrFmt();
	StrFmt *se=owner->getStrFmt("sequence");
	Res *r;
	Atm *a;

	int n=0;
	int m=0;
	for(r=ss->chn->res;r;r=r->next) {
		r->id0=n++;
		if(root->onlyrefine==0) r->oid=r->id0+se->start;
		for(a=r->atm;a;a=a->next) {
			a->id0=m++;
			a->id=a->id0;
			if(root->onlyrefine==0) a->oid=a->id+1;
		}		
	}
	ss->chn->start=se->start;
}

void Mutate::reindex() {
	StrFmt *root=owner->getrootStrFmt();
	StrFmt *se=owner->getStrFmt("sequence");
	Res *r;
	Atm *a;

	int n=0;
	int m=0;
	for(r=mpdb->chn->res;r;r=r->next) {
		r->id0=n++;
		if(root->onlyrefine==0) r->oid=r->id0+se->start;
		for(a=r->atm;a;a=a->next) {
			a->id0=m++;
			a->id=a->id0;
		}		
	}
	mpdb->chn->start=se->start;
}

void Mutate::actmutate() {

	if(segen) delete segen;segen=0;
        segen=new Segen();
	segen->flex=1;
	segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
	segen->randcoil=0;
	//TRES.smoothclash=1;
	segen->smoothclash=1; 
 
	 
        //refine: 0
   	  
	refine=0;
	FILE *fp=stderr;
	if(TRES.logg)mpdb->write("start.pdb");
        StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

        int nlen=strlen(sqnto);
        mpdb->setlastresorder();
	
 	int m=0;
	int *stem=0;
	int n=actsidechainmutation(100);
	//mpdb->chn->setbackbonetorsion(mpdb->chn->res,10000);
	
	int insert=0;//int newres=0;
	int isbrk=0; int midr=0;
	 
	while(1) {
		Res **tmp;	
		segen->secd='-';			
		m=findbadsite(); //find the best insertion and deletion region
		if(m==-1) break;	
		midr=m;	
		int tlen=getlength(m);	tlen=2;	
		isbrk=0; 	 
		if(isbreak(m)) { tlen=1000; isbrk=1;}
		tmp=setinitialstructure(m,tlen);	
		tmp=reordertmp(tmp);			 
        	mpdb->setlastresorder();
		linkstructureviaboth(tmp);
		//updatepdbcopy();		
		setreliable(tmp);				
		insert=calcinsert(); 
 		if(tmp==0)   stem=getsegmentworkon(m);		
		else 	     stem=getsegmentworkon(tmp,isbrk);
		if(tmp) delete [] tmp;tmp=0;
		//int ntot=0; 
		while(1) {
			
			int m1,m2;
			m1=match[stem[0]];
			m2=match[stem[1]];
			
			if(m1==-1||m2==-1) break;
			Res *r1=resn[match[stem[0]]];
			Res *r2=resn[match[stem[1]]];
		 	if(r1&&r2&&insert==0&&mpdb->chn->isalllinked(r1,r2->id0)==1) break;  		
			 
			//newres=calcnewres(r1,r2->id0); 			

			if(segen->chiangle) {
                        	delete segen->chiangle;
                        	segen->chiangle=0;
                	}
                	else if(segen->bound) {
                        	delete segen->bound;
                        	segen->bound=0;
                	}						
			//setchiangle(stem);
			if(TRES.logg)printalign(fp,stem);
			int nn=getreslength(stem); 
			if(TRES.logg) printsegment(stem);
			delsidechain(stem);
 
			if(insert>=3) segen->arbt=20;
                        else if(insert>=2) segen->arbt=10;
                        else if(insert>=1) segen->arbt=5;
                        else if(r2->id0-r1->id0+1>5) segen->arbt=5;
                        else segen->arbt=1;

			if(TRES.logg)mpdb->chn->write("semp.pdb");
			m=segen->fixsegment(stem); 
			if(TRES.logg)mpdb->chn->write("temp.pdb");
			if(stem) {delete [] stem; stem=0;}
			if(m) {
				//mpdb->write("s1");
				actsidechainmutation(100);
				//mpdb->write("s2");
				//mpdb->configure();
				//mpdb->chn->setallnear(r1,r2->id0);				
				//minsegment(r1,r2);
				//updatepdbcopy();
				break;	
			}
			stem=findloosesegment(midr,r1,r2);			 
		}
	}	
        actsidechainmutation(-99999);
	minimize(mpdb->chn->res,100000,-100000);       
	char line[100];
	sprintf(line,"%s.pdb",code);
	/*
	char c=mpdb->chn->id;
        StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        mpdb->chn->id=se->cid;
        if(mpdb->chn->id=='0') mpdb->chn->id=' ';
	*/
	mpdb->write(line);
	//mpdb->chn->id=c;
}	
 

void Mutate::setlocaldssp(int *stem) {
	
	Res *r1=resn[match[stem[0]]];
        Res *r2=resn[match[stem[1]]];

	int n;
	for(Res *r=r1;r;r=r->next) {
		if(r->id0>r2->id0) break;
		n=compare[r->id0];
		r->sec=dssp[n];
	}
}

void Mutate::setlocaldssp() {
 	 
	int n;
	for(Res *r=mpdb->chn->res;r;r=r->next) {		 
		n=compare[r->id0];
		r->sec=dssp[n];
	}
}

int *Mutate::checkreslength(int *stem) {
	 
	Res *r1=resn[match[stem[0]]];
        Res *r2=resn[match[stem[1]]];

	int n=0;
	for(Res *r=r1;r;r=r->next) {
		if(r->id0>r2->id0) break;
		n++;
	}
	if(TRES.logg>3) cerr<<"the total residue in the stem: "<<n<<endl;

	if(n>3) return stem;
	if(stem) delete [] stem;stem=0;
	int nn=(r1->id0+r2->id0)/2;
	return findloosesegment(nn,r1,r2);	
}

int Mutate::getreslength(int *stem) {
	 
	Res *r1=resn[match[stem[0]]];
        Res *r2=resn[match[stem[1]]];

	int n=0;
	for(Res *r=r1;r;r=r->next) {
		if(r->id0>r2->id0) break;
		n++;
	}
	if(TRES.logg>3) cerr<<"the total residue in the stem: "<<n<<endl;
	return n;
}

int *Mutate::findloosesegment(int *stem) {

	Res *r1=resn[match[stem[0]]];
        Res *r2=resn[match[stem[1]]];
        int mm=stem[2];
	delete [] stem;stem=0;
	return findloosesegment(r1,r2,mm);	
}

int Mutate::findmidresidue(int *stem) {

	Res *r=getmidresidue(stem);

	if(r) return r->id0;
	else  return -1;
}
 
int Mutate::findmidresidue(Res *a,Res *b) {

        Res *r=getmidresidue(a,b);

        if(r) return r->id0;
        else  return -1;
}

Res *Mutate::getmidresidue(Res *r1,Res *r2) {

        if(r1->last==0||r2->next==0) return 0;

        Chn *c=r1->chn;
        Res *r;
        for(r=r1;r;r=r->next) {
                if(r->id0>=r2->id0) break;
                if(c->islinked(r,r->next)==0) return r;
        }

        int n=(r1->id0+r2->id0)/2;

        for(r=r1;r;r=r->next) {
                if(r->id0>r2->id0) break;
                if(r->id0==n) return r;
        }

        return 0;
}

Res *Mutate::getmidresidue(int *stem) {

        int n1=stem[0];
        int n2=stem[1];

        Res *r1=resn[n1];
        Res *r2=resn[n2];

	if(r1->last==0||r2->next==0) return 0;

	Chn *c=r1->chn;
	Res *r;
        for(r=r1;r;r=r->next) {
                if(r->id0>=r2->id0) break;
		if(c->islinked(r,r->next)==0) return r;
        }

	int n=(r1->id0+r2->id0)/2;

	for(r=r1;r;r=r->next) {
                if(r->id0>r2->id0) break;
                if(r->id0==n) return r;
        }

        return 0;
}



void Mutate::linkstructure(Res **tmp) {
 
	if(tmp==0||tmp[0]==0) return; //delete case
 
	int mm=0;
	while(tmp[mm]) mm++;
 
	Rotate rot;
	 
	Res *r,*r1;
 
	for(int ii=0;ii<mm-1;ii++){		
		r=tmp[ii];
		r1=tmp[ii+1];		
		rot.link(r,r1,r1->id0,1);
	}

	r=tmp[0];r1=tmp[mm-1];

	if(r->last) {	 
		rot.link(r->last,r,r1->id0,1);
	} 	
	else {
		rot.link(r1->next,r,r->id0,0);	
	}
	 
	if(TRES.logg)mpdb->chn->write("mutate4.pdb");
}

void Mutate::linkstructureviaboth(Res **tmp) {
 
	if(tmp==0||tmp[0]==0) return; //delete case
 
	int mm=0;
	while(tmp[mm]) mm++;

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;
 	

	Rotate rot;
	 
	Res *a,*b;
 	Res *r,*r1;
	int n1,n2;
	
	a=tmp[0];b=tmp[mm-1];

	int both=1;

	if(a->last==0&&b->next==0) return;

	if(a->last==0) {
		for(int ii=0;ii<mm-1;ii++){		
			r=tmp[ii];
			r1=tmp[ii+1];		
			rot.link(r,r1,r1->id0,1);
		}
		rot.link(b,b->next,a->id0,0);	
		if(TRES.logg)mpdb->chn->write("mutate4.pdb");
		b=b->next;	 	
		if(TRES.logg>3) cerr<<" the mid residue :"<<a->id0<<" "<<b->id0<<" "<<segen->findmidresidue(a,b)<<endl;
		return;
	}
	else if(b->next==0) {
		for(int ii=0;ii<mm-1;ii++){		
			r=tmp[ii];
			r1=tmp[ii+1];		
			rot.link(r,r1,r1->id0,1);
		}
		rot.link(a->last,a,b->id0,1);	
		if(TRES.logg)mpdb->chn->write("mutate4.pdb");
		a=a->last;		 
		if(TRES.logg>3)cerr<<" the mid residue :"<<a->id0<<" "<<b->id0<<" "<<segen->findmidresidue(a,b)<<endl;
		return;
	}
	else if(1) {

		for(int ii=0;ii<mm-1;ii++){		
			r=tmp[ii];
			r1=tmp[ii+1];		
			rot.link(r,r1,r1->id0,1);
		}

		n1=compare[a->last->id0];		
		n2=compare[b->next->id0];
		int mid=b->id0-a->id0+1; 
		int mcut;
			
		if(parent->sitescore[n1]>parent->sitescore[n2]) {
			if(mid==1||both==0) {
				rot.link(a->last,a,b->id0,1);	
			}	
			else {
				mcut=(int)(a->id0-1+mid/2.+0.6);
				rot.link(a->last,a,mcut,1); 
				rot.link(b,b->next,mcut+1,0);
			}	 
		}
		else {
			if(mid==1||both==0) {
				rot.link(b,b->next,a->id0,0);	
			}	
			else {
				mcut=(int)(a->id0-1+mid/2.-0.4);
				rot.link(a->last,a,mcut,1); 
				rot.link(b,b->next,mcut+1,0);
			}			 	 
		}	
		if(TRES.logg)mpdb->chn->write("mutate4.pdb");
		a=a->last;b=b->next;
		if(TRES.logg>3)mpdb->chn->write("s",a,b->id0); 
		if(TRES.logg>3)cerr<<" the mid residue :"<<a->id0<<" "<<b->id0<<" "<<segen->findmidresidue(a,b)<<endl;
		return;
	}
 

 	a=tmp[0];b=tmp[mm-1];

        int i=0; 
        while(i<mm) {				
		n1=compare[a->last->id0];		
		n2=compare[b->next->id0];					
		if(parent->sitescore[n1]>parent->sitescore[n2]) {
			rot.link(a->last,a,a->id0,1);
			a=a->next;	
			i++;
		}
		else {
			rot.link(b,b->next,b->id0,0);	
			b=b->last;
			i++;
		}			
        }
	a=tmp[0]->last;b=tmp[mm-1]->next;
	 
	cerr<<" the mid residue :"<<n1<<" "<<n2<<" "<<segen->findmidresidue(a,b)<<endl;
	if(TRES.logg)mpdb->chn->write("mutate4.pdb");
}


int Mutate::tunebadsegment(int mm) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return mm;

	int nlen=strlen(se->seqngap);

	if(sqnto[mm]!='-'&&se->seqngap[mm]=='-') return mm; //delete
		
	int m1,m2;

	for(m1=mm-1;m1>=0;m1--) {
		if(sqnto[m1]!='-'&&se->seqngap[m1]!='-') break;									 
	}
		
	for(m2=mm+1;m2<nlen;m2++) {
		if(sqnto[m2]!='-'&&se->seqngap[m2]!='-') break;	 		 	 	
	}	 
 
	Res *ga,*gb;

	if(m1==-1)   ga=mpdb->chn->res;
	else	     ga=resn[match[m1]];
	
	if(m2==nlen) gb=mpdb->chn->lastres();
	else	     gb=resn[match[m2]];

	int i;
	int ii=findloose(ga,gb,mm);
	if(ii==1) i=m2-1;
	else  	  i=m1+1;	

	return i;
}



int Mutate::tunebadsite(int mm) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return mm;

	int nlen=strlen(se->seqngap);

	if(sqnto[mm]!='-'&&se->seqngap[mm]=='-') return mm; //delete
		
	int m1,m2;

	for(m1=mm-1;m1>=0;m1--) {
		if(sqnto[m1]!='-'&&se->seqngap[m1]!='-') break;									 
	}
		
	for(m2=mm+1;m2<nlen;m2++) {
		if(sqnto[m2]!='-'&&se->seqngap[m2]!='-') break;	 		 	 	
	}	 
 
	return findbadsite(m1,m2); 
}

void Mutate::printalign(FILE *fp,int *stem) {

        StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

	fprintf(fp,"\n\n\n\n");

        fprintf(fp,"%s\n",dssp);

        fprintf(fp,"%s\n",sqnto);
	
        fprintf(fp,"%s\n",se->seqngap);
	fprintf(fp,"%s\n",owner->seqngap);
        if(TRES.logg)printsegment(fp,stem);
}

void Mutate::printalign(char *s,int *stem) {

        StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

        FILE *fp=fopen(s,"w");

        fprintf(fp,"%s\n",dssp);

        fprintf(fp,"%s\n",sqnto);

        fprintf(fp,"%s\n",se->seqngap);

	if(TRES.logg) printsegment(fp,stem);
        fclose(fp);
}

void Mutate::printalign(char *s) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

	FILE *fp=fopen(s,"w");

	fprintf(fp,"%s\n",dssp);

	fprintf(fp,"%s\n",sqnto);

	fprintf(fp,"%s\n",se->seqngap);

	fclose(fp);
}

Res **Mutate::setinitialstructure(int mm,int nst) {

	Res **tmp=0;

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	if(TRES.logg) {
	cerr<<dssp<<endl;
	cerr<<sqnto<<endl;
	cerr<<se->seqngap<<endl;
	cerr<<".................."<<endl;
	}

	int nlen=strlen(se->seqngap);

	int m1,m2;

	for(m1=mm-1;m1>=0;m1--) {		 
		if(sqnto[m1]!='-'&&se->seqngap[m1]!='-') break;									 
	}
		
	for(m2=mm+1;m2<nlen;m2++) {		 
		if(sqnto[m2]!='-'&&se->seqngap[m2]!='-') break;	 		 	 	
	}	 

	if(se->seqngap[mm]!='-'&&owner->seqngap[mm]=='-') {
		tmp=new Res*[m2-m1+1];
		for(int j=0;j<m2-m1+1;j++) tmp[j]=0;
	}
 
	int i,jj;//,k;	
	int ii;

	//start using new so that three residues, one between

	int newres[100];
	int nrew=0;
	int both=1;  
	for(jj=m1+1;jj<m2;jj++) {
		
		if(se->seqngap[jj]=='-'&&sqnto[jj]=='-') continue;
		if(sqnto[jj]=='-'||se->seqngap[jj]=='-') {
			if(nrew>90) continue;
			newres[nrew]=jj;
			nrew++;
		}
	}
	if(nst+1==1&&nrew>2&&both==1) {
		i=newres[0];
		ii=nrew/2;		
		newres[0]=newres[ii];
		newres[ii]=i;
	}
	if(nst+1==2&&nrew>2&&both==1) {
		i=newres[1];
		newres[1]=newres[nrew-1];
		newres[nrew-1]=i;
	}
	else if(nst+1==3&&nrew>3&&both==1) {
		i=newres[1];
		newres[1]=newres[nrew-1];
		newres[nrew-1]=i;
		//
		ii=nrew/2;
		i=newres[2];
		newres[2]=newres[ii];
		newres[ii]=i;
	}
	
	// end
 
	int tt=0;
	int nres=0;
	//int med=m1+1;
	//int ned=m2-1;  
	
	
	
	//comment out
	/*
	for(jj=m1+1;jj<m2;jj++) {

		if(tt>nst) break;
		if(se->seqngap[jj]=='-'&&owner->seqngap[jj]=='-') {
			continue;
		}
		else if(se->seqngap[jj]!='-'&&owner->seqngap[jj]=='-'&&m2-m1-1>1) { //insert
			if(med>ned) continue;
			if(m1==-1||m2==nlen||both==0) {				 				
				for(k=med;k<=ned;k++) {
					if(se->seqngap[k]!='-'&&owner->seqngap[k]=='-') break;
				}				
				i=k;
				med=k+1;
			}
			else if(tt%2==0) {				 				 
				for(k=med;k<=ned;k++) {
					if(se->seqngap[k]!='-'&&owner->seqngap[k]=='-') break;
				}
				i=k;
				med=k+1;
				 
			}
			else {				
				for(k=ned;k>=med;k--) {
					if(se->seqngap[k]!='-'&&owner->seqngap[k]=='-') break;
				}				
				i=k; 
				ned=k-1;			 
			}				
			tt++;	
		}
		else {  //delete	
			i=jj;
		}
	*/
	//end
	
	//start new ways
	for(ii=0;ii<nrew;ii++) {

		jj=newres[ii];

		if(tt>nst) break;
		
		if(se->seqngap[jj]=='-'&&owner->seqngap[jj]=='-') {
			continue;
		}
		else if(se->seqngap[jj]!='-'&&owner->seqngap[jj]=='-'&&m2-m1-1>1) { //insert
			
			i=jj;
			tt++;	
		}
		else {  //delete	
			i=jj;
		}
				
	//end
		if(sqnto[i]=='-'&&se->seqngap[i]!='-'&&tmp) {//insert
			int n ;
			for(n=i-1;n>=0;n--) {
				if(sqnto[n]!='-') break; 
			}
			int m;
			for(m=i+1;m<nlen;m++) {
				if(sqnto[m]!='-') break; 
			}
			Res *r0=0,*r=0;
			if(n!=-1) {			
		 		int j=match[n];
				r0=resn[j];
			}	
			if(m!=nlen) {
				int j=match[m];
				r=resn[j];
			}
			Tres *tres_temp;
			tres_temp=TRES[se->seqngap[i]];
			Res *rr=new Res(tres_temp);
			rr->nemp=9;
			rr->temp=new float[9];
			int h;
			for(h=0;h<9;h++) rr->temp[h]=tres_temp->head[h];

				
			if(r0&&r) {				
				r0->next=rr;
				rr->next=r;
				rr->chn=r0->chn;
				rr->id=r0->id+1;
				if(rr->id==rr->next->id) {
					for(Res *t=r;t;t=t->next) t->id++;
				}
			}	
			else if(r0) {				
				r0->next=rr;				 
				rr->chn=r0->chn;
				rr->id=r0->id+1;
			}		 
			else if(r) {				
				mpdb->chn->res=rr;
				rr->next=r;				 		 
				rr->chn=r->chn;
				rr->id=r->id;
				for(Res *t=r;t;t=t->next) t->id++;
			}
 
			sqnto[i]=se->seqngap[i];
			mpdb->configure();
			setmatch();
			tmp[nres++]=rr;
		}
		else if(sqnto[i]!='-'&&se->seqngap[i]=='-'&&tmp==0) {//delete
			int n;
			for(n=i-1;n>=0;n--) {
				if(sqnto[n]!='-') break; 
			}
			int m;
			for(m=i+1;m<nlen;m++) {
				if(sqnto[m]!='-') break; 
			}
			Res *r0=0,*r=0,*r1=0;
			if(n!=-1) {			
		 		int j=match[n];
				r0=resn[j];
			}	
			if(m!=nlen) {
				int j=match[m];
				r=resn[j];
			}
			if( sqnto[i]!='-') {
				int j=match[i];
				r1=resn[j];
			}
			if(r0&&r) {				 
				r0->next=r;
				if(r1) {
					r1->next=0;
					delete r1;r1=0;
				}
				int ni=r->id-r0->id-1;
				Res *t;
				for(t=r;t;t=t->next) t->id-=ni;
			}			 
			else if(r0) {
				r0->next=0;
				if(r1) {
					r1->next=0;
					delete r1;r1=0;
				} 
			}		 
			else if(r) {
				mpdb->chn->res=r;				
				if(r1) {
					r1->next=0;
					delete r1;r1=0;
				} 
				for(Res *t=r;t;t=t->next) t->id--; 
			}
			sqnto[i]=se->seqngap[i];
			mpdb->configure();
			setmatch();	
		}
		
	}
	
	mpdb->configure();
	setmatch();	
	if(TRES.logg) {
	cerr<<dssp<<endl;
	cerr<<sqnto<<endl;
	cerr<<se->seqngap<<endl;
	}
	if(TRES.logg)mpdb->chn->write("mutate3.pdb");
	
	return tmp;
}

int *Mutate::getsegmentworkon(Res **tmp,int isbrk) {
 
	if(tmp==0||tmp[0]==0) return 0;

	//insert case

	int nn=0;
	while(tmp[nn]) nn++;		
	
	Res *r0=tmp[0];
	Res *r1=tmp[nn-1];
	int midr=compare[r0->id0];
	Res *t0=r0->last;
	Res *t1=r1->next;
	//int nlen=strlen(sqnto); 		
	//int mm=compare[r0->id0];

	if(t0==0&&t1) {
		return findloosesegment(midr,r0,t1);
	}	
	else if(t0&&t1==0) {
		return findloosesegment(midr,t0,r1);
	}
	else if(t0==0&&t1==0) {
		return findloosesegment(midr,r0,r1);
	}	
	else if(isbrk) {
		int *stem=new int[2];
		stem[0]=compare[t0->id0];
		stem[1]=compare[t1->id0];
		return stem; 
	} 
	else {				
		return findloosesegment(midr,t0,t1);		
	}
}

int *Mutate::getsegmentworkon(int mm ) {

	//delete case

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);

	int m1,m2;
	
	for(m1=mm-1;m1>=0;m1--) {
		if(sqnto[m1]!='-') break;
		//if(owner->seqngap[m1]!='-'&&se->seqngap[m1]!='-') break;									 
	}
		
	for(m2=mm+1;m2<nlen;m2++) {
		if(sqnto[m2]!='-') break;
		//if(owner->seqngap[m2]!='-'&&se->seqngap[m2]!='-') break;	 		 	 	
	}	 
	
	if(m1==-1) m1=0;
	if(m2==nlen) m2=nlen-1;

	Res *r1,*r2; 
 	 
	if(match[m1]!=-1)  r1=resn[match[m1]];	
	else		   r1=mpdb->chn->res; 
	if(match[m2]!=-1)  r2=resn[match[m2]];
	else		   r2=mpdb->chn->lastres();  
 
	return findloosesegment(mm,r1,r2);		 		
}

int *Mutate::findloosesegment(int midr,Res *r1,Res *r2) {
	
	//find which side is less
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;
 	
	int nlen=strlen(se->seqngap);

	int m1,m2;
 
	m1=compare[r1->id0];
	m2=compare[r2->id0];

	int del=0;

	if(owner->seqngap[midr]!='-'&&se->seqngap[midr]=='-') del=1;
	else del=0;

	int t1,t2;
	for(t2=midr;t2<nlen;t2++){
		if(t2>m2) break; 
		if(del==1&&owner->seqngap[t2]!='-'&&se->seqngap[t2]=='-') continue;
		if(del==0&&owner->seqngap[t2]=='-'&&se->seqngap[t2]!='-') continue;
		
		break;
	}

	for(t1=midr;t1>=0;t1--){
		if(t1<m1) break;
		if(del==1&&owner->seqngap[t1]!='-'&&se->seqngap[t1]=='-') continue;
		if(del==0&&owner->seqngap[t1]=='-'&&se->seqngap[t1]!='-') continue;
		break;
	}

	int p1,p2;
	int ii;
	Res *rr;
	p1=p2=0;
	for(rr=r1;rr;rr=rr->next) {
		if(rr->id0>r2->id0) break;
		ii=compare[rr->id0];
		if(ii>=m1&&ii<=t1) p1++;
		if(ii>=t2&&ii<=m2) p2++;
	}
  

        int w1,w2;
        if(r1->last) {
	    w1=compare[r1->last->id0];
        }
	else w1=-1;

        if(r2->next) {
	    w2=compare[r2->next->id0];
	}
  	else w2=-1; 

	int ge=-1;
	if(w1!=-1&&w2!=-1) {
	    if(dssp[w1]!='-'&&dssp[w2]=='-') {
		ge=1;	
	    }		
	    else if(dssp[w1]=='-'&&dssp[w2]!='-') {
		ge=0;
	    }	
	    else if(dssp[w1]=='h'&&dssp[w2]=='e') {
		ge=0;
	    } 
            else if(dssp[w1]=='e'&&dssp[w2]=='h') {
                ge=1;
            }
	}
	ge=-1;
	
	if(p1>p2) ge=1;
	else if(p1<p2) ge=0;

	float a1=0,a2=0;
 
	a1=parent->sitescore[m1];
	a2=parent->sitescore[m2]; 
	
	if(r1->last==0) {
		m1=compare[r1->id0];
                m2=compare[r2->next->id0];
	}
	else if(r2->next==0) {
		m1=compare[r1->last->id0];
                m2=compare[r2->id0];
	}
        else if(ge==0&&fabs(a1-a2)<0.1) {
		m1=compare[r1->last->id0];
                m2=compare[r2->id0];

	}
        else if(ge==1&&fabs(a1-a2)<0.1) {
		m1=compare[r1->id0];
                m2=compare[r2->next->id0];   
	}
	else if(a1<a2) {
		m1=compare[r1->last->id0];
		m2=compare[r2->id0];
	}
	else {
		m1=compare[r1->id0];
		m2=compare[r2->next->id0];
	}

	int *out=new int[2];
	out[0]=m1;out[1]=m2;
	return out;
}

int *Mutate::findloosesegment(Res *r1,Res *r2,int n) {
	 
	int i;
	Res *t1=r1;
	Res *t2=r2;
	for(i=0;i<n;i++) {		
		int *stem=findloosesegment(n,t1,t2);
		if(i==n-1) return stem;
		int n1=match[stem[0]];
		int n2=match[stem[1]];
		t1=resn[n1];t2=resn[n2];
		delete [] stem;	stem=0;
	}
	int *stem=new int[2];
	stem[0]=compare[r1->id0];
	stem[1]=compare[r2->id0];
	return stem;
}



int *Mutate::findloosesegment0(Res *r1,Res *r2,int mm) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	Res *rr;
	int m1,m2;
	int nlen=strlen(sqnto); 
	
	float a1=0,a2=0;
	int n=0;
	for(rr=r1;rr;rr=rr->last) {
		if(rr->sec!=r1->sec) break;
		if(rr->sec=='-') a1+=10000;
		else if(rr->sec=='h') a1+=100;
		else if(rr->sec=='e') a1+=1;
		n=compare[rr->id0];
		if(se->seqngap[n]=='-'&&owner->seqngap[n]!='-') {
			a1+=20000;
		}
		else if(se->seqngap[n]!='-'&&owner->seqngap[n]=='-') {
			a1+=10000;
                }
		else {
			a1+=100;
		}

	} 
	if(rr==0) a1+=100000;

	for(rr=r2;rr;rr=rr->next) {
		if(rr->sec!=r2->sec) break;
		if(rr->sec=='-') a2+=10000;
		else if(rr->sec=='h') a2+=100;
		else if(rr->sec=='e') a2+=1;
		
		n=compare[rr->id0];
                if(se->seqngap[n]=='-'&&owner->seqngap[n]!='-') {
                        a2+=20000;
                }
                else if(se->seqngap[n]!='-'&&owner->seqngap[n]=='-') {
                        a2+=10000;
                }
                else {
                        a2+=100;
                }
	} 
	if(rr==0) a2+=100000;

	int n1=compare[r1->id0];
	int n2=compare[r2->id0];

	if(parent->sitescore[n1]>parent->sitescore[n2]) {
		a2+=(parent->sitescore[n1]-parent->sitescore[n2])*10000;	
	}
	else {
		a1+=(parent->sitescore[n2]-parent->sitescore[n1])*10000;
	}

	if(r1->last==0) {
		m1=compare[r1->id0];
                m2=compare[r2->next->id0];
	}
	else if(r2->next==0) {
		m1=compare[r1->last->id0];
                m2=compare[r2->id0];
	}
	else if(a1>a2) {
		m1=compare[r1->last->id0];
		m2=compare[r2->id0];
	}
	else {
		m1=compare[r1->id0];
		m2=compare[r2->next->id0];
	}

	int *out=new int[4];
	out[0]=m1;out[1]=m2;out[2]=mm;out[3]=-nlen; 
	return out;
}

int Mutate::findloose(Res *r1,Res *r2,int mm) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	Res *rr;
	//int m1,m2;
	int nlen=strlen(sqnto); 
	
	float a1=0,a2=0;

	if(r1->last==0) return 1;
	if(r2->next==0) return 0;

	for(rr=r1;rr;rr=rr->last) {
		if(rr->sec!=r1->sec) break;
		if(rr->sec=='-') a1+=10000;
		else if(rr->sec=='h') a1+=100;
		else if(rr->sec=='e') a1+=1;
	} 

	for(rr=r2;rr;rr=rr->next) {
		if(rr->sec!=r2->sec) break;
		if(rr->sec=='-') a2+=10000;
		else if(rr->sec=='h') a2+=100;
		else if(rr->sec=='e') a2+=1;
	} 

	if(a1>a2) {
		return 0;
	}
	else {
		return 1;
	}
}

int Mutate::findloose(Res *r1,Res *r2) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	Res *rr;
	//int m1,m2;
	int nlen=strlen(sqnto); 
	
	float a1=0,a2=0;

	if(r1->last==0) return 1;
	if(r2->next==0) return 0;

	for(rr=r1;rr;rr=rr->last) {
		if(rr->sec!=r1->sec) break;
		if(rr->sec=='-') a1+=10000;
		else if(rr->sec=='h') a1+=100;
		else if(rr->sec=='e') a1+=1;
	} 

	for(rr=r2;rr;rr=rr->next) {
		if(rr->sec!=r2->sec) break;
		if(rr->sec=='-') a2+=10000;
		else if(rr->sec=='h') a2+=100;
		else if(rr->sec=='e') a2+=1;
	} 

	if(a1>a2) {
		return 0;
	}
	else {
		return 1;
	}
}
//stop.....

int *Mutate::findsegmentworkon(int mm ) {
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);

 	int m1,m2;
	int n1,n2;
	
	for(m1=mm-1;m1>=0;m1--) {
		if(owner->seqngap[m1]!='-'&&se->seqngap[m1]!='-') break;									 
	}
		
	for(m2=mm+1;m2<nlen;m2++) {
		if(owner->seqngap[m2]!='-'&&se->seqngap[m2]!='-') break;	 		 	 	
	}	 
 
	int n=mm;
	for(n1=n-1;n1>=0;n1--) {
		if(dssp[n1]==dssp[n]) continue;
		if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;									 
	}
		
	for(n2=n+1;n2<nlen;n2++) {
		if(dssp[n2]==dssp[n]) continue;
		if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;	 		 	 	
	}	
 	
	
	int p1=0;
	int p2=0; 
	int i=0;
		
	for(i=n1+1;i<n2;i++) {			
		if(owner->seqngap[i]!='-') p1++;
		if(se->seqngap[i]!='-') p2++;
	}
	
	if(abs(p1-p2)<=2&&(p2>0&&p1>0&&p2<6)&&dssp[n]=='-') {		
		int *out=new int[4];
		out[0]=n1;out[1]=n2;out[2]=n;out[3]=-nlen; 
		return out;				 		 		
	}
	
	p1=p2=0;
	

	for(m1=n-1;m1>=n1;m1--) {
		if(owner->seqngap[m1]!='-'&&se->seqngap[m1]!='-') break;									 
	}
		
	for(m2=n+1;m2<n2;m2++) {
		if(owner->seqngap[m2]!='-'&&se->seqngap[m2]!='-') break;	 		 	 	
	}	 	

	for(i=m1+1;i<m2;i++) {			
		if(sqnto[i]!='-') p1++;
		if(se->seqngap[i]!='-') p2++;
	}

	if(p2-p1>0) {//insert
		int *out=new int[4];		
		out[0]=m1;out[1]=m2;out[2]=n;out[3]=0; 
		return out;	
	}
	else {	//delete
		int *out=new int[4];		
		out[0]=m1;out[1]=m2;out[2]=n;out[3]=0; 
		return out;			
	}

	return 0;	 
}


int Mutate::findbadsegment() {

//this subroutine is to find those unlinked on the loop region. 
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);

	int n1,n2;
	int p1,p2;
	int m1,m2,m,i;
 
	float dd=-1;
	float e=0;
 	m=-1;
	for(int n=0;n<nlen;n++) { 

		if(sqnto[n]=='-'&&se->seqngap[n]=='-') continue;
		if(sqnto[n]!='-'&&se->seqngap[n]!='-') continue;
 
		//find nearest secondary structure
		e=0;			
	 
		if(dssp[n]=='-') e+=0;
		else if(dssp[n]=='h') e+=10000;
		else if(dssp[n]=='e') e+=20000;
	  
		if(sqnto[n]!='-'&&se->seqngap[n]=='-') e+=50000; //deletion last
		 

		for(n1=n-1;n1>=0;n1--) {
			if(dssp[n1]==dssp[n]) continue;
			if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;			 
		}
		
		for(n2=n+1;n2<nlen;n2++) {
			if(dssp[n2]==dssp[n]) continue;
			if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;	 		 	 	
		}	 
		//end
		
		if(n1==-1||n2==nlen) e+=30000;
 		else 	e+=fabs((n2-n1)-5)*1000;

		//find nearest residue	  
		for(m1=n-1;m1>=n1;m1--) {			 
			if(sqnto[m1]!='-'&&se->seqngap[m1]!='-') break;									 
		}
		
		for(m2=n+1;m2<=n2;m2++) {
			if(sqnto[m2]!='-'&&se->seqngap[m2]!='-') break;	 		 	 	
		}	 
		//end
		
		e+=10*fabs((m2-m1)-3);
			
		//detect the number of sequence in the region before and after mutation	
		p1=0;p2=0;
 		
		for(i=m1+1;i<m2;i++) {			
			if(sqnto[i]!='-') p1++;
			if(se->seqngap[i]!='-') p2++;
		}
		//end
		if(p2>p1) e+=(p2-p1)*50;	
		else 	  e+=(p1-p2)*100;
		if(dd<0||dd>e) {
			dd=e;
			m=n;
		}		

		//start
		p1=0;p2=0; 		
		for(i=n1+1;i<n2;i++) {			
			if(owner->seqngap[i]!='-') p1++;
			if(se->seqngap[i]!='-') p2++;
		}

		if(abs(p1-p2)<=2&&(p2>0&&p1>0&&p2<6)&&e<10000) {			 
			return n;		 		
		}					
		//end
	} 
 	
	return m;
}

int Mutate::findbadsite(int n1,int n2) {

//this subroutine is to find those unlinked on the loop region.

        StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

        int m;

        float dd=0;
        float e=0;
        m=-1;
        for(int n=n1;n<=n2;n++) {

                if(sqnto[n]=='-'&&se->seqngap[n]=='-') continue;
                if(sqnto[n]!='-'&&se->seqngap[n]!='-') continue;

                e=parent->sitescore[n];

                if(m<0) {
                        dd=e;
                        m=n;
                }
                else if(dd<e) {
                        dd=e;
                        m=n;
                }
        }
	if(TRES.logg>3) cerr<<"the score at site: "<<m<<" is: "<<dd<<endl;
        return m;
}

int Mutate::isbreak(int m) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);
 
	int n1,n2;
	for(n1=m-1;n1>=0;n1--) {
		if(sqnto[n1]!='-'&&se->seqngap[n1]!='-') break;
	}

	for(n2=m+1;n2<nlen;n2++) {
		if(sqnto[n2]!='-'&&se->seqngap[n2]!='-') break;
	}
	if(n1==-1||n2==nlen) return 0;

	n1=match[n1];n2=match[n2];

	Res *r1=resn[n1];
	Res *r2=resn[n2];

	if(r1->next==r2&&r2->id-r1->id>1&&r1->chn->islinked(r1,r2)==0) return 1;
	else return 0;
}

void Mutate::findrealstem(int m,int *stem) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) {stem[0]=-1;stem[1]=-1;return ;}

	int nlen=strlen(sqnto);
 
	int n1,n2;
	for(n1=m-1;n1>=0;n1--) {
		if(sqnto[n1]!='-'&&se->seqngap[n1]!='-') break;
	}

	for(n2=m+1;n2<nlen;n2++) {
		if(sqnto[n2]!='-'&&se->seqngap[n2]!='-') break;
	}
	if(n1==-1||n2==nlen) {
		stem[0]=-1;
		stem[1]=-1;
		return;
	}
	 
	int tot2=0;

	int i;
	for(i=n1;i<=n2;i++) {		 
		if(se->seqngap[i]!='-')tot2++;
	}

	int d1,d2;

	d1=match[n1];d2=match[n2];

	Res *r1=resn[d1];
	Res *r2=resn[d2];

	if(r2->id-r1->id+1==tot2) {
		stem[0]=n1;
		stem[1]=n2;
		return;
	}
	else {
		stem[0]=-1;
		stem[1]=-1;
		return;
	}
}

int Mutate::findbadsite() {

//this subroutine is to find those unlinked on the loop region. 
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);
 
	int m;
 
	Res *r;

	for(r=mpdb->chn->res;r;r=r->next) {
		if(r->next==0) continue;
		if(r->next->id-r->id>1) {
			Res *r2=r->next;
			Res *r1=r;
			int n2=compare[r2->id0];
			int n1=compare[r1->id0];
			int n;
			int i=0;
			for(n=n1+1;n<n2;n++) {
				if(se->seqngap[n]!='-') {
					m=n;
					i++;
				}
			}
			if(i) return m;
		}
	}

	float dd=0;
	float e=0;
 	m=-1;
	for(int n=0;n<nlen;n++) { 

		if(sqnto[n]=='-'&&se->seqngap[n]=='-') continue;
		if(sqnto[n]!='-'&&se->seqngap[n]!='-') continue;
 
		e=parent->sitescore[n];
 		
		if(m<0) {
			dd=e;
			m=n;
		}
		else if(dd<e) { //find the most conserved
			dd=e;
			m=n;
		}		
	} 
 	if(TRES.logg>3) cerr<<"the score at site: "<<m<<" is: "<<dd<<endl;
	return m;
}

int Mutate::findbadsite(int *stem) {

//this subroutine is to find those unlinked on the loop region. 
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);
 
	int m;
 
	Res *r;

	for(r=mpdb->chn->res;0&&r;r=r->next) {
		if(r->next==0) continue;
		if(r->next->id-r->id>1) {
			Res *r2=r->next;
			Res *r1=r;
			int n2=compare[r2->id0];
			int n1=compare[r1->id0];
			int n;
			int i=0;
			for(n=n1+1;n<n2;n++) {
				if(se->seqngap[n]!='-') {
					m=n;
					i++;
				}
			}
			if(i) return m;
		}
	}

	float dd=0;
	float e=0;
 	m=-1;
	for(int n=stem[0];n<=stem[1];n++) { 

		if(sqnto[n]=='-'&&se->seqngap[n]=='-') continue;
		if(sqnto[n]!='-'&&se->seqngap[n]!='-') continue;
 
		e=parent->sitescore[n];
 		
		if(m<0) {
			dd=e;
			m=n;
		}
		else if(dd<e) { //find the most conserved
			dd=e;
			m=n;
		}		
	} 
 	if(TRES.logg>3) cerr<<"the score at site: "<<m<<" is: "<<dd<<endl;
	return m;
}

void  Mutate::setdloopstem(int m,int *stem) {

//this subroutine is to find those unlinked on the loop region. 
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

	int nlen=strlen(sqnto);
 
	int n1,n2;

	for(n1=m;n1>=0;n1--) {		
		if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;
	}

	if(n1==-1) n1=0;
	for(n2=m;n2<nlen;n2++) {		
		if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;
	}
	if(n2==nlen) n2=nlen-1;
	
	stem[0]=n1;stem[1]=n2;
}

int  Mutate::calctotalinsert(int m) {

//this subroutine is to find those unlinked on the loop region. 
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);
 
	int n1,n2;

	for(n1=m;n1>=0;n1--) {		
		if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;
	}

	if(n1==-1) n1=0;
	for(n2=m;n2<nlen;n2++) {		
		if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;
	}
	if(n2==nlen) n2=nlen-1;
	
	int nn=0;
	int i;
	for(i=n1;i<=n2;i++) {
		if(owner->seqngap[i]=='-'&&se->seqngap[i]=='-') continue;
		if(owner->seqngap[i]!='-'&&se->seqngap[i]!='-') continue;
		if(se->seqngap[i]!='-') nn++;
		if(owner->seqngap[i]!='-') nn--;
	}
	return nn;
}
int Mutate::checkdloop(int m) {
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);
 
	int stem[2];
	setdloopstem(m,stem);

	m=0;
	for(int n=stem[0];n<=stem[1];n++) { 

		if(sqnto[n]=='-'&&se->seqngap[n]=='-') continue;
		if(sqnto[n]!='-'&&se->seqngap[n]!='-') continue;
 		return 1;
	} 
 	//cerr<<"the score at site: "<<m<<" is: "<<dd<<endl;
	return 0;
}




int Mutate::actsinglemutation(Res *r,char c,float f) {
  
 	float xyz[100];
 
	int tot=0;
	int i=compare[r->id0];

 	mpdb->setflgr(-99999);
	r->transfer(xyz,0);	 	
	r->mutateResidue(c);
	
	mpdb->configure();
	if(strchr("AGP",r->name)) {					
		if(TRES.logg) cerr<<"residue:"<<sqnto[i]<<r->id<<" is mutated to:"<<c<<endl;
		sqnto[i]=c;
		tot++;	
		return tot;
	}
	r->flag=10000; 
	Scap scprd; 	 		
	//charge
	strcpy(scprd.force,"124");
	TRES.setdonaldcharge();
	scprd.dielectric=15;
	//	
	scprd.pdb=mpdb;
		 
	float d=scprd.mutate();
	if(f<-9999||d<f) { 
                if(TRES.logg) cerr<<"residue:"<<sqnto[i]<<r->id<<" is mutated to:"<<c<<" with energy:"<<d<<endl;
                sqnto[i]=c;
		mpdb->configure();
                tot++;
        }
	else { //no change at all
		r->mutateResidue(sqnto[i]);
		mpdb->configure();
		r->transfer(xyz,1);
	}
	r->flag=-99999;	
	r->setfreenextonmore();	
	if(r->more) delete r->more;r->more=0;
	//restore charge	
	TRES.switchcharge(1);
	//
	return tot;			 					 
} 



int Mutate::actsidechainmutation(float f) {
 
	//get the sequence class of this 

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return -1;

	int nlen=strlen(se->seqngap);

	if(TRES.logg>3)mpdb->chn->write("mutate0.pdb");

	int i;
 	//float xyz[100];

	mpdb->setflgr(-99999);
	int tot=0;
	char c;
	//int tot0=0;
        for(i=0;i<nlen;i++) {
		if(sqnto[i]=='-'||se->seqngap[i]=='-') continue;
		if(se->seqngap[i]==sqnto[i])  {			
			continue;
		}
		
		//int n=owner->match[i];
		int n=match[i];
		c=se->seqngap[i];
		Res *r=resn[n];
		
		if(r==0) {
			cerr<<"the residue does not exist:"<<n<<endl;
			continue;
		}		 
		tot+=actsinglemutation(r,c,f);				 				
	}
	setmatch();
	if(TRES.logg>3) {
 		cerr<<getsecondary()<<endl;
		cerr<<sqnto<<endl;
		cerr<<se->seqngap<<endl;
		if(TRES.logg)mpdb->chn->write("mutate.pdb");
	}	
	return tot;
} 

void Mutate::fullsidechain(int n) {
	
	mpdb->setflgr(10000);
	mpdb->setflgr('A',-99999);
	mpdb->setflgr('P',-99999);
	mpdb->setflgr('G',-99999);	

	Scap scprd; 
	//08/16/2002
	TRES.setdonaldcharge();
	strcpy(scprd.force,"124");
	scprd.dielectric=15;
	//end
	scprd.includeself=1;
	scprd.singletorsion=1;
	scprd.colonyline=1;
	scprd.colony=2;
	scprd.ncolony=1;
	scprd.nncolony=1;
	scprd.nummore=100;
	scprd.bmax=2;
	scprd.tormax=2;
	scprd.ring=1;	
	scprd.pdb=mpdb;
	scprd.scpred(n); 
	
	for(Chn *c=mpdb->chn;c;c=c->next) c->setfreenextonmore();
	//restore charge
	TRES.switchcharge(1);
	//
	if(TRES.logg)mpdb->write("side");
}

void Mutate::domutate(){
 	setdssp();
	dosidechainmutation(2); //easy part;
	createinitialstructure(); //easy part
	linkinitialstructure();   //easy part
	zipgappedstructure();    //difficult
	//dooutlierinsertion();
	//dooutlierdeletion();	 
	//linkgappedstructure();
	//doeditloop(115);	
	//doeditresidue(8);
}
int Mutate::findlongestreliable(Res *r) {
	
	Res *a,*t;

	t=r;
	for(a=r;a;a=a->next) {
		if(a==0) continue;
		if(rely[a->id0]==0) return t->id0;
		t=a;
	}

	return t->id0;
}

int Mutate::findminuslongestreliable(Res *r) {
 
	Res *a,*t;

	int n=r->id0;
	t=r;
	for(int i=n;i>=0;i--) {
		a=r->chn->isres0(i);
		if(a==0) continue;
		if(rely[a->id0]==0) return t->id0;
		t=a;
	}

	return t->id0;
}


void Mutate::delsidechain() {
 
	Res *r;
	for(r=mpdb->chn->res;r;r=r->next) {
		
		if(rely[r->id0]==1) continue;		
		mpdb->chn->transform(r,r->id0,3);
	}
	setmatch();
	mpdb->configure();
	if(TRES.logg)mpdb->write("s5");
}
void Mutate::delsidechain(Res **tmp) {

	if(tmp==0||tmp[0]==0) return;
	int mm=0;
	while(tmp[mm])mm++;
	for(int i=0;i<mm;i++) { 
		Res *r=tmp[i]; 
		mpdb->chn->transform(r,r->id0,3);
	}
	mpdb->configure();
	if(TRES.logg)mpdb->write("del.pdb");
}


void Mutate::delsidechain(int *stem) {
	 
	int m1=stem[0];
	int m2=stem[1];
	
	int n1=match[m1];
	int n2=match[m2];

	Res *r1=resn[n1];
	Res *r2=resn[n2];

	mpdb->chn->transform(r1,r2->id0,3); 

	mpdb->configure();
	if(TRES.logg)mpdb->write("del.pdb");
}


void Mutate::delsidechain(int mm) { 
	Res *r;
	int i=match[mm];
	r=resn[i];
	if(r==0) return;
	mpdb->chn->transform(r,r->id0,3);
	mpdb->configure();
	if(TRES.logg)mpdb->write("s5");
}

void Mutate::zipgappedstructure(){

	if(segen) delete segen;segen=0;
	segen=new Segen();;
	
  	segen->pdb=mpdb;
	segen->cid=mpdb->chn->id;
	segen->mutate=this; 

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

	int nlen=strlen(sqnto);
	mpdb->setlastresorder();
	setreliable();
	delsidechain();

	int *stem=0;

	while(1) {
		 
		//find easy loops

		if(segen->chiangle) {
			delete segen->chiangle;
			segen->chiangle=0;
		}		
		else if(segen->bound) {
			delete segen->bound;
			segen->bound=0;
		}

		stem=findlooptoworkon();	
		 
		int n=segen->fixsegment(stem);	
		 
		if(stem) {
			delete [] stem; stem=0;	
		}	 
		 
	}
	
}


int *Mutate::findlooptoworkon() {

//this subroutine is to find those unlinked on the loop region. 

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	//cflg=0;

	int nlen=strlen(sqnto);
	
	Res *r;
	
	//int m1=0,m2=0,m3=0;

	//float d=10000;

	//mpdb->write("ss");

	//tt=0: 1 or 2 residue difference, no secondary structure frame shift,
	//	less than 7 residues
	//tt=1: 1 or 2 residue difference, no secondary structure frame shift,
	//	more than 7 residues, stepwise
	//tt=2: 1 or 2 residue difference, zero occurs
	//tt=3: more than 3 residue, consider frame shift, less than 7 residue
	//tt=4: more than 3 residue, consider frame shift, more than 7 residue, stepwise
	//tt=5: deletion/insertion occur in secondary structure
	//tt=6: occur in head and tail.

	int tt=0;
	int itern=0;
	re200:
	//int dis=1000;
	//m1=-1;m2=-1;m3=-1;
	itern=0;
	Res *rr=0;
	for(r=mpdb->chn->res;r;r=r->next) {

		if(r->next==0) continue;
		
		int n=compare[r->id0];	
		
		int nn=compare[r->next->id0];
		
		if(dssp[n]!='-'&&dssp[nn]!='-') continue;
		 	
		if(mpdb->chn->islinked(r,r->next)==1)continue;	
 		
 		itern++;

		int n1,n2;
		int p1=0,p2=0;

		//the following code to detect if zig-zag alignment exists
		int window;
		for(window=1;window<5;window++) {
				
			p1=0;p2=0;
			for(n1=n-window;n1<=n+window+1;n1++) {
		 		if(owner->seqngap[n1]!='-') p1++;
				if(sqnto[n1]!='-')p2++; 
			}
			if(p1==p2) {
				int *out=new int[4];
				out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 				 
				return out;
			}		 
		}		
		//end
		
		//find nearest secondary structure
	  
		for(n1=n-1;n1>=0;n1--) {
			if(dssp[n1]=='-') continue;
			//if(dssp[n1]==dssp[n]) continue;
			if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;									 
		}
		
		for(n2=n+1;n2<nlen;n2++) {
			if(dssp[n2]=='-') continue;
			//if(dssp[n2]==dssp[n]) continue;
			if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;	 		 	 	
		}
		 
		//end

		//not exist, tails
		if(n1==-1) continue;		
		if(n2==nlen) continue;		
			
		//detect the number of sequence in the region before and after mutation	
		p1=0;p2=0;
		
		int i;
		for(i=n1+1;i<n2;i++) {
			if(owner->seqngap[i]!='-') p1++;
			if(sqnto[i]!='-') p2++;
		}
		//end

		//simple case,zig and zag
		if(p1==p2) {
			int *out=new int[4];
			out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 
			return out;
		}		
		//end
		
		//only one difference
		if(abs(p1-p2)==1&&(p2>0&&p1>0&&p2<6)) {			 
			int *out=new int[4];
			out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 			 
			return out;					 		
		}	
		else if(tt==0) { //only one difference
			rr=r;	
			continue;
		}			
		//only two difference		
		else if(abs(p1-p2)==2&&(p2>0&&p1>0&&p2<6)) {													 
			int *out=new int[4];
			out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 
			return out;					 
		} 
		else if(tt==1) { //only two difference
			rr=r;	
			continue;
		}
		else if(abs(p1-p2)==1&&(p2>0&&p1>0&&p2<6)) {
						
		}
		rr=r;		 
	}
	if(rr==0) return 0;
	if(tt<5) {
		tt++;
		goto re200;
	}
 	
	return 0;
}


int *Mutate::findhelixtoworkon() {
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);
	
	Res *r;
	
	//int m1=0,m2=0,m3=0;

	//float d=10000;

	if(TRES.logg)mpdb->write("ss");

	//tt=0: 1 or 2 residue difference, no secondary structure frame shift,
	//	less than 7 residues
	//tt=1: 1 or 2 residue difference, no secondary structure frame shift,
	//	more than 7 residues, stepwise
	//tt=2: 1 or 2 residue difference, zero occurs
	//tt=3: more than 3 residue, consider frame shift, less than 7 residue
	//tt=4: more than 3 residue, consider frame shift, more than 7 residue, stepwise
	//tt=5: deletion/insertion occur in secondary structure
	//tt=6: occur in head and tail.

	//int tt=0;
	//int itern=0;
	//re200:
	//int dis=1000;
	//m1=-1;m2=-1;m3=-1;itern=0;
	for(r=mpdb->chn->res;r;r=r->next) {

		if(r->next==0) continue;
		int n=compare[r->id0];	
		
		if(dssp[n]!='h') continue;

		if(mpdb->chn->islinked(r,r->next)==1)continue;	
 
 		//itern++;

		int n1,n2;
		int p1=0,p2=0;

		//the following code to detect if zig-zag alignment exists
		int window;
		for(window=1;window<5;window++) {
				
			p1=0;p2=0;
			for(n1=n-window;n1<=n+window+1;n1++) {
		 		if(owner->seqngap[n1]!='-') p1++;
				if(sqnto[n1]!='-')p2++; 
			}
			if(p1==p2) {
				int *out=new int[4];
				out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 				 
				return out;
			}		 
		}
		
		//end
		
		//find nearest secondary structure
	 	//p1=0;		
		for(n1=n-1;n1>=0;n1--) {
			//if(dssp[n1]=='h'&&owner->seqngap[n1]!='h') p1++;
			if(dssp[n1]=='h') continue;
			if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;									 
		}
		
		//p2=0;
		for(n2=n+1;n2<nlen;n2++) {
			//if(dssp[n2]=='h'&&owner->seqngap[n2]!='h') p2++;
			if(dssp[n2]=='h') continue;		 
			if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;	 		 	 	
		} 
		
		//end
 
		//not exist, tails
		if(n1==-1) continue;		
		if(n2==nlen) continue;		
			
		//detect the number of sequence in the region before and after mutation	
		p1=0;p2=0;
		
		int i;
		for(i=n1+1;i<n2;i++) {
			if(owner->seqngap[i]!='-') p1++;
			if(sqnto[i]!='-') p2++;
		}
		//end

		//simple case,zig and zag
		if(p1==p2) {
			int *out=new int[4];
			out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 
			return out;
		}		
		//end
		else {
			int *out=new int[4];
			out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 			 
			return out;	
		}		
	}
 
	return 0;
}


int *Mutate::findloopselftoworkon() {
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);
	Res *r;
	//int m1=0,m2=0,m3=0;

	//float d=10000;

	//mpdb->write("ss");

	//tt=0: 1 or 2 residue difference, no secondary structure frame shift,
	//	less than 7 residues
	//tt=1: 1 or 2 residue difference, no secondary structure frame shift,
	//	more than 7 residues, stepwise
	//tt=2: 1 or 2 residue difference, zero occurs
	//tt=3: more than 3 residue, consider frame shift, less than 7 residue
	//tt=4: more than 3 residue, consider frame shift, more than 7 residue, stepwise
	//tt=5: deletion/insertion occur in secondary structure
	//tt=6: occur in head and tail.

	//int tt=0;
	//int itern=0;
	//re200:
	//int dis=1000;
	//m1=-1;m2=-1;m3=-1;itern=0;
	for(r=mpdb->chn->res;r;r=r->next) {

		if(r->next==0) continue;
		int n=compare[r->id0];	
		int nn=compare[r->next->id0];
		if(dssp[n]!='-'&&dssp[nn]!='-') continue;

		if(mpdb->chn->islinked(r,r->next)==1)continue;	
 
 		//itern++;

		int n1,n2;
		int p1=0,p2=0;

		//the following code to detect if zig-zag alignment exists
		int window;
		for(window=1;window<5;window++) {
				
			p1=0;p2=0;
			for(n1=n-window;n1<=n+window+1;n1++) {
		 		if(owner->seqngap[n1]!='-') p1++;
				if(sqnto[n1]!='-')p2++; 
			}
			if(p1==p2) {
				int *out=new int[4];
				out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 				 
				return out;
			}		 
		}
		//end
		
		//find nearest secondary structure
				
		int t1=n-1;
		int t2=n+1;	  	

		while(1) {

			for(n1=t1;n1>=0;n1--) {
				if(dssp[n1]!='-') break;			 
				if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;									 
			}
		
			for(n2=t2;n2<nlen;n2++) {
				if(dssp[n2]!='-') break;
				if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;	 		 	 	
			}
		 	
			//end

			//not exist, tails
			if(n1==-1) continue;		
			if(n2==nlen) continue;		
			
			//detect the number of sequence in the region before and after mutation	
			p1=0;p2=0;
		
			int i;
			for(i=n1+1;i<n2;i++) {
				if(i<0||i>=nlen) continue;
				if(owner->seqngap[i]!='-') p1++;
				if(sqnto[i]!='-') p2++;
			}
			//end

			//simple case,zig and zag
			if(p1==p2) {
				int *out=new int[4];
				out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 
				return out;
			}					
			//end		 
			//only one difference
			if(abs(p1-p2)==1&&(p2>0&&p1>0&&p2<6)) {			 
				int *out=new int[4];
				out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 			 
				return out;					 		
			}		
			//only two difference		
			else if(abs(p1-p2)==2&&(p2>0&&p1>0&&p2<6)) {													 
				int *out=new int[4];
				out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 
				return out;					 
			} 
			if(p1>6&&p2>6) {
				break;
			}			
			else if(n1<=0&&n2>=nlen) {
				break;
			}
			else {
				if(t1>0) n1=t1-1;
				else 	 n1=t1;
				if(t2<nlen) n2=t2+1;
				else     n2=t2;
				continue;
			}
		}
		 
	}
 
	return 0;
}



void Mutate::setreliable() {
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;
	
	int nlen=strlen(sqnto);
 	
	if(rely) delete [] rely;
	rely=new int[nlen];
	
	int i;
	for(i=0;i<nlen;i++)  rely[i]=0;
 
	Res *r;
	for(r=mpdb->chn->res;r;r=r->next) {
		i=compare[r->id0];
		if(se->seqngap[i]!='-'&&owner->seqngap[i]!='-')rely[r->id0]=1; 
		else  rely[r->id0]=0;
	}
}

void Mutate::setnewunchanged() {
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;
	
	int nlen=strlen(sqnto);
 	
	if(rely) delete [] rely;
	rely=new int[nlen];
	
	int i;
	for(i=0;i<nlen;i++)  rely[i]=0;
 
	Res *r;
	for(r=pdbold->chn->res;r;r=r->next) {
		i=compare[r->id0];
		int j=owner->match[i];
		if(j==-1) continue;
		Res *r0=owner->resn[j];
		if(r0==0) continue;
		float d=r->directrmsdanyway(r0,0,3);
		if(d<0.5)   rely[r->id0]=2;
		else 	  rely[r->id0]=1;		
	}
}

void Mutate::setunchanged() {

        StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

        int nlen=strlen(sqnto);

        if(rely) delete [] rely;
        rely=new int[nlen];

        int i;
        for(i=0;i<nlen;i++)  rely[i]=1;

        Res *r;
        for(r=pdbold->chn->res;r;r=r->next) {
                i=compare[r->id0];
                int j=owner->match[i];
                if(j==-1) continue;
                Res *r0=owner->resn[j];
                if(r0==0) continue;
                float d=r->directrmsdanyway(r0,0,3);
                if(d<0.1) rely[r->id0]=0;
        }
}



void Mutate::setreliable(Res **tmp) {
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;
	
	int nlen=strlen(sqnto);
  
	if(rely==0) rely=new int[nlen];
	
	int i;
	for(i=0;i<nlen;i++)  rely[i]=0;
 
	Res *r;
	for(r=mpdb->chn->res;r;r=r->next) {
		i=compare[r->id0];
		if(se->seqngap[i]!='-'&&owner->seqngap[i]!='-')rely[r->id0]=2; 
		else  rely[r->id0]=1;
	}

	int n=0;
	while(tmp&&tmp[n]) {
		rely[tmp[n]->id0]=0;
		n++;
	}
}

void Mutate::printoutdistance(){

        int i;

        float mid,low,high;
        float x=1/segen->bound->npower;
        for(i=0;i<segen->bound->ndist;i++) {
                if(segen->bound->npower==1) {
                        mid=segen->bound->dist[3*i];
                        low=segen->bound->dist[3*i+1];
                        high=segen->bound->dist[3*i+2];
                }
                else {
                        mid=pow(segen->bound->dist[3*i],x);
                        low=pow(segen->bound->dist[3*i+1],x);
                        high=pow(segen->bound->dist[3*i+2],x);
                }
                int n1=segen->bound->atmpair[2*i];
                int n2=segen->bound->atmpair[2*i+1];
                float d=TRES.distance(segen->bound->atoms[n1],segen->bound->atoms[n2]);
		int e1=segen->bound->atoms[n1]->res->id0;
		int e2=segen->bound->atoms[n2]->res->id0;
                cerr<<segen->bound->atoms[n1]->id0<<" "<<rely[e1]<<" "<<segen->bound->atoms[n2]->id0<<" "<<rely[e2]<<" bound:"<<d<<" "<<mid<<" "<<low<<" "<<high<<" "<<endl;
        }
}



void Mutate::fixgappedstructure() {
 
	Segen segen;
	
  	segen.pdb=mpdb;
	segen.cid=mpdb->chn->id;
	 

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;
	
	int nlen=strlen(sqnto);
	mpdb->setlastresorder();
	
	int i;
	if(rely) delete [] rely;
	rely=new int[nlen];
	for(i=0;i<nlen;i++) rely[i]=0;

	Res *r,*r1,*r2;
	for(r=mpdb->chn->res;r;r=r->next) {
		i=compare[r->id0];	
		char c1,c2;
		c1=owner->seqngap[i];c2=se->seqngap[i]; 
		if(c1!='-'&&c2!='-'&&c1==c2) rely[i]=2;//same
		else if(c1!='-'&&c2!='-') rely[i]=1; //mutate
		else if(c1=='-'&&c2=='-') rely[i]=0; //same empty
		else if(c1=='-'&&c2!='-') rely[i]=-1; //insert
		else if(c1!='-'&&c2=='-') rely[i]=-2; //delete
		if(rely[i]>0) r->setflg(1);
		else 	      r->setflg(0);
	}

	
	float *desc=new float[4*nlen];

	char *sqa=new char[nlen];
	char *sqb=new char[nlen];
	int *ida=new int[nlen];
	int *idb=new int[nlen];
	sqa[0]='\0';sqb[0]='\0';
	char secd='-';
	int  edge=0;
 	int j1,j2,p1,p2,n,n1,n2;
	for(r=mpdb->chn->res;r;r=r->next) {
		
		if(r->next==0) continue;	
		if(mpdb->chn->islinked(r,r->next)==1)continue;	 
		//if(r->id0==12) mpdb->chn->write("s3");
		
		//find the two fixed points of n1 and n2
 		
		n=compare[r->id0];

 		if(dssp[n]!=secd) continue;
  
		//find the end segment position

		for(n1=n-1;n1>=0;n1--) {				
			if(dssp[n1]!=secd&&rely[n1]>0) break;
		}
		
		for(n2=n+1;n2<nlen;n2++) {
			if(dssp[n2]!=secd&&rely[n2]>0) break;
		}
			
		if(n1==-1&&edge==0) break;
		if(n2==nlen&&edge==0) break;
			
		p1=0;
		p2=0;
	 		
		for(n=n1+1;n<n2;n++) {
			if(owner->seqngap[n]!='-') {
				sqa[p1]=owner->seqngap[n];	
				ida[p1]=n;
				p1++;
			}		
			if(sqnto[n]!='-') {
				sqb[p2]=sqnto[n];
				idb[p2]=n;
				p2++;
			}
		}
		//sqa: original sequence
		//sqb: existing sequence

		//num=0;
				
		if(p1==p2&&p2) { //the sequence just equals, then do only mutation;
			mpdb->chn->setflgr(-9999);
			for(n=0;n<p1;n++) {
				j1=ida[n];
				j2=idb[n];
				j1=owner->match[j1];
				j2=match[j2];
				r1=resn[j1];
				r2=resn[j2];
				if(sqa[n]==sqb[n]) {						
					r2->transfer(r1);
				}
				else {						
					r2->transferbackbone(r1);
					r2->hooksidechain();
				}	
				r2->flag=10000;	
				r2->setflg(1);			
			}			 
	 		Scap scprd; 
			//modified 08/16/2002 for charge energy
			strcpy(scprd.force,"124");
			TRES.setdonaldcharge();
			scprd.dielectric=15;
			//
			scprd.pdb=mpdb;
			TRES.smt=TRES.nsmt-1;
			scprd.mutate();	
			//restore charge
			TRES.switchcharge(1);
			//			 	 	
			continue;							 				 
		}
		if(p1==0&&p2) {//insert, there is no between and insert into it.
			//decides who needs this inserted residues
			
			if(dssp[n]=='-'&&dssp[n1]=='h'&&dssp[n2]=='h') {
				
			}
			else if(dssp[n]=='-'&&dssp[n1]=='e'&&dssp[n2]=='h') {
				
			}
			else if(dssp[n]=='-'&&dssp[n1]=='h'&&dssp[n2]=='e') {
				
			}
			else if(dssp[n]=='-'&&dssp[n1]=='e'&&dssp[n2]=='e') {
				
			}
			else if(dssp[n]=='e'&&dssp[n1]=='h'&&dssp[n2]=='h') {

			}
			else if(dssp[n]=='e'&&dssp[n1]=='e'&&dssp[n2]=='h') {

			}
			else if(dssp[n]=='e'&&dssp[n1]=='h'&&dssp[n2]=='e') {

			}
			else if(dssp[n]=='e'&&dssp[n1]=='e'&&dssp[n2]=='e') {

			}
			else if(dssp[n]=='h'&&dssp[n1]=='h'&&dssp[n2]=='h') {

			}
			else if(dssp[n]=='h'&&dssp[n1]=='e'&&dssp[n2]=='h') {

			}
			else if(dssp[n]=='h'&&dssp[n1]=='h'&&dssp[n2]=='e') {

			}
			else if(dssp[n]=='h'&&dssp[n1]=='e'&&dssp[n2]=='e') {

			}


				
			int ii=0;
			for(n=n1;n<=n2;n++) {
				if(owner->seqngap[n]!='-') {
					int jj=owner->match[n];
					Res *rr=owner->resn[jj];
					Atm *ca=rr->isatm(" CA ");
					Atm *c=rr->isatm(" C  ");						
					 desc[4*ii]= ca->chi;
					desc[4*ii+1]= 10;
					desc[4*ii+2]= c->chi;
					desc[4*ii+3]= 10;
					ii++;
				}						 							
			}	
			if(n2-n1+1!=ii) {
				cerr<<"strange...not equal!"<<endl;
			}				 				 					 				
			 		 
		}
		else if(p1&&p2==0) {//delete
				 
		}
		else if(p1&&p2) {
 
			 	
		}
		else {
			//impossible
			cerr<<"impossible case happened.."<<endl;
		}
		 
		 
						
		 	
		//int jj=segen.fixsegment(n1,n2,desc);
 
		//if(jj) break;
 		 
		
	
		/*
		//performing segment minimization

		for(n1=n-1;n1>=0;n1--) {
			if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;			
		}		
		
		n=compare[r->next->id0];
		for(n2=n+1;n2<nlen;n2++) {
			if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;
		}
			
		//find the number of residues between n1 and n2, not including n1 and n2	
		
 
		while(1) {
			int i1,i2;
			
			//find boundary
			for(i1=n1-1;i1>=0;i1--) {
				if(owner->seqngap[i1]!='-'&&se->seqngap[i1]!='-') break;
			} 
			for(i2=n2+1;i2<nlen;i2++) {
				if(owner->seqngap[i2]!='-'&&se->seqngap[i2]!='-') break;
			}
			
			if(i1==-1) n1=-1;
			if(i2==nlen) n2=nlen;
			//if(n1==-1) n1=0;
			//if(n2==nlen) n2=nlen-1;

			//find the original and new sequence
			int p1=0;
			int p2=0; 
			for(n=n1+1;n<n2;n++) {
				if(owner->seqngap[n]!='-') p1++;
				if(sqnto[n]!='-') p2++;
			}

			if(n1==-1) n1=0;
			if(n2==nlen) n2=nlen-1;
			if(n1==0||n2==nlen-1) {
				goto re300;
			}
			else if(p1==0&&p2==0) {
				break;
			}
			else if(p1==0&&p2) {			
				goto re200; //insert
			}
			else if(p1&&p2==0) {
				goto re200; //delete 
			}
			
			re300:
			i1=match[n1];
			i2=match[n2];
			if(i1==-1) i1=mpdb->chn->res->id0;
			if(i2==-1) i2=mpdb->chn->lastres()->id0;			

			int out;
			if((p1&&p2)||n1==0||n2==nlen-1) out=segen.linkviaboth(i1,i2,r->id0);
  			else out=0;

			if(out==1) break;

			cerr<<"can not close residue:"<<resn[i1]->name<<resn[i1]->id;
			cerr<<" "<<resn[i2]->name<<resn[i2]->id<<endl;	
		
			re200:
						
			
			if(dssp[n1]!='-'&&dssp[n2]!='-') {
				int e1=0,e2=0;

				for(p1=n1-1;p1>=0;p1--) {
					if(dssp[p1]=='-') break;
					e1++;	
				}
				for(p2=n2+1;p2<nlen;p2++) {
					if(dssp[p2]=='-') break;
					e2++;	 	
				}	
				
				for(p1=n1-1;p1>=0;p1--) {
					if(owner->seqngap[p1]!='-'&&se->seqngap[p1]!='-') break;	
				}
				for(p2=n2+1;p2<nlen;p2++) {
					if(owner->seqngap[p2]!='-'&&se->seqngap[p2]!='-') break;	
				}	
				if(p1==-1) {
					n1=p1;
				}
				else if(p2==nlen) {
					n2=p2;
				}
				else if(dssp[n1]=='h'&&dssp[n1]=='e') {
					n1=p1;
				}
				else if(dssp[n1]=='e'&&dssp[n2]=='h') {
					n2=p2;
				}
				else if(e1>e2) {
					n1=p1;
				}
				else if(e2>e1) {
					n2=p2;
				}
				else  {
					n1=p1;
					n2=p2;
				}	
			}
			else if(dssp[n1]=='-'&&dssp[n2]!='-') {
				for(p1=n1-1;p1>=0;p1--) {
					if(owner->seqngap[p1]!='-'&&se->seqngap[p1]!='-') break;	
				}
				
				n1=p1;				 
			}
			else if(dssp[n1]!='-'&&dssp[n2]=='-') {
				
				for(p2=n2+1;p1<nlen;p2++) {
					if(owner->seqngap[p2]!='-'&&se->seqngap[p2]!='-') break;	
				}
				n2=p2;				 
			}  
			else {
				for(p1=n1-1;p1>=0;p1--) {
					if(owner->seqngap[p1]!='-'&&se->seqngap[p1]!='-') break;	
				}
				for(p2=n2+1;p1<nlen;p2++) {
					if(owner->seqngap[p2]!='-'&&se->seqngap[p2]!='-') break;	
				}
				if(dssp[p1]=='-'&&dssp[n2]=='-') {
					n1=p1;n2=p2;
				}
				else if(dssp[p1]!='-'&&dssp[n2]=='-') {
					n2=p2;
				}
				else if(dssp[p1]=='-'&&dssp[n2]!='-') {
					n1=p1;
				}
				else {
					n1=p1;n2=p2;
				}
			}  
			 
		}
		*/
 
	}	  
 
 
	if(TRES.logg)mpdb->chn->write("mutate5.pdb");
}
 

void Mutate::linkgappedstructure() {
 
	Segen segen;
	
  	segen.pdb=mpdb;
	segen.cid=mpdb->chn->id;
	 

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;
	
	int nlen=strlen(sqnto);
		
	mpdb->setlastresorder();

	Res *r;
	for(r=mpdb->chn->res;r;r=r->next) {
		
		if(r->next==0) continue;	
		if(mpdb->chn->islinked(r,r->next)==1)continue;	 
		if(TRES.logg)if(r->id0==12) mpdb->chn->write("s3");
		
		//find the two fixed points of n1 and n2
		int n1,n2;
		
		int n=compare[r->id0];
 
		for(n1=n-1;n1>=0;n1--) {
			if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;			
		}		
		
		n=compare[r->next->id0];
		for(n2=n+1;n2<nlen;n2++) {
			if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;
		}
			
		
		//find the number of residues between n1 and n2, not including n1 and n2	
		
 
		while(1) {
			int i1,i2;
			
			//find boundary
			for(i1=n1-1;i1>=0;i1--) {
				if(owner->seqngap[i1]!='-'&&se->seqngap[i1]!='-') break;
			} 
			for(i2=n2+1;i2<nlen;i2++) {
				if(owner->seqngap[i2]!='-'&&se->seqngap[i2]!='-') break;
			}
			
			if(i1==-1) n1=-1;
			if(i2==nlen) n2=nlen;
			//if(n1==-1) n1=0;
			//if(n2==nlen) n2=nlen-1;

			//find the original and new sequence
			int p1=0;
			int p2=0; 
			for(n=n1+1;n<n2;n++) {
				if(owner->seqngap[n]!='-') p1++;
				if(sqnto[n]!='-') p2++;
			}

			if(n1==-1) n1=0;
			if(n2==nlen) n2=nlen-1;
			if(n1==0||n2==nlen-1) {
				goto re300;
			}
			else if(p1==0&&p2==0) {
				break;
			}
			else if(p1==0&&p2) {			
				goto re200; //insert
			}
			else if(p1&&p2==0) {
				goto re200; //delete 
			}
			
			re300:
			i1=match[n1];
			i2=match[n2];
			if(i1==-1) i1=mpdb->chn->res->id0;
			if(i2==-1) i2=mpdb->chn->lastres()->id0;			

			int out;
			if((p1&&p2)||n1==0||n2==nlen-1) out=segen.linkviaboth(i1,i2,r->id0);
  			else out=0;

			if(out==1) break;

			cerr<<"can not close residue:"<<resn[i1]->name<<resn[i1]->id;
			cerr<<" "<<resn[i2]->name<<resn[i2]->id<<endl;	
		
			re200:
						
			
			if(dssp[n1]!='-'&&dssp[n2]!='-') {
				int e1=0,e2=0;

				for(p1=n1-1;p1>=0;p1--) {
					if(dssp[p1]=='-') break;
					e1++;	
				}
				for(p2=n2+1;p2<nlen;p2++) {
					if(dssp[p2]=='-') break;
					e2++;	 	
				}	
				
				for(p1=n1-1;p1>=0;p1--) {
					if(owner->seqngap[p1]!='-'&&se->seqngap[p1]!='-') break;	
				}
				for(p2=n2+1;p2<nlen;p2++) {
					if(owner->seqngap[p2]!='-'&&se->seqngap[p2]!='-') break;	
				}	
				if(p1==-1) {
					n1=p1;
				}
				else if(p2==nlen) {
					n2=p2;
				}
				else if(dssp[n1]=='h'&&dssp[n1]=='e') {
					n1=p1;
				}
				else if(dssp[n1]=='e'&&dssp[n2]=='h') {
					n2=p2;
				}
				else if(e1>e2) {
					n1=p1;
				}
				else if(e2>e1) {
					n2=p2;
				}
				else  {
					n1=p1;
					n2=p2;
				}	
			}
			else if(dssp[n1]=='-'&&dssp[n2]!='-') {
				for(p1=n1-1;p1>=0;p1--) {
					if(owner->seqngap[p1]!='-'&&se->seqngap[p1]!='-') break;	
				}
				
				n1=p1;				 
			}
			else if(dssp[n1]!='-'&&dssp[n2]=='-') {
				
				for(p2=n2+1;p1<nlen;p2++) {
					if(owner->seqngap[p2]!='-'&&se->seqngap[p2]!='-') break;	
				}
				n2=p2;				 
			}  
			else {
				for(p1=n1-1;p1>=0;p1--) {
					if(owner->seqngap[p1]!='-'&&se->seqngap[p1]!='-') break;	
				}
				for(p2=n2+1;p1<nlen;p2++) {
					if(owner->seqngap[p2]!='-'&&se->seqngap[p2]!='-') break;	
				}
				if(dssp[p1]=='-'&&dssp[n2]=='-') {
					n1=p1;n2=p2;
				}
				else if(dssp[p1]!='-'&&dssp[n2]=='-') {
					n2=p2;
				}
				else if(dssp[p1]=='-'&&dssp[n2]!='-') {
					n1=p1;
				}
				else {
					n1=p1;n2=p2;
				}
			}  
			 
		}
 
	}	  
 
 
	if(TRES.logg)mpdb->chn->write("mutate5.pdb");
}

 
void Mutate::linkinitialstructure() {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

	int nlen=strlen(se->seqngap);

	int *addon=new int[nlen];
	int i=0;
 
	for(i=0;i<nlen;i++) {
		if(se->seqngap[i]!='-'&&owner->seqngap[i]=='-')addon[i]=1; 
		else if(se->seqngap[i]!='-'&&owner->seqngap[i]!='-')addon[i]=0;
		else if(se->seqngap[i]=='-'&&owner->seqngap[i]!='-')addon[i]=-1; 
	}


	Rotate rot;
	int j,n1,n2;
	Res *r,*r0;
	 
	while(1) {
	   	int tot=0;
		for(i=0;i<nlen;i++) {
			//if(sqnto[i]=='-') continue;
			if(addon[i]!=1) continue;
			//if(se->seqngap[i]!='-'&&owner->seqngap[i]!='-')continue;
		 
			j=0;
			n1=0;
			n2=0;

			for(j=i-1;j>=0;j--) {			 
				//if(se->seqngap[j]!='-'&&owner->seqngap[j]!='-') {
				if(addon[j]==0) {
					n1++;
				}
				else {
					break;
				}
			}

			for(j=i+1;j<nlen;j++) {
				//if(se->seqngap[j]!='-'&&owner->seqngap[j]!='-') {
				if(addon[j]==0) {
					n2++;
				}
				else {
					break;
				}
			}

			if(n1==0&&n2==0) continue;
			r=resn[match[i]];
			if(n2>n1) {
				r0=r->next;
				rot.link(r,r0,r->id0,0);
				addon[i]=2;tot++;
			}
			else {
				r0=resn[match[i-1]];
				rot.link(r0,r,r->id0,1);
				addon[i]=2;tot++;
			} 
		}

		cerr<<"the total number of residues assembled:"<<tot<<endl;
		if(tot==0) {
			for(i=0;i<nlen;i++) {
				if(addon[i]==2) {
					addon[i]=0;
					tot++;
				}
			} 
		}		
		if(tot==0) break;
	}
	delete [] addon;addon=0;

	if(TRES.logg)mpdb->chn->write("mutate4.pdb");
}

void Mutate::createinitialstructure() {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

	int nlen=strlen(se->seqngap);
	int i;
	for(i=0;i<nlen;i++) {
		
		if(sqnto[i]=='-'&&se->seqngap[i]!='-') {//insert
			int n;
			for(n=i-1;n>=0;n--) {
				if(sqnto[n]!='-') break; 
			}
			int m;
			for(m=i+1;m<nlen;m++) {
				if(sqnto[m]!='-') break; 
			}
			Res *r0=0,*r=0;
			if(n!=-1) {			
		 		int j=match[n];
				r0=resn[j];
			}	
			if(m!=nlen) {
				int j=match[m];
				r=resn[j];
			}
			Tres *tres_temp;
			tres_temp=TRES[se->seqngap[i]];
			Res *rr=new Res(tres_temp);
			rr->nemp=9;
			rr->temp=new float[9];
			int h;
			for(h=0;h<9;h++) rr->temp[h]=tres_temp->head[h];
			if(r0&&r) {				
				r0->next=rr;
				rr->next=r;
				rr->chn=r0->chn;
				rr->id=r0->id+1;
				for(Res *t=r;t;t=t->next) t->id++;
			}	
			else if(r0) {				
				r0->next=rr;				 
				rr->chn=r0->chn;
				rr->id=r0->id+1;
			}		 
			else if(r) {				
				mpdb->chn->res=rr;
				rr->next=r;				 		 
				rr->chn=r->chn;
				rr->id=r->id;
				for(Res *t=r;t;t=t->next) t->id++;
			}
			sqnto[i]=se->seqngap[i];
			mpdb->configure();
			setmatch();
		}
		else if(sqnto[i]!='-'&&se->seqngap[i]=='-') {//delete
			int n,m;
			for(n=i-1;n>=0;n--) {
				if(sqnto[n]!='-') break; 
			}
			for(m=i+1;m<nlen;m++) {
				if(sqnto[m]!='-') break; 
			}
			Res *r0=0,*r=0,*r1=0;
			if(n!=-1) {			
		 		int j=match[n];
				r0=resn[j];
			}	
			if(m!=nlen) {
				int j=match[m];
				r=resn[j];
			}
			if( sqnto[i]!='-') {
				int j=match[i];
				r1=resn[j];
			}
			if(r0&&r) {				 
				r0->next=r;
				if(r1) {
					r1->next=0;
					delete r1;r1=0;
				}
				for(Res *t=r;t;t=t->next) t->id--;
			}			 
			else if(r0) {
				r0->next=0;
				if(r1) {
					r1->next=0;
					delete r1;r1=0;
				} 
			}		 
			else if(r) {
				mpdb->chn->res=r;				
				if(r1) {
					r1->next=0;
					delete r1;r1=0;
				} 
				for(Res *t=r;t;t=t->next) t->id--; 
			}
			sqnto[i]=se->seqngap[i];
			mpdb->configure();
			setmatch();	
		}
		
	}
	mpdb->configure();
	setmatch();	
	cerr<<dssp<<endl;
	cerr<<sqnto<<endl;
	cerr<<se->seqngap<<endl;
	if(TRES.logg)mpdb->chn->write("mutate3.pdb");
}

 


void Mutate::dosidechainmutation(int f) {
 
	//get the sequence class of this 

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

	int nlen=strlen(se->seqngap);

	if(TRES.logg)mpdb->chn->write("mutate0.pdb");

	int i;

	//float **xyz =new float[nlen];

	//for(i=0;i<nlen;i++) xyz[i]=0;

	mpdb->setflgr(-99999);
	int tot=0;
        for(i=0;i<nlen;i++) {
		if(sqnto[i]=='-'||se->seqngap[i]=='-') continue;
		if(se->seqngap[i]==sqnto[i])  continue;
		
		int n=owner->match[i];
		 
		Res *r=mpdb->chn->isres(n);

		if(r==0) {
			cerr<<"the residue does not exist:"<<n<<endl;
			continue;
		}

		float xyz[100];
		r->transfer(xyz,0);
		char c=se->seqngap[i];		
		r->mutateResidue(c);
		mpdb->configure();
		if(strchr("AGP",r->name)) {					
			cerr<<"residue:"<<sqnto[i]<<r->id<<" is mutated to:"<<c<<endl;
			sqnto[i]=c;
			tot++;	
			continue;
		}
		r->flag=10000; 
		Scap scprd;	 	 		
		//08/16/2002 for charge
		strcpy(scprd.force,"124");
		TRES.setdonaldcharge();
		scprd.dielectric=15;	
		//
		scprd.pdb=mpdb;
		TRES.smt=TRES.nsmt-1;
		float d=scprd.mutate();
		TRES.switchcharge(1);
		if(d>5&&f==0) { //no change at all
			r->mutateResidue(sqnto[i]);
			mpdb->configure();
			r->transfer(xyz,1);
		}
		else if(d>5&&f==1) { //too many clashes,delete sidechain 
			mpdb->chn->transform(r,r->id0,3);
			mpdb->configure();
		} 
		else {								 
			cerr<<"residue:"<<sqnto[i]<<r->id<<" is mutated to:"<<c<<" with energy:"<<d<<endl;
			sqnto[i]=c;	
			tot++; 
		}
		r->flag=-99999;
		r->setfreenextonmore();	
		if(r->more) delete r->more;r->more=0;				 				
	}
	setmatch();
 	cerr<<getsecondary()<<endl;
	cerr<<sqnto<<endl;
	cerr<<se->seqngap<<endl;
	if(TRES.logg)mpdb->chn->write("mutate.pdb");
	 
} 


void Mutate::dooutlierdeletion(){

	 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return; 

	int nlen=strlen(se->seqngap);

	int i=0;
 	int n=0;
	for(i=0;i<nlen;i++) {
		 
		if(sqnto[i]!='-'&&se->seqngap[i]=='-') {
			Res *r=mpdb->chn->res;
			mpdb->chn->res=r->next;
			r->chn=0;r->next=0;
			cerr<<"residue:"<<r->name<<r->id0<<" is deleted"<<endl;	
			delete r;r=0;
			sqnto[i]='-';
			n++;		
										
		}
		else if(sqnto[i]=='-'&&se->seqngap[i]=='-') {
			continue;
		}
		else if(sqnto[i]=='-'&&se->seqngap[i]!='-') {
			break;
		}
		else if(sqnto[i]!='-'&&se->seqngap[i]!='-') {
			break;
		}
	}


	for(i=nlen-1;i>=0;i--){

		if(sqnto[i]!='-'&&se->seqngap[i]=='-') {
			int n=owner->match[i];
			Res *r=mpdb->isres(n-1);
			Res *rr=r->next;
			r->next=0;
			rr->next=0;
			cerr<<"residue:"<<rr->name<<rr->id0<<" is deleted"<<endl;	
			delete rr;	rr=0;
			sqnto[i]='-';		
			n++;			
		}
		else if(sqnto[i]=='-'&&se->seqngap[i]=='-') {
			continue;
		}
		else if(sqnto[i]=='-'&&se->seqngap[i]!='-') {
			break;
		}
		else if(sqnto[i]!='-'&&se->seqngap[i]!='-') {
			break;
		} 		 
	}

	mpdb->configure();

	setmatch();
 	if(TRES.logg)mpdb->write("mutate1.pdb");
	cerr<<sqnto<<endl;
	cerr<<se->seqngap<<endl;

}

void Mutate::doeditloop(int len){

	
	//in the window, the new sequence is window+/-diff

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	int nlen=strlen(se->seqngap);
	
	int i=0;
 	char sqa[100],sqb[100];
 
	for(i=0;i<nlen;i++) {
		
		if(dssp[i]=='e'||dssp[i]=='h') continue; //into loops

		//find loop stems
		int n1,n2,j;

		for(j=i;j>=0;j--) { 
			if(dssp[j]=='e'||dssp[j]=='h') break;
		}
		n1=j;

		for(j=i;j<nlen;j++) { 
			if(dssp[j]=='e'||dssp[j]=='h') break;
		}
		n2=j;
 
		//find original and new sequence of loops
		int m=0;
		for(j=n1+1;j<n2;j++) {
			if(sqnto[j]=='-') continue;
			sqa[m++]=sqnto[j];	
		}
		sqa[m]='\0';

		m=0;
		for(j=n1+1;j<n2;j++) {
			if(se->seqngap[j]=='-') continue;
			sqb[m++]=se->seqngap[j];	
		}
		sqb[m]='\0';

		//skip too long loops
		if(strlen(sqb)>len)  { i=n2-1; continue; }

		//same loop length, do only mutation
		if(strlen(sqa)==strlen(sqb)) {
			//edit the alignment
			for(j=n1+1;j<n2;j++) sqnto[j]='-';
			m=0;
			for(j=n1+1;j<n2;j++) {
				if(se->seqngap[j]=='-') continue;
				sqnto[j]=sqa[m++];
			}
			setmatch();
			dosidechainmutation(1);
		}
		
		else if(strlen(sqa)!=0) {
			//int m1=-1;
			//int m2=-1;
			
			mpdb->chn->replace(sqb,n1,n2);
			for(j=n1+1;j<n2;j++) {
				sqnto[j]=se->seqngap[j];
			}
			mpdb->configure();
			setmatch();
			if(TRES.logg)mpdb->write("mutate1.pdb");
			cerr<<dssp<<endl;
        		cerr<<sqnto<<endl;
        		cerr<<se->seqngap<<endl;

		}
	}

	mpdb->configure();
	setmatch();
 	if(TRES.logg)mpdb->write("mutate1.pdb");
	cerr<<sqnto<<endl;
	cerr<<se->seqngap<<endl;
} 

void Mutate::doeditresidue(int fff){

	 
} 


void Mutate::dooutlierinsertion(){

	char outer[1000];

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	int nlen=strlen(se->seqngap);

	int i=0;
 	int n=0;
	outer[0]='\0';
	int m=0;
	for(i=0;i<nlen;i++) {
		 
		if(sqnto[i]=='-'&&se->seqngap[i]!='-') {
			
			outer[m++]=se->seqngap[i];n++;
			sqnto[i]=se->seqngap[i];															
		}
		else if(sqnto[i]=='-'&&se->seqngap[i]=='-') {
			continue;
		}
		else if(sqnto[i]!='-'&&se->seqngap[i]=='-') {
			break;
		}
		else if(sqnto[i]!='-'&&se->seqngap[i]!='-') {
			break;
		}
	}

	outer[m]='\0';
	//build the header sequence
	inserthead(outer);
	
	//setmatch();
	//mpdb->configure();	
	outer[0]='\0';m=0;
	for(i=nlen-1;i>=0;i--){

		if(sqnto[i]=='-'&&se->seqngap[i]!='-') {
			outer[m++]=se->seqngap[i];n++;
			sqnto[i]=se->seqngap[i];		
		}
		else if(sqnto[i]=='-'&&se->seqngap[i]=='-') {
			continue;
		}
		else if(sqnto[i]!='-'&&se->seqngap[i]=='-') {
			break;
		}
		else if(sqnto[i]!='-'&&se->seqngap[i]!='-') {
			break;
		} 		 
	}
	outer[m]='\0';
	appendtail(outer);
	
	//mpdb->configure();
	//setmatch();
 	if(TRES.logg)mpdb->write("mutate2.pdb");
	cerr<<dssp<<endl;
	cerr<<sqnto<<endl;
	cerr<<se->seqngap<<endl;
}

void Mutate::inserthead(char *outer){

	if(strlen(outer)==0) return;
	cerr<<"insert sequence at head:"<<outer<<endl;	
	mpdb->chn->insertresidues(outer,mpdb->chn->res->id0,-1);
	
	setmatch();
}

void Mutate::appendtail(char *outer){

	if(strlen(outer)==0) return;
	cerr<<"append sequence at head:"<<outer<<endl;
	Res *r=mpdb->chn->lastres();
	mpdb->chn->insertresidues(outer,r->id0,1);

	setmatch();
}


void Mutate::setmatch() {
//reset the compare table from alignment to structure

        if(sqnto==0) return;
 
        int n=strlen(sqnto);
 
        int i;

        for(i=0;i<n;i++) {
		match[i]=-1;
 		compare[i]=-1;
		resn[i]=0;
		seqn[i]='\0';
	}

        int m=0;

        for(i=0;i<n;i++) {
                if(sqnto[i]=='-') continue;                
                match[i]=m;
		compare[m]=i;
             	resn[m]=mpdb->chn->isres0(m);
		seqn[m]=sqnto[i];
                m++;
        }
 	seqn[m]=0;
	resn[m]=0;
	
        return;
}


char * Mutate::getsecondary() {
    
	int n=strlen(sqnto);
	char *out=new char[n+2];
    
	for(int i=0;i<n;i++) {
		if(sqnto[i]=='-') {
			out[i]=' ';
			continue;
		}
		int n=match[i];
		Res *r=mpdb->isres(n);
		if(r==0) {
			out[i]=' ';
			continue;
		}
		out[i]=r->sec;
	}
	out[n]='\0';
	return out;
}
void Mutate::setdssp() {
	setdssp(0);
}
void Mutate::setdssp(int ff) {

 	
	StrFmt *parent=owner->getparentStrFmt();
	StrFmt *str=parent->findstructurefmt();
	StrFmt *se=parent->findsequencefmt();

	int n=strlen(sqnto);

	if(dssp) delete [] dssp; dssp=0;
	dssp=new char[n+1];

	int i;

	for(i=0;i<n;i++) dssp[i]='?';
	dssp[n]='\0'; 	

	str->pdb->setthreestatesec();

	for(i=0;i<n;i++) {
		if(ff==0) {
			int j=str->match[i];
			if(j==-1) continue;
			Res *r=str->resn[j];
			if(r==0) continue;		 
			dssp[i]=r->sec;
		}
		else {
			int j=match[i];
			if(j==-1) continue;
			Res *r=resn[j];
			if(r==0) continue;		 
			dssp[i]=r->sec;
		}
	}

	for(i=0;i<n;i++) {
		if(dssp[i]!='?') continue;
		char a,b;
		int m;
		for(m=i-1;m>=0;m--) {        
			if(dssp[m]!='?') {
				a=dssp[m];break;
			}
                }  
		for(m=i+1;m<n;m++) {
			if(dssp[m]!='?') {
                               b=dssp[m];break;
			}
		}	
		if(a=='?'||b=='?') dssp[i]='-';
		else if(a!=b) dssp[i]='-';
		else dssp[i]=a;
	}
	dssp[n]='\0';
	
	if(TRES.logg>3) cerr<<dssp<<endl;
	if(TRES.logg>3) cerr<<sqnto<<endl;
	if(TRES.logg>3) cerr<<se->seqngap<<endl;
	 
}

int *Mutate::findbadsegment(int n,int tt) {
 		
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;
	int nlen=strlen(sqnto);
 
	int n1,n2,m;
	int p1=0,p2=0;
	
	if(tt==1) {
		
		if(dssp[n]!='-') return 0;
		
		int window;
		for(window=1;window<5;window++) {				
			p1=0;p2=0;
			for(m=n-window;m<=n+window+1;m++) {
		 		if(owner->seqngap[m]!='-') p1++;
				if(se->seqngap[m]!='-')p2++; 
			}
			if(p1==p2) {
				int *out=new int[4];
				out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 				 
				return out;
			}		 
		}		

		//find nearest secondary structure	  
		for(n1=n-1;n1>=0;n1--) {
			if(dssp[n1]=='-') continue;
			if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;									 
		}
		
		for(n2=n+1;n2<nlen;n2++) {
			if(dssp[n2]=='-') continue;
			if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;	 		 	 	
		}	 
		//end

		//not exist, tails
		if(n1==-1) return 0;		
		if(n2==nlen) return 0;		
			
		//detect the number of sequence in the region before and after mutation	
		p1=0;p2=0;
		
		int i;
		for(i=n1+1;i<n2;i++) {
			if(owner->seqngap[i]!='-') p1++;
			if(se->seqngap[i]!='-') p2++;
		}
		//end
	 
		//simple case,zig and zag
		if(p1==p2) {
			int *out=new int[4];
			out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 
			return out;
		}		
		//end
 
		//only one difference
		if(abs(p1-p2)==1&&(p2>0&&p1>0&&p2<6)) {			 
			int *out=new int[4];
			out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 			 
			return out;					 		
		}		 
	}
	else if(tt==2) {
		
		if(dssp[n]!='-') return 0;
		
		int window;
		for(window=1;window<5;window++) {				
			p1=0;p2=0;
			for(m=n-window;m<=n+window+1;m++) {
		 		if(owner->seqngap[m]!='-') p1++;
				if(se->seqngap[m]!='-')p2++; 
			}
			if(p1==p2) {
				int *out=new int[4];
				out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 				 
				return out;
			}		 
		}		

		//find nearest secondary structure	  
		for(n1=n-1;n1>=0;n1--) {
			if(dssp[n1]=='-') continue;
			if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;									 
		}
		
		for(n2=n+1;n2<nlen;n2++) {
			if(dssp[n2]=='-') continue;
			if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;	 		 	 	
		}	 
		//end

		//not exist, tails
		if(n1==-1) return 0;		
		if(n2==nlen) return 0;		
			
		//detect the number of sequence in the region before and after mutation	
		p1=0;p2=0;
		
		int i;
		for(i=n1+1;i<n2;i++) {
			if(owner->seqngap[i]!='-') p1++;
			if(se->seqngap[i]!='-') p2++;
		}
		//end
	 
		//simple case,zig and zag
		if(p1==p2) {
			int *out=new int[4];
			out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 
			return out;
		}		
		//end
  
		//only two difference		
		if(abs(p1-p2)==2&&(p2>0&&p1>0&&p2<6)) {													 
			int *out=new int[4];
			out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 
			return out;					 
		}  
	}
	else if(tt==3) {
		
		if(dssp[n]!='-') return 0;
		
		int window;
		for(window=1;window<5;window++) {				
			p1=0;p2=0;
			for(m=n-window;m<=n+window+1;m++) {
		 		if(owner->seqngap[m]!='-') p1++;
				if(se->seqngap[m]!='-')p2++; 
			}
			if(p1==p2) {
				int *out=new int[4];
				out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 				 
				return out;
			}		 
		}		

		//find nearest secondary structure	  
		for(n1=n-1;n1>=0;n1--) {
			if(dssp[n1]=='-') continue;
			if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;									 
		}
		
		for(n2=n+1;n2<nlen;n2++) {
			if(dssp[n2]=='-') continue;
			if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;	 		 	 	
		}	 
		//end

		//not exist, tails
		if(n1==-1) return 0;		
		if(n2==nlen) return 0;		
			
		//detect the number of sequence in the region before and after mutation	
		p1=0;p2=0;
		
		int i;
		for(i=n1+1;i<n2;i++) {
			if(owner->seqngap[i]!='-') p1++;
			if(se->seqngap[i]!='-') p2++;
		}
		//end
	 
		//simple case,zig and zag
		if(p1==p2) {
			int *out=new int[4];
			out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 
			return out;
		}		
		//end
  
		//only two difference		
		if(abs(p1-p2)==2&&(p2>0&&p1>0&&p2<6)) {													 
			int *out=new int[4];
			out[0]=n1;out[1]=n2;out[2]=n;out[3]=0; 
			return out;					 
		}  
	} 
	 
	return 0;
}
void Mutate::createinitialstructure(int mm) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

	int nlen=strlen(se->seqngap);

	
	int m1,m2;

	for(m1=mm-1;m1>=0;m1--) {
		//if(owner->seqngap[m1]!='-'&&se->seqngap[m1]!='-') break;
		if(sqnto[m1]!='-'&&se->seqngap[m1]!='-') break;									 
	}
		
	for(m2=mm+1;m2<nlen;m2++) {
		//if(owner->seqngap[m2]!='-'&&se->seqngap[m2]!='-') break;
		if(sqnto[m2]!='-'&&se->seqngap[m2]!='-') break;	 		 	 	
	}	 
	
	Res *ga,*gb;

	if(m1==-1)   ga=mpdb->chn->res;
	else	     ga=resn[match[m1]];
	
	if(m2==nlen) gb=mpdb->chn->lastres();
	else	     gb=resn[match[m2]];

	int i;
	int tt=0;
	int jj;
	for(jj=m1+1;jj<m2;jj++) {

		if(tt>1) break;
		
		if(se->seqngap[jj]!='-'&&owner->seqngap[jj]=='-'&&m2-m1-1>1) { //insert
			tt++;	
			int ii=findloose(ga,gb,mm);
			if(ii==1) i=m2-1;		 
			else  	  i=m1+1;		
		}
		else {  //delete
			i=jj;
		}

		if(sqnto[i]=='-'&&se->seqngap[i]!='-') {//insert
			int n;
			for(n=i-1;n>=0;n--) {
				if(sqnto[n]!='-') break; 
			}
			int m;
			for(m=i+1;m<nlen;m++) {
				if(sqnto[m]!='-') break; 
			}
			Res *r0=0,*r=0;
			if(n!=-1) {			
		 		int j=match[n];
				r0=resn[j];
			}	
			if(m!=nlen) {
				int j=match[m];
				r=resn[j];
			}
			Tres *tres_temp;
			tres_temp=TRES[se->seqngap[i]];
			Res *rr=new Res(tres_temp);
			rr->nemp=9;
			rr->temp=new float[9];
			int h;
			for(h=0;h<9;h++) rr->temp[h]=tres_temp->head[h];
			if(r0&&r) {				
				r0->next=rr;
				rr->next=r;
				rr->chn=r0->chn;
				rr->id=r0->id+1;
				for(Res *t=r;t;t=t->next) t->id++;
			}	
			else if(r0) {				
				r0->next=rr;				 
				rr->chn=r0->chn;
				rr->id=r0->id+1;
			}		 
			else if(r) {				
				mpdb->chn->res=rr;
				rr->next=r;				 		 
				rr->chn=r->chn;
				rr->id=r->id;
				for(Res *t=r;t;t=t->next) t->id++;
			}
			sqnto[i]=se->seqngap[i];
			mpdb->configure();
			setmatch();
		}
		else if(sqnto[i]!='-'&&se->seqngap[i]=='-') {//delete
			int n,m;
			for(n=i-1;n>=0;n--) {
				if(sqnto[n]!='-') break; 
			}
			for(m=i+1;m<nlen;m++) {
				if(sqnto[m]!='-') break; 
			}
			Res *r0=0,*r=0,*r1=0;
			if(n!=-1) {			
		 		int j=match[n];
				r0=resn[j];
			}	
			if(m!=nlen) {
				int j=match[m];
				r=resn[j];
			}
			if( sqnto[i]!='-') {
				int j=match[i];
				r1=resn[j];
			}
			if(r0&&r) {				 
				r0->next=r;
				if(r1) {
					r1->next=0;
					delete r1;r1=0;
				}
				for(Res *t=r;t;t=t->next) t->id--;
			}			 
			else if(r0) {
				r0->next=0;
				if(r1) {
					r1->next=0;
					delete r1;r1=0;
				} 
			}		 
			else if(r) {
				mpdb->chn->res=r;				
				if(r1) {
					r1->next=0;
					delete r1;r1=0;
				} 
				for(Res *t=r;t;t=t->next) t->id--; 
			}
			sqnto[i]=se->seqngap[i];
			mpdb->configure();
			setmatch();	
		}
		
	}
	mpdb->configure();
	setmatch();	
	cerr<<dssp<<endl;
	cerr<<sqnto<<endl;
	cerr<<se->seqngap<<endl;
	if(TRES.logg)mpdb->chn->write("mutate3.pdb");
}
void Mutate::linkinitialstructure(int mm) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return ;

	int nlen=strlen(se->seqngap);

	int *addon=new int[nlen];
	int i=0;
 
	

	int m1,m2;

	for(m1=mm-1;m1>=0;m1--) {
		if(owner->seqngap[m1]!='-'&&se->seqngap[m1]!='-') break;									 
	}
		
	for(m2=mm+1;m2<nlen;m2++) {
		if(owner->seqngap[m2]!='-'&&se->seqngap[m2]!='-') break;	 		 	 	
	}	 
	
	for(i=0;i<nlen;i++) {
		if(se->seqngap[i]!='-'&&owner->seqngap[i]=='-')addon[i]=1; 
		else if(se->seqngap[i]!='-'&&owner->seqngap[i]!='-')addon[i]=0;
		else if(se->seqngap[i]=='-'&&owner->seqngap[i]!='-')addon[i]=-1; 
	}

	Rotate rot;
	int j,n1,n2;
	Res *r,*r0;
	 
	while(1) {
	   	int tot=0;
		for(i=m1+1;i<m2;i++) {
			//if(sqnto[i]=='-') continue;
			if(addon[i]!=1) continue;
			//if(se->seqngap[i]!='-'&&owner->seqngap[i]!='-')continue;
		 
			j=0;
			n1=0;
			n2=0;

			for(j=i-1;j>=0;j--) {			 
				//if(se->seqngap[j]!='-'&&owner->seqngap[j]!='-') {
				if(addon[j]==0) {
					n1++;
				}
				else {
					break;
				}
			}

			for(j=i+1;j<nlen;j++) {
				//if(se->seqngap[j]!='-'&&owner->seqngap[j]!='-') {
				if(addon[j]==0) {
					n2++;
				}
				else {
					break;
				}
			}

			if(n1==0&&n2==0) continue;
			r=resn[match[i]];
			if(n2>n1) {
				r0=r->next;
				rot.link(r,r0,r->id0,0);
				addon[i]=2;tot++;
			}
			else {
				r0=resn[match[i-1]];
				rot.link(r0,r,r->id0,1);
				addon[i]=2;tot++;
			} 
		}

		cerr<<"the total number of residues assembled:"<<tot<<endl;
		if(tot==0) {
			for(i=m1+1;i<m2;i++) {
				if(addon[i]==2) {
					addon[i]=0;
					tot++;
				}
			} 
		}		
		if(tot==0) break;
	}
	delete [] addon;addon=0;

	if(TRES.logg)mpdb->chn->write("mutate4.pdb");
}
void Mutate::linkstructure(int mm) {
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;
	
	//...
	cerr<<getsecondary()<<endl;
	cerr<<sqnto<<endl;
	cerr<<se->seqngap<<endl;
	//...

	if(se->seqngap[mm]=='-'&&owner->seqngap[mm]!='-') return; //delete case, no link required

	int nlen=strlen(se->seqngap);
		 
	int i=0;	

	int m1,m2;

	for(m1=mm-1;m1>=0;m1--) {
		if(owner->seqngap[m1]!='-'&&se->seqngap[m1]!='-') break;									 
	}
		
	for(m2=mm+1;m2<nlen;m2++) {
		if(owner->seqngap[m2]!='-'&&se->seqngap[m2]!='-') break;	 		 	 	
	}	 
 
	Rotate rot;
	int j;//,n1,n2;
	Res *r,*r0,*r1;
 	Chn *cc=mpdb->chn;
	while(1) {
	   	int tot=0;
		for(i=m1+1;i<m2;i++) {
			if(sqnto[i]=='-') continue;
			if(i!=mm) continue;			 
		 	j=match[i];
			r=resn[j];
			if(r==0) continue;			 
			r0=r->last;
			r1=r->next;
			if(cc->islinked(r0,r)==1||cc->islinked(r,r1)==1) continue;
			if(r0==0) {
				if(mpdb->chn->islinked(r1,r1->next)==1) {
					rot.link(r,r1,r->id0,0);
					tot++;
				}
			}			
			else if(r1==0) {
				if(mpdb->chn->islinked(r0->last,r0)==1) {
					rot.link(r0,r,r->id0,1);
					tot++;
				}
			}
			else {
				int n1=compare[r0->id0];
				int n2=compare[r1->id0];
				if(se->seqngap[n1]!='-'&&owner->seqngap[n1]!='-') {
					rot.link(r0,r,r->id0,1); 
					tot++;
				}
				else if(se->seqngap[n2]!='-'&&owner->seqngap[n2]!='-') {
					rot.link(r,r1,r->id0,0); 
					tot++;
				}  
				else if(mpdb->chn->islinked(r0->last,r0)==1) {
					rot.link(r0,r,r->id0,1);
					tot++;
				}
				else if(mpdb->chn->islinked(r1,r1->next)==1) {
					rot.link(r,r1,r->id0,0);
					tot++;
				} 						
			}			 
		}
		if(tot==0) break;		 
	}
	if(TRES.logg)mpdb->chn->write("mutate4.pdb");
}

void Mutate::setseqnres(){

	//remove residue does not exist in structure
	int i;

	//remove residue does not exist in sequence

	Res *r,*r0;

	int m;

	int n=strlen(owner->seqn);

	r=mpdb->chn->res;
	r0=0;
	while(r) 
	{
		m=0; 
		for(i=0;i<n;i++) {
			if(owner->resn[i]==0) continue;
			 
			if(r->id==owner->resn[i]->id) {
				m=1;
				break;
			}
		}
		if(m==0) {
			if(r0) {
				r0->next=r->next;
				r->next=0;
				delete r;r=0;
				r=r0->next;
			}
			else   {
				mpdb->chn->res=r->next;
				r->next=0;
				delete r;r=0;	
				r=mpdb->chn->res;
			}	
		}
		else {
			r0=r;
			r=r->next;
		}
	}
}

void Mutate::actsinglemutation(Res *r) {
   

 	mpdb->setflgr(-99999);
	
	r->flag=10000; 
	Scap scprd; 	 		
	scprd.pdb=mpdb;
	//08/16/2002
	strcpy(scprd.force,"124");
	TRES.setdonaldcharge();
	scprd.dielectric=15;	 
	//end
	float d=scprd.mutate();
	//08/16/2002
	TRES.switchcharge(1);
	//
	r->flag=-99999;	
	r->setfreenextonmore();	
	if(r->more) delete r->more;r->more=0;	
				 					 
} 

float Mutate::minimize(float en) {
	return minimize(mpdb->chn->res,100000,en);
}

 

float Mutate::minimize(Res *start,int end,float en) {
	
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);
	int *status=new int[nlen];
	
	mpdb->chn->allnearbond(3);
	setmatch();
        mpdb->setlastresorder();
	sethbonddssp();	
	setlocaldssp();
		 
	refine=-1;
 
	updatepdbold();
	updatepdbcopy(); 
	setunchanged();
	
	int tot=0;
	Res *r; 
	int i; 
	 
	float entot=0,enpre=0;
	float emax=0;
	while(1) {						
		 		
		setcontact(start,end);
		mpdb->chn->setresenergy(start,end);
		mpdb->chn->setresenergy(start,end,5);
		for(i=0;i<nlen;i++) status[i]=0;
		Res **tmp=getlist(start,end);
		int n=0;while(tmp[n])n++;
		entot=0;
		for(i=0;i<n;i++) entot+=tmp[i]->energy;
		if(tot%3!=2)  segen->onlysidechain=1;			 
		else  segen->onlysidechain=0;
		emax=tmp[0]->energy;
		for(i=0;i<n;i++) {
			if(tmp[i]->id0<start->id0) continue;
			if(tmp[i]->id0>end) continue;
			if(tmp[i]==0||status[tmp[i]->id0]==1) continue;			
			if(tmp[i]->energy<en) continue;
			 
			Res *rr=tmp[i];
			 
			segen->onlysidechain=1;
			Res **stem=finddsspsegment(rr);	
			if(stem==0) continue;
			if(TRES.logg>3) printsegment(stem); 
			for(r=stem[0];r;r=r->next) {
				if(r->id0>stem[1]->id0) break;
				status[r->id0]=1; 
			}
			segen->readyminfix(stem[0],stem[1]->id0);
			if(xyzself) delete [] xyzself;xyzself=0;
                        xyzself=getpdboldtransfer(stem[0],stem[1]); 
			float d=segen->predtmin();

			if(rely[rr->id0]==0&&segen->onlysidechain==0) {
				delete [] stem;	stem=0;	
				continue;		
			}
			//minimize backbone
			segen->onlysidechain=0;					
			stem=findchanged(rr,stem);		
			if(stem==0) continue;	

			if(stem[0]->id0<start->id0)stem[0]=start;
                        if(stem[1]->id0>end) stem[1]=mpdb->chn->findsmallres(end);
                        if(stem[1]==0) continue;

			 	
			if(TRES.logg>3) printsegment(stem); 	
			segen->readyminfix(stem[0],stem[1]->id0);
 			if(xyzself) delete [] xyzself;xyzself=0;
			xyzself=getpdboldtransfer(stem[0],stem[1]); 
			d=segen->predtmin();			 
			delete [] stem;	stem=0;				
		}
		 
		if(TRES.logg>3)mpdb->chn->write("out.pdb");
		if(tmp) delete [] tmp;tmp=0;	
		 
		if(TRES.logg>3) cerr<<"the difference between energy..."<<enpre<<"  "<<entot<<" "<<entot-enpre<<endl;
		if(fabs(entot-enpre)<0.1&&fabs(entot-enpre)<fabs(0.01*entot)&&tot%2==1) break;
		
		tot++;
		if(tot>10) break;
		if(fapr==5&&tot>=2) break;
		enpre=entot;
	}

	return emax;
}

float Mutate::minimizeold(Res *start,int end,float en) {
	
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nlen=strlen(sqnto);
	int *status=new int[nlen];
	
	mpdb->chn->allnearbond(3);
	setmatch();
        mpdb->setlastresorder();
	sethbonddssp();	
	setlocaldssp();
		 
	refine=-1;
 
	updatepdbold();
	updatepdbcopy(); 
	setunchanged();
	
	int tot=0;
	Res *r; 
	int i; 
	 
	float entot=0,enpre=0;
	float emax=0;
	while(1) {						
		//t1=0; t2=0; erg=0;			
		setcontact(start,end);
		mpdb->chn->setresenergy(start,end);
		mpdb->chn->setresenergy(start,end,5);
		for(i=0;i<nlen;i++) status[i]=0;
		Res **tmp=getlist(start,end);
		int n=0;while(tmp[n])n++;
		entot=0;
		for(i=0;i<n;i++) entot+=tmp[i]->energy;
		if(tot%3!=2)  segen->onlysidechain=1;			 
		else  segen->onlysidechain=0;
		emax=tmp[0]->energy;
		for(i=0;i<n;i++) {
			if(tmp[i]->id0<start->id0) continue;
			if(tmp[i]->id0>end) continue;
			if(tmp[i]==0||status[tmp[i]->id0]==1) continue;			
			if(tmp[i]->energy<en) continue;
			 
			Res *rr=tmp[i];
			 
			
			Res **stem=finddsspsegment(rr);	
			if(TRES.logg>3) printsegment(stem); 
			if(segen->onlysidechain==0) stem=findchanged(rr,stem);		
			if(stem==0) continue;	
			 	
			if(TRES.logg>3) printsegment(stem); 	
			for(r=stem[0];r;r=r->next) {
				if(r->id0>stem[1]->id0) break;
				status[r->id0]=1; 
			}
 			if(segen->onlysidechain==0) {
				if(rr->sec=='-'){
					if(stem[0]->id0<start->id0) stem[0]=start;
					if(stem[1]->id0>end) stem[1]=mpdb->chn->findsmallres(end);
					if(stem[1]==0) continue;
				}
				else {
					if(stem[0]->id0+6<start->id0) stem[0]=0;
					if(stem[1]->id0-6>end) stem[1]=0;
					if(stem[1]==0) continue;
				}
			}
			if(stem==0||stem[0]==0||stem[1]==0) continue;
			segen->readyminfix(stem[0],stem[1]->id0);
			if(xyzself) delete [] xyzself;xyzself=0;
                        xyzself=getpdboldtransfer(stem[0],stem[1]); 

			//disregard
			if(rely[rr->id0]==0&&segen->onlysidechain==0) continue;
			
 			//minimize			 
			float d=segen->predtmin();
			 
			delete [] stem;	stem=0;				
		}
		 
		if(TRES.logg>3)mpdb->chn->write("out.pdb");
		if(tmp) delete [] tmp;tmp=0;	
		if(TRES.logg>3) cerr<<"the difference between energy..."<<enpre<<"  "<<entot<<" "<<entot-enpre<<endl;
		if(fabs(entot-enpre)<0.1&&fabs(entot-enpre)<fabs(0.01*entot)&&tot%2==1) break;
		
		tot++;
		if(tot>10) break;
		if(fapr==5&&tot>=2) break;
		enpre=entot;
	}

	return emax;
}
void Mutate::minimize0(Res *start,int end,float en) {
	
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	int nlen=strlen(sqnto);
	int *status=new int[nlen];
	
	mpdb->chn->allnearbond(3);
	setmatch();
        mpdb->setlastresorder();	
	setlocaldssp();
	
	setboundscore();
	for(int m=0;m<nlen;m++) {
		boundscore[m]=sqrt(boundscore[m]);
		int i=m;
		cerr<<i<<"new bound ..."<<boundscore[i]<<" "<<se->seqngap[i]<<" "<<sqnto[i]<<" "<<endl;
	}

	refine=0;
	Res *r,*t1,*t2;
	int i,i1,i2;
	float erg=0;
	char line[100];
	sprintf(line,"%s_ini.pdb",code);
	mpdb->write(line);
 	 
	int tot=0;
	updatepdbcopy(); 
	setunchanged();
	float entot=0,enpre=0;
	while(1) {
						
		t1=0; t2=0; erg=0;			
		setcontact(start,end);
		mpdb->chn->setresenergy(start,end);
		mpdb->chn->setresenergy(start,end,5);
		for(i=0;i<nlen;i++) status[i]=0;
		Res **tmp=getlist(start,end);
		int n=0;while(tmp[n])n++;
		entot=0;
		for(i=0;i<n;i++) entot+=tmp[i]->energy;
		if(tot%3==0)  segen->onlysidechain=1;			 
		else  segen->onlysidechain=0;

		for(i=0;i<n;i++) {
			if(tmp[i]==0||status[tmp[i]->id0]==1) continue;			
			if(tmp[i]->energy<en) continue;
			r=tmp[i];
			Res *rt=r;
			int ie=compare[rt->id0]; 
			
			Res **stem=findminsegment(r);	
			if(TRES.logg) printsegment(stem); 
			if(segen->onlysidechain==0) stem=findchanged(r,stem);		
			if(stem==0) continue;	
			 	
			if(TRES.logg) printsegment(stem); 	
			for(r=stem[0];r;r=r->next) {
				if(r->id0>stem[1]->id0) break;
				status[r->id0]=1; 
			}
 
			segen->readyminfix(stem[0],stem[1]->id0); 
			if(rely[rt->id0]==0&&segen->onlysidechain==0) continue;
			if(dssp[ie]!='-'&&segen->onlysidechain==0) continue;
			float d=segen->predtmin();
			 
			if(t1==0||t2==0||d>erg) {
				t1=stem[0];t2=stem[1];erg=d;
			}
			delete [] stem;	stem=0;				
		}
		int it=0;
		while(erg>0&&t1&&t2&&tot%2==1) {
			int nn=(t1->id0+t2->id0)/2;
			int *stem=findloosesegment(nn,t1,t2);
			i1=match[stem[0]];i2=match[stem[1]];
			t1=resn[i1];t2=resn[i2];
			if(TRES.logg) printsegment(stem); 
			segen->readyminfix(t1,t2->id0); 
			float ert=segen->predtmin();	
			if(it>1||ert<0) break;
			it++;			
		}
		if(TRES.logg)mpdb->chn->write("out.pdb");
		if(tmp) delete [] tmp;tmp=0;	
		//if(erg<=0&&tot%2==1) break;
		cerr<<"the difference between energy..."<<enpre<<"  "<<entot<<" "<<entot-enpre<<endl;
		if(fabs(entot-enpre)<fabs(0.01*entot)&&tot%2==1) break;
		
		tot++;
		if(tot>10) break;
		if(fapr==5&&tot>=2) break;
		enpre=entot;
	}
}
 
Res **Mutate::findchanged(Res *rr,Res **stem){

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;
	
	int nlen=strlen(sqnto);

	//Res *r;
	
	int n=compare[rr->id0];

	if(dssp[n]!='-') return stem;
	if(stem==0) return stem;

	Res *r1,*r2;
	Res *r10,*r20;
	r10=0;
	for(r1=rr->last;r1;r1=r1->last) {
		if(r1->id0<stem[0]->id0) break;
		if(rely[r1->id0]==0) break;
		r10=r1;
	}
	r20=0;
	for(r2=rr->next;r2;r2=r2->next) {
		if(r2->id0>stem[1]->id0) break;
		if(rely[r2->id0]==0) break;
		r20=r2;
	}
	r1=r10;
	r2=r20;
	
	if(r1==0) r1=stem[0];
	if(r2==0) r2=stem[1];
	if(r1->last&&r1->id0>stem[0]->id0) r1=r1->last;
	if(r2->next&&r2->id0<stem[1]->id0) r2=r2->next;
	/*
	if(r2->id0-r1->id0+1<=3) {
		if(r1->last) r1=r1->last;
		if(r2->next) r2=r2->next;
	}
	*/
	stem[0]=r1;
	stem[1]=r2;
	return stem;
}

Res **Mutate::finddsspsegment(Res *rr) {

	Res *r1,*r2;

	for(r1=rr->last;r1;r1=r1->last) {
		if(r1->sec!=rr->sec) break;
	}
	if(r1==0) r1=rr->chn->res;

	for(r2=rr->next;r2;r2=r2->next) {
		if(r2->sec!=rr->sec) break;
	}
	if(r2==0) r2=rr->chn->lastres();

	Res *v1,*v2;

	if(rr->sec!='-') {
		v1=r1;
      		//find boundary of v1
          	if(v1->sec=='-') {
                 	while(v1->last&&v1->last->sec=='-'&&r1->id0-v1->id0<2) v1=v1->last;
          	}
 
            	v2=r2;
        	//find boundary of v2
           	if(v2->sec=='-') {
               		while(v2->next&&v2->next->sec=='-'&&v2->id-r2->id0<2) v2=v2->next;
           	}
	}
	else {
		if(r2->id0-r1->id0+1>8) {
			v1=mpdb->chn->isres0(rr->id0-4);
			if(v1==0||v1->id0<r1->id0) v1=r1;
			v2=mpdb->chn->isres0(rr->id0+4);
			if(v2==0||v2->id0>r2->id0) v2=r2;		  
		}
		else {
			v1=r1;
			v2=r2;
		}
	}
	Res **stem=new Res*[2];
	stem[0]=v1;
	stem[1]=v2;
	return stem;
}

Res **Mutate::findminsegment(Res *rr) {
	
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;
	
	int nlen=strlen(sqnto);

	//Res *r;
	
	int n=compare[rr->id0];

	int n1,n2;

	for(n1=n-1;n1>=0;n1--) {
                if(dssp[n1]==dssp[n]) continue;
                if(owner->seqngap[n1]!='-'&&se->seqngap[n1]!='-') break;
        }
           
	
        for(n2=n+1;n2<nlen;n2++) {
                if(dssp[n2]==dssp[n]) continue;
                if(owner->seqngap[n2]!='-'&&se->seqngap[n2]!='-') break;
        }
	
	if(dssp[n]!='-') {

		int i1,i2;	

		for(i1=n1;i1>=0;i1--) {
			if(dssp[i1]==dssp[n1]) continue;
                	if(owner->seqngap[i1]!='-'&&se->seqngap[i1]!='-') break;
		}

		for(i2=n2;i2<nlen;i2++) {
                	if(dssp[i2]==dssp[n2]) continue;
                	if(owner->seqngap[i2]!='-'&&se->seqngap[i2]!='-') break;
       		} 
	
		if(i1<n1&&dssp[i1]=='-') n1--;
		if(i2>n2&&dssp[i2]=='-') n2++;		
	}
  	else {
		n1++;
		n2--;
	}

	
	Res *r1,*r2;

	int i1,i2;
	
	while(1)  {
		if(n1<0)  {
			i1=mpdb->chn->res->id0;
		}
		else i1=match[n1];
	 	if(i1!=-1) break;
		else n1--;
	}

	while(1) {
		if(n2>=nlen) {
			i2=mpdb->chn->lastres()->id0;
		}
		else i2=match[n2];
		if(i2!=-1) break;
		else n2++;
	}	

	r1=resn[i1];
	r2=resn[i2];

	n1=compare[r1->id0];
	n2=compare[r2->id0];

	int p1=0;
	int p2=0;
        int i;
	for(i=n1;i<=n2;i++) {	
		if(owner->seqngap[i]!='-') p1++;
		if(se->seqngap[i]!='-') p2++;	
 	}
 
	if(r2->id0-r1->id0+1==2) {
		int *stem=findloosesegment(r1,r2,2);
		int i1=match[stem[0]];
		int i2=match[stem[1]];
		r1=resn[i1];r2=resn[i2];
		delete [] stem;		stem=0;
	}

	if(dssp[n]=='-') {
		Res **tmp=new Res*[2];
		tmp[0]=r1;tmp[1]=r2;
		return tmp;
	}
 	else if(dssp[n]=='h') { 
		Res **tmp=new Res*[2];
		tmp[0]=r1;tmp[1]=r2;
		return tmp;  
	}
	else if(dssp[n]=='e') {
		Res **tmp=new Res*[2];
		tmp[0]=r1;tmp[1]=r2;
		return tmp; 
	}
	return 0;
}

Res *Mutate::findsegmaxclash(Res *r1,int n) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

        Res *r;
        Res *s0=0;
        float d=-999999;
	
	int nn1=compare[r1->id0];
	r=mpdb->chn->isres0(n); 
	if(r==0) r=mpdb->chn->lastres();
	int nn2=compare[r->id0];
		
	int i;
	for(i=nn1;i<=nn2;i++) {
		if(owner->seqngap[i]=='-'&&se->seqngap[i]=='-') continue;
		if(owner->seqngap[i]!='-'&&se->seqngap[i]!='-') continue;
		int n=match[i];
		if(n==-1) continue;
		r=resn[n];
		if(s0==0||r->energy>d) {
		 	s0=r;
                        d=r->energy;	
		}	
	}
       
        return s0;
}

Res **Mutate::getlist() {

	return getlist(mpdb->chn->res,100000);
}

Res *Mutate::getres(Res *s) {
	
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int n=compare[s->id0];
	int i=owner->match[n];
	if(i==-1) return 0;
	return owner->resn[i];
}

char Mutate::getmodelresname(Res *s) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return '-';

	int n=owner->findpost(s);
	if(n==-1) return '-';
	return se->seqngap[n];
}

Res *Mutate::getmodelres(Res *s) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int n=owner->findpost(s);
	if(n==-1) return 0;
	int i=match[n];
	if(i==-1) return 0;
	return resn[i];
}



Res **Mutate::getlist(Res *start,int end) {
 		
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	int nres=strlen(sqnto);
	
	Res **ren=new Res*[nres];
	float *temp=new float[nres];
	int   *order=new int[nres];

	int n=0;
	Res *r;
	float d=0;
	for(r=start;r;r=r->next) {
		if(r->id0>end) break;		 
		ren[n]=r;		
		temp[n]=-r->energy;
		if(r==start||r->energy>d) {
			d=r->energy;
		}		 
		n++;
	}
 	if(d<0) d=0;

	int i;
	for(i=0;i<n;i++) {
		r=ren[i];
		int j=compare[r->id0];
		if(owner->seqngap[j]=='-'||se->seqngap[j]=='-') {
			if(r->energy>0) temp[i]=-d*2;
		}
	}

	Qsort cc;

	cc.sort(temp,n,order);	

	
	Res **ret=new Res*[nres+100];	 
	for(i=0;i<n;i++) {		
		int j=order[i];
		r=ren[j];
		ret[i]=r;
	}
	ret[n]=0;
	ret[n+1]=0;
	delete [] temp;temp=0;
	delete [] order;order=0;
	delete [] ren;ren=0;
	return ret;
}

void Mutate::minsc() {

	setcontact(mpdb->chn->res,100000);
	mpdb->chn->setresenergy(mpdb->chn->res,100000,4);
		
	int nres=mpdb->chn->lastres()->id0+100;
	
	Res **ren=new Res*[nres];

	float *temp=new float[nres];
	int   *order=new int[nres];

	int n=0;Res *r;
	for(r=mpdb->chn->res;r;r=r->next) {
		if(strchr("PGA",r->name)) continue;
		ren[n]=r;
		temp[n]=-r->energy;
		n++;
	}

	Qsort cc;

	cc.sort(temp,n,order);	

	int i;

	for(i=0;i<n;i++) {		
		int j=order[i];
		r=ren[j];
		if(r->energy<0) continue;
		actsinglemutation(r);
	}
	delete [] temp;temp=0;
	delete [] order;order=0;
	delete [] ren;ren=0;
}


void Mutate::minloop(Res *rr,int nn) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
        segen->flex=1;
        segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
        segen->randcoil=0;
	//TRES.smoothclash=1;
        segen->smoothclash=1;
        segen->part=0;
	segen->randcoil=2;
 
	segen->cutoff=8.;
	 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	mpdb->setlastresorder();
	mpdb->chn->dihedral((FILE*)0);
	int nlen=strlen(sqnto);
        int *status=new int[nlen];
		
	setnewunchanged();
	 
	Res *r;//,*t;
	int i;
	for(i=0;i<nlen;i++) status[i]=1;

	refine=100; //loop prediction

	while(1) {
		int tote=0; int tota=0;
		for(r=rr;r;r=r->next) {
			if(r->id0>nn) break;
			tota=0;
			//if(status[r->id0]==0) continue;
			if(rely[r->id0]==2) continue;
			//if(r->sec!='-') continue;
			//if(rely[r->id0]==1) continue;
			Res *r1,*r2;
			for(r1=r;r1;r1=r1->last){
				//if(r1->sec!='-') break;
				//if(status[r1->id0]==0) break;
				if(rely[r1->id0]>1) break;			 
			}
			if(r1==0) r1=mpdb->chn->res;
			//else	  r1=r1->next;
			for(r2=r;r2;r2=r2->next){
				//if(r2->sec!='-') break;
				//if(status[r2->id0]==0) break;	
				if(rely[r2->id0]>1) break;						 
			}
			if(r2==0) r2=mpdb->chn->lastres();
			//else      r2=r2->last;

			if(r1==0) r1=mpdb->chn->res;
			if(r2==0) r2=mpdb->chn->lastres();		

			Res *t;
			Res *t1,*t2;
				
			t=findminscore(status,r1,r2->id0); 		
			if(t==0) goto re200;

			t1=r1;
			t2=r2;
			if(aloop==0) {
				for(t1=t->last;t1;t1=t1->last) {
					//if(t1->id0<=r1->id0) break;
					if(rely[t1->id0]>0) break;
				}

				for(t2=t->next;t2;t2=t2->next) {
					//if(t2->id0>=r2->id0) break;
					if(rely[t2->id0]>0) break;
				}
		 	}			 
			if(t1==0) t1=mpdb->chn->res;
			if(t2==0) t2=mpdb->chn->lastres();

			int n1;n1=calcinsert(t1,t2->id0);
			int n2;n2=calcdelete(t1,t2->id0);
			
			if(n1==0&&n2==0) goto re200;			

			 
			if(t2->id0-t1->id0+1<5&&n1==0) {				
				if(t2->next) t2=t2->next;
				if(t1->last) t1=t1->last;
				/*
				int n1=compare[t1->id0];
				int n2=compare[t2->id0];
				
				if(parent->sitescore[n1]<parent->sitescore[n2]) {
					if(t1->last) t1=t1->last;
				}
				else {
					if(t2->next) t2=t2->next;
				}	
				*/
			}

			if(t2->id0-t1->id0+1<5&&n2==0) {
				if(t2->next) t2=t2->next;
				if(t1->last) t1=t1->last;				
			}
			 

			if(t2->id0-t1->id0+1<5) goto re200;	 
			
			int no;no=0;

			for(t=t1;t;t=t->next) {
				if(t->id0>t2->id0) break;
				if(status[t->id0]==1) {
					no=1;break;
				}
			}

			if(no==0) goto re200;

			if(t2->id0-t1->id0+1>=10) {
				segen->arbt=200;
			}
			else if(t2->id0-t1->id0+1>=8) {
				segen->arbt=100;
			}
			else if(t2->id0-t1->id0+1>=6) {
				segen->arbt=50;
			}
			else if(t2->id0-t1->id0+1>=4) {
				segen->arbt=30;
			}
			else{
				segen->arbt=20;
			}
 			
			//if(n1==0) segen->arbt=segen->arbt/2;
			if(test) segen->arbt=10; 
			if((t2->id0-t1->id0+1-n1)>6&&aloop==0) goto re200;
  
			int *stem;stem=new int[2];
			stem[0]=compare[t1->id0];
			stem[1]=compare[t2->id0];
			float *xyzorg;xyzorg=mpdb->chn->gettransfer(t1,t2->id0);
			if(xyzself) delete [] xyzself;xyzself=0;
			xyzself=mpdb->chn->gettransfer(t1,t2->id0);
			tote++;tota++;
			mpdb->chn->dihedral();
			if(TRES.logg) printsegment(stem);
			float *xyzout;xyzout=segen->myfixsegment(stem);		
			if(xyzout) {				 
				mpdb->chn->transfer(xyzout,t1,t2->id0,1);
				segen->disc->setupall(mpdb,10,1);
				strcpy(segen->disc->force,"uDd");
				float d1;d1=segen->disc->clash(mpdb->chn->res,100000);
				//
				mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
				segen->disc->setupall(mpdb,10,1);
				strcpy(segen->disc->force,"uDd");
				float d2;d2=segen->disc->clash(mpdb->chn->res,100000);	
				if(d1<d2) mpdb->chn->transfer(xyzout,t1,t2->id0,1);
			}
			else {
				mpdb->chn->transfer(xyzorg,t1,t2->id0,1);				 
			}

			if(xyzout) delete [] xyzout;xyzout=0;
			if(xyzorg) delete [] xyzorg;xyzorg=0;

			re200:
			 
			for(t=t1;t&&tota;t=t->next) {
				if(t->id0>t2->id0) break;
				status[t->id0]=0;
			}	
				 							
		} 

		if(tote==0) break;
	}

	

	if(status) delete [] status;status=0;
}


void Mutate::minwindowloop(Res *rr,int nn) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
        segen->flex=1;
        segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
        segen->randcoil=0;
	//TRES.smoothclash=1;
        segen->smoothclash=1;
        segen->part=0;
	segen->randcoil=2;
	segen->cutoff=8.;
	 

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	refine=100; //loop prediction

	mpdb->setlastresorder();
	mpdb->chn->dihedral((FILE*)0);
	int nlen=strlen(sqnto);
        int *status=new int[nlen];
		
	setnewunchanged();
	 
	Res *r;//,*t;
	int i;
	for(i=0;i<nlen;i++) status[i]=1;
 
	
	int wsize=3;
	while(1) {
		int tote=0; int tota=0;
		wsize+=2;
		if(wsize==7) wsize=6;
		if(wsize>6) break;
		for(i=0;i<nlen;i++) status[i]=1;
		for(r=rr;r;r=r->next) {
			if(r->id0>nn) break;
			tota=0;
			//if(status[r->id0]==0) continue;
			if(rely[r->id0]==2) continue;
			//if(r->sec!='-') continue;
			//if(rely[r->id0]==1) continue;
			Res *r1,*r2;
			for(r1=r;r1;r1=r1->last){
				//if(r1->sec!='-') break;
				//if(status[r1->id0]==0) break;
				if(rely[r1->id0]>1) break;			 
			}
			if(r1==0) r1=mpdb->chn->res;
			//else	  r1=r1->next;
			for(r2=r;r2;r2=r2->next){
				//if(r2->sec!='-') break;
				//if(status[r2->id0]==0) break;	
				if(rely[r2->id0]>1) break;						 
			}
			if(r2==0) r2=mpdb->chn->lastres();
			//else      r2=r2->last;
			Res *t;
			if(r1==0) r1=mpdb->chn->res;
			if(r2==0) r2=mpdb->chn->lastres();		
			int n1=calcinsert(r1,r2->id0);
			int n2=calcdelete(r1,r2->id0);			
			if(n1==0&&n2==0) {
				for(t=r1;t;t=t->next) {
					if(t->id0>r2->id0) break;
					status[t->id0]=0;
				}	
				continue;
			}
		
			
			Res *t1,*t2;
				
			t=findminscore(status,r1,r2->id0); 		
			if(t==0) goto re200;
			
			t1=r1;
			t2=r2;
			int nt1;nt1=r1->id0;
			int nt2;nt2=r2->id0;
			int nt;nt=0;
			for(nt=nt1;nt<=nt2;nt++) { 	
				 		
				t1=mpdb->chn->isres0(nt);				
				int ntx;ntx=nt+wsize;
				if(ntx>=nt2) {
					ntx=nt2;
					nt=nt2;
				}				
				t2=mpdb->chn->isres0(ntx);

				if(t1==0) t1=r1;
				if(t2==0) t2=r2;

				int n1;n1=calcinsert(t1,t2->id0);
				int n2;n2=calcdelete(t1,t2->id0);
			
				//if(n1==0&&n2==0) goto re200;			
			 
				//if(t2->id0-t1->id0+1<5&&n1==0) {				
				//	if(t2->next) t2=t2->next;
				//	if(t1->last) t1=t1->last;				
				//}

				//if(t2->id0-t1->id0+1<5&&n2==0) {
				//	if(t2->next) t2=t2->next;
				//	if(t1->last) t1=t1->last;				
				//}			 

				//if(t2->id0-t1->id0+1<5) goto re200;	 
			
				int no;no=0;

				for(t=t1;t;t=t->next) {
					if(t->id0>t2->id0) break;
					if(status[t->id0]==1) {
						no=1;break;
					}
				}

				if(no==0) goto re200;

				if(t2->id0-t1->id0+1>=10) {
					segen->arbt=200;
				}
				else if(t2->id0-t1->id0+1>=8) {
					segen->arbt=100;
				}
				else if(t2->id0-t1->id0+1>=6) {
					segen->arbt=50;
				}
				else if(t2->id0-t1->id0+1>=4) {
					segen->arbt=30;
				}
				else{
					segen->arbt=20;
				}
 			
				if(n1==0) segen->arbt=segen->arbt/2;
				if(test) segen->arbt=5; 
				//if((t2->id0-t1->id0+1-n1)>6&&aloop==0) goto re200;
  
				int *stem;stem=new int[2];
				stem[0]=compare[t1->id0];
				stem[1]=compare[t2->id0];
				float *xyzorg;xyzorg=mpdb->chn->gettransfer(t1,t2->id0);
				if(xyzself) delete [] xyzself;xyzself=0;
				xyzself=mpdb->chn->gettransfer(t1,t2->id0);
				tote++;tota++;
				mpdb->chn->dihedral();
				if(TRES.logg) printsegment(stem);
				float *xyzout;xyzout=segen->myfixsegment(stem);		
				if(xyzout) {				 
					mpdb->chn->transfer(xyzout,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d1;d1=segen->disc->clash(mpdb->chn->res,100000);
					//
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d2;d2=segen->disc->clash(mpdb->chn->res,100000);	
					if(d1<d2) mpdb->chn->transfer(xyzout,t1,t2->id0,1);
				}
				else {
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);				 
				}

				if(xyzout) delete [] xyzout;xyzout=0;
				if(xyzorg) delete [] xyzorg;xyzorg=0;

				re200:
			 
				for(t=t1;t&&tota;t=t->next) {
					if(t->id0>t2->id0) break;
					status[t->id0]=0;
				}	
			}	 							
		} 

		if(tote==0) break;
	}

 

	if(status) delete [] status;status=0;
}

void Mutate::refineloop(Res *rr,int nn) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
        segen->flex=1;
        segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
        segen->randcoil=0;
        segen->smoothclash=0;
        segen->part=0;
	segen->randcoil=2;
	segen->cutoff=6.;
	segen->onlyenergy=0;
	
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	refine=100; //loop prediction

	mpdb->header();
	mpdb->setlastresorder();
	mpdb->chn->dihedral((FILE*)0);
	int nlen=strlen(sqnto);
        int *status=new int[nlen];
	
	//calculate unchanged regions.
	
	setnewunchanged();
	 
	Res *r; 

	//status indicates if minimized.
	int i;
	for(i=0;i<nlen;i++) status[i]=1; 
	
	//window size
	int wsize=3;
 	
	int ih;	  
	for(ih=0;ih<1;ih++) {
		int tote=0; int tota=0;	 
		//wsize+=2;		
		//if(wsize>9) break;
		for(i=0;i<nlen;i++) status[i]=1;
		for(r=rr;r;r=r->next) {
			if(r->id0>nn) break;			
			if(rely[r->id0]==2&&r->sec!='-') continue;//only in secondary and intact, go away 				
			if(rely[r->id0]==2&&sharp==2) continue;	//only intact go away			 				 				
			if(status[r->id0]==0) continue;
			tota=0;		

			//find boundary in two cases: 
			//a) only inserted or deleted regions
			//b) in all loop regions
							 
			Res *r1,*r2;
			for(r1=r;r1;r1=r1->last){
				if(rely[r1->id0]==2&&r1->sec!='-') break;				
				if(rely[r1->id0]>1&&sharp==2) break;						 
			}
			if(r1==0) r1=mpdb->chn->res;
			 
			for(r2=r;r2;r2=r2->next){
				if(rely[r2->id0]==2&&r2->sec!='-') break;				
				if(rely[r2->id0]>1&&sharp==2) break;				 							 
			}
			if(r2==0) r2=mpdb->chn->lastres();
			 
			Res *t;
			if(r1==0) r1=mpdb->chn->res;
			if(r2==0) r2=mpdb->chn->lastres();

			//check boundaries	
			if(r1->id0<rr->id0) r1=rr;
			if(r2->id0>nn) r2=mpdb->chn->findsmallres(nn);
			if(r2==0) continue;

			//check 
			int n1=calcinsert(r1,r2->id0);
			int n2=calcdelete(r1,r2->id0);

			//no insertion and deletion case			
			if(n1==0&&n2==0&&sharp==2) {
				for(t=r1;t;t=t->next) {
					if(t->id0>r2->id0) break;
					status[t->id0]=0;
				}	
				continue;
			}

			if(n1==0&&n2==0&&r2->id0-r1->id0+1<4) {
				status[r1->id0]=0;
				for(t=r1;t;t=t->next) {				
					if(t->id0>r2->id0) break;					
					else status[t->id0]=0;
				}
				continue;
			}	

			if(r1->last&&r2->next&&r2->id0-r1->id0+1<3) {
				status[r1->id0]=0;
				for(t=r1;t;t=t->next) {				
					if(t->id0>r2->id0) break;					
					else status[t->id0]=0;
				}
				continue;				
			}

			Res *t1,*t2;
				
			//find the most unconserved region
			t=findminscore(status,r1,r2->id0); 	
	
			if(t==0) {
				status[r1->id0]=0;
				for(t=r1;t;t=t->next) {				
					if(t->id0>r2->id0) break;					
					else status[t->id0]=0;
				}
				continue;
			}
 
			t1=r1;
			t2=r2;
			int nt1=r1->id0;
			int nt2=r2->id0;
			int nt=0;int ntx=0;
			
			//minimize the region of length wsize
			//for(nt=nt1;nt<=nt2;nt+=(wsize/2+0.7)) {
			//
			wsize=5;
			cerr<<"refine regions between: "<<r1->name<<r1->id0<<"---"<<r2->name<<r2->id0<<endl;
			while(1) {  
				int ns;ns=0;
				int y1;y1=0;
				for(y1=r1->id0;y1<r2->id0;y1++)  if(status[y1]==1) ns++;
				if(ns==0&&wsize<r2->id0-r1->id0&&wsize<7) {
					wsize+=4;
					for(y1=r1->id0;y1<r2->id0;y1++) status[y1]=1;
				}
				
				t=findminscore(status,r1,r2->id0);
				if(t==0) break;
				if(wsize+1<r2->id0-r1->id0+1) {
					nt=t->id0-(wsize+1)/2;
					if(nt<nt1) nt=nt1;
					ntx=t->id0+(wsize+1)/2;
					if(ntx>nt2) ntx=nt2;	
				}
				else {
					nt=r1->id0;
					ntx=r2->id0;
				}
				Res *tsave;tsave=t;										
			//  			
				t1=mpdb->chn->isres0(nt);	
				/*			
				int ntx=nt+wsize;
				if(ntx>=nt2) {
					ntx=nt2;
					nt=nt2;
				}
				*/				
				t2=mpdb->chn->isres0(ntx);

				if(t1==0) t1=r1;
				if(t2==0) t2=r2;

				int n1;n1=calcinsert(t1,t2->id0);
				int n2;n2=calcdelete(t1,t2->id0);
 
				int no;no=0;

				for(t=t1;t;t=t->next) {
					if(t->id0>t2->id0) break;
					if(status[t->id0]==1) {
						no=1;break;
					}
				}
				
				//all have been done!
				if(no==0) goto re200;
				if(mpdb->chn->isalllinked(t1->id0-1,t2->id0+1)==0) {
					cerr<<"there is chain breaker in the segment: "<<t1->name<<t1->id0<<"---"<<t2->name<<t2->id0<<endl;
					cerr<<"ignore this segment..."<<endl;
					goto re200;
				}
				float coeff;coeff=1;

				if(fapr==5) 	 coeff=0.25;
				else if(fapr==4) coeff=0.5;
				else if(fapr==3) coeff=1;
				else if(fapr==2) coeff=2;
				else if(fapr==1) coeff=4;
				else if(fapr==0) coeff=8;
				
				//
				coeff=2*coeff;
				//
				if(t2->id0-t1->id0+1>=10) {					
					segen->arbt=(int)(200*coeff);
				}
				else if(t2->id0-t1->id0+1>=8) {
					segen->arbt=(int)(100*coeff);
				}
				else if(t2->id0-t1->id0+1>=6) {
					segen->arbt=(int)(50*coeff);
				}
				else if(t2->id0-t1->id0+1>=4) {
					segen->arbt=(int)(30*coeff);
				}
				else{
					segen->arbt=(int)(20*coeff);
				} 			
				if(n1==0) segen->arbt=segen->arbt/2;
				
				if(test) segen->arbt=5; 
				
				//if((t2->id0-t1->id0+1-n1)>6&&aloop==0) goto re200;
  				
				int *stem;stem=new int[2];
				stem[0]=compare[t1->id0];
				stem[1]=compare[t2->id0];
				float *xyzorg;xyzorg=mpdb->chn->gettransfer(t1,t2->id0);
				if(xyzself) delete [] xyzself;xyzself=0;
				xyzself=getpdboldtransfer(t1,t2);
				 
				tote++;tota++;
				segen->updatexyzsave(stem);
				updatepdbcopy();				 
				delsidechain(stem);
				mpdb->chn->dihedral();
				if(TRES.logg>3) printsegment(stem);	

				//
				//int sa1=calcinsert(t1,t2->id0);
				//int sa2=calcdelete(t1,t2->id0);	
				//float rmsd0=rmsd;
				//if(sa1||sa2) rmsd+=1;						
				//rmsd=rmsd0; 
				//float *xyzout=mpdb->chn->gettransfer(t1,t2->id0); //test...
				float *xyzout;xyzout=segen->myfixsegment(stem);
				if(stem) delete [] stem;stem=0;		
				if(xyzout==0){					 
					segen->hooksidechain(); 
					mpdb->configure();					
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1); 
					segen->setallnear();					
				}						
				if(xyzout) {	
					Res *ti=mpdb->chn->isres0(segen->start-3);
					if(ti==0) ti=mpdb->chn->res;	
					int endo=segen->end+3;	 
					mpdb->chn->transfer(xyzout,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d1=segen->disc->clash(ti,endo);
					//
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d2=segen->disc->clash(ti,endo);	
					if(d1<d2) {
						mpdb->chn->transfer(xyzout,t1,t2->id0,1);
						mpdb->chn->header(t1->id0-3,t2->id0+3);
					}
				}
				else {
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);				 
				}

				if(xyzout) delete [] xyzout;xyzout=0;
				if(xyzorg) delete [] xyzorg;xyzorg=0;

				re200:
				if(tsave) status[tsave->id0]=0;
			 	 		
				for(t=t1;t&&tota;t=t->next) {				
					if(t->id0>t2->id0) break;					
					else status[t->id0]=0;
				}	
				 
			}	 							
		} 

		if(tote==0) break;
	}

 

	if(status) delete [] status;status=0;
}

//minimize loop regions
void Mutate::refineminloop(Res *rr,int nn) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
        segen->flex=1;
        segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
        segen->randcoil=0;
        segen->smoothclash=0;
        segen->part=0;
	segen->randcoil=2;
	segen->cutoff=6.;
	segen->onlyenergy=0; 
	
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	refine=100; //loop prediction

	mpdb->setlastresorder();
	mpdb->chn->dihedral((FILE*)0);
	int nlen=strlen(sqnto);
        int *status=new int[nlen];
	
	//calculate unchanged regions.
	
	setnewunchanged();
	 
	Res *r; 

	//status indicates if minimized.
	int i;
	for(i=0;i<nlen;i++) status[i]=1; 
	
	//window size
	int wsize=3;
 	
	int ih;	  
	for(ih=0;ih<1;ih++) {
		int tote=0; int tota=0;	 
		//wsize+=2;		
		//if(wsize>9) break;
		for(i=0;i<nlen;i++) status[i]=1;
		for(r=rr;r;r=r->next) {
			if(r->id0>nn) break;			
			if(rely[r->id0]==2&&r->sec!='-') continue;//only in secondary and intact, go away 				
			if(rely[r->id0]==2&&sharp==2) continue;	//only intact go away			 				 				
			if(status[r->id0]==0) continue;
			tota=0;		

			//find boundary in two cases: 
			//a) only inserted or deleted regions
			//b) in all loop regions
							 
			Res *r1,*r2;
			for(r1=r;r1;r1=r1->last){
				if(rely[r1->id0]==2&&r1->sec!='-') break;				
				if(rely[r1->id0]>1&&sharp==2) break;						 
			}
			if(r1==0) r1=mpdb->chn->res;
			 
			for(r2=r;r2;r2=r2->next){
				if(rely[r2->id0]==2&&r2->sec!='-') break;				
				if(rely[r2->id0]>1&&sharp==2) break;				 							 
			}
			if(r2==0) r2=mpdb->chn->lastres();
			 
			Res *t;
			if(r1==0) r1=mpdb->chn->res;
			if(r2==0) r2=mpdb->chn->lastres();
	
			//check boundaries	
			if(r1->id0<rr->id0) r1=rr;
			if(r2->id0>nn) r2=mpdb->chn->findsmallres(nn);
			if(r2==0) continue;

			int n1=calcinsert(r1,r2->id0);
			int n2=calcdelete(r1,r2->id0);

			//no insertion and deletion case			
			if(n1==0&&n2==0&&sharp==2) {
				for(t=r1;t;t=t->next) {
					if(t->id0>r2->id0) break;
					status[t->id0]=0;
				}	
				continue;
			}

			if(n1==0&&n2==0&&r2->id0-r1->id0+1<4) {
				status[r1->id0]=0;
				for(t=r1;t;t=t->next) {				
					if(t->id0>r2->id0) break;					
					else status[t->id0]=0;
				}
				continue;
			}	

			if(r1->last&&r2->next&&r2->id0-r1->id0+1<3) {
				status[r1->id0]=0;
				for(t=r1;t;t=t->next) {				
					if(t->id0>r2->id0) break;					
					else status[t->id0]=0;
				}
				continue;				
			}

			Res *t1,*t2;
				
			//find the most unconserved region
			t=findminscore(status,r1,r2->id0); 	
	
			if(t==0) {
				status[r1->id0]=0;
				for(t=r1;t;t=t->next) {				
					if(t->id0>r2->id0) break;					
					else status[t->id0]=0;
				}
				continue;
			}
 
			t1=r1;
			t2=r2;
			int nt1=r1->id0;
			int nt2=r2->id0;
			int nt=0;int ntx=0;

			//minimize the region of length wsize
			//for(nt=nt1;nt<=nt2;nt+=(wsize/2+0.7)) {
			//
			wsize=5;
			cerr<<"minimize regions between: "<<r1->name<<r1->id0<<"---"<<r2->name<<r2->id0<<endl;
			while(1) {  
				int ns;ns=0;
				int y1;y1=0;
				for(y1=r1->id0;y1<r2->id0;y1++)  if(status[y1]==1) ns++;
				if(ns==0&&wsize<r2->id0-r1->id0&&wsize<7) {
					wsize+=4;
					for(y1=r1->id0;y1<r2->id0;y1++) status[y1]=1;
				}
				
				t=findminscore(status,r1,r2->id0);
				if(t==0) break;
				if(wsize+1<r2->id0-r1->id0+1) {
					nt=t->id0-(wsize+1)/2;
					if(nt<nt1) nt=nt1;
					ntx=t->id0+(wsize+1)/2;
					if(ntx>nt2) ntx=nt2;	
				}
				else {
					nt=r1->id0;
					ntx=r2->id0;
				}
				Res *tsave;tsave=t;										
			//  			
				t1=mpdb->chn->isres0(nt);	
				/*			
				int ntx=nt+wsize;
				if(ntx>=nt2) {
					ntx=nt2;
					nt=nt2;
				}
				*/				
				t2=mpdb->chn->isres0(ntx);

				if(t1==0) t1=r1;
				if(t2==0) t2=r2;

				int n1;n1=calcinsert(t1,t2->id0);
				int n2;n2=calcdelete(t1,t2->id0);
 
				int no;no=0;

				for(t=t1;t;t=t->next) {
					if(t->id0>t2->id0) break;
					if(status[t->id0]==1) {
						no=1;break;
					}
				}
				
				//all have been done!
				if(no==0) goto re200;
				if(mpdb->chn->isalllinked(t1->id0-1,t2->id0+1)==0) {
					cerr<<"chain breakers found in the segment: "<<t1->name<<t1->id0<<"---"<<t2->name<<t2->id0<<endl;
					cerr<<"ignore this segment..."<<endl;
					goto re200;
				}
				float coeff;coeff=1;

				if(fapr==5) 	 coeff=0.25;
				else if(fapr==4) coeff=0.5;
				else if(fapr==3) coeff=1;
				else if(fapr==2) coeff=2;
				else if(fapr==1) coeff=4;
				else if(fapr==0) coeff=8;
				
				//
				coeff=2*coeff;
				//
				if(t2->id0-t1->id0+1>=10) {					
					segen->arbt=(int)(200*coeff);
				}
				else if(t2->id0-t1->id0+1>=8) {
					segen->arbt=(int)(100*coeff);
				}
				else if(t2->id0-t1->id0+1>=6) {
					segen->arbt=(int)(50*coeff);
				}
				else if(t2->id0-t1->id0+1>=4) {
					segen->arbt=(int)(30*coeff);
				}
				else{
					segen->arbt=(int)(20*coeff);
				} 			
				if(n1==0) segen->arbt=segen->arbt/2;
				
				if(test) segen->arbt=5; 

				//if((t2->id0-t1->id0+1-n1)>6&&aloop==0) goto re200;
  
				int *stem;stem=new int[2];
				stem[0]=compare[t1->id0];
				stem[1]=compare[t2->id0];
				float *xyzorg;xyzorg=mpdb->chn->gettransfer(t1,t2->id0);
				if(xyzself) delete [] xyzself;xyzself=0;
				xyzself=getpdboldtransfer(t1,t2);
				 
				tote++;tota++;
				segen->updatexyzsave(stem);
				updatepdbcopy();				 
				//delsidechain(stem);
				mpdb->chn->dihedral();
				if(TRES.logg>3) printsegment(stem);	
 	
				segen->readyminfix(t1,t2->id0);	

				//minimize sidechain only
				segen->onlysidechain=1;
				float d;d=segen->predtmin();
				//mpdb->chn->transfer(xyzorg,t1,t2->id0,0);

 				//minimize backbone
				segen->onlysidechain=0;
				segen->predtmin();
				float *xyzout;xyzout=mpdb->chn->gettransfer(t1,t2->id0);
				 
				if(stem) delete [] stem;stem=0;		
				if(xyzout==0){					 
					segen->hooksidechain(); 
					mpdb->configure();					
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1); 
					segen->setallnear();					
				}						
				if(xyzout) {	
					Res *ti;ti=mpdb->chn->isres0(segen->start-3);
					if(ti==0) ti=mpdb->chn->res;	
					int endo;endo=segen->end+3;	 
					mpdb->chn->transfer(xyzout,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d1;d1=segen->disc->clash(ti,endo);
					//
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d2;d2=segen->disc->clash(ti,endo);	
					if(d1<d2) mpdb->chn->transfer(xyzout,t1,t2->id0,1);
				}
				else {
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);				 
				}

				if(xyzout) delete [] xyzout;xyzout=0;
				if(xyzorg) delete [] xyzorg;xyzorg=0;

				re200:
				if(tsave) status[tsave->id0]=0;
			 	//int nmid=(t1->id0+t2->id0)/2;				
				for(t=t1;t&&tota;t=t->next) {				
					if(t->id0>t2->id0) break;
					//else if(fabs(t->id0-nmid)<=(wsize/2-1)) status[t->id0]=0;
					//else if(t->id-t1->id0<2)  continue;
					//else if(t2->id0-t->id0<2) continue;
					else status[t->id0]=0;
				}	
				 
			}	 							
		} 

		if(tote==0) break;
	}

 

	if(status) delete [] status;status=0;
}

//minimize loop regions
void Mutate::refinesmoothloop(Res *rr,int nn) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
        segen->flex=1;
        segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
        segen->randcoil=0;
        segen->smoothclash=TRES.smoothclash;
        segen->part=0;
	segen->randcoil=2;
	segen->cutoff=6.;
	segen->onlyenergy=0; 
	
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	refine=-100; //loop prediction

	mpdb->setlastresorder();
	mpdb->chn->dihedral((FILE*)0);
	int nlen=strlen(sqnto);
        int *status=new int[nlen];
	
	//calculate unchanged regions.
	
	setnewunchanged();
	 
	Res *r; 

	//status indicates if minimized.
	int i;
	for(i=0;i<nlen;i++) status[i]=1; 
	
	//window size
	int wsize=3;
 	
	int ih;	  
	for(ih=0;ih<1;ih++) {
		int tote=0; int tota=0;	 
		//wsize+=2;		
		//if(wsize>9) break;
		for(i=0;i<nlen;i++) status[i]=1;
		for(r=rr;r;r=r->next) {
			if(r->id0>nn) break;			
			if(rely[r->id0]==2&&r->sec!='-') continue;//only in secondary and intact, go away 				
			if(rely[r->id0]==2&&sharp==2) continue;	//only intact go away			 				 				
			if(status[r->id0]==0) continue;
			tota=0;		

			//find boundary in two cases: 
			//a) only inserted or deleted regions
			//b) in all loop regions
							 
			Res *r1,*r2;
			for(r1=r;r1;r1=r1->last){
				if(rely[r1->id0]==2&&r1->sec!='-') break;				
				if(rely[r1->id0]>1&&sharp==2) break;						 
			}
			if(r1==0) r1=mpdb->chn->res;
			 
			for(r2=r;r2;r2=r2->next){
				if(rely[r2->id0]==2&&r2->sec!='-') break;				
				if(rely[r2->id0]>1&&sharp==2) break;				 							 
			}
			if(r2==0) r2=mpdb->chn->lastres();
			 
			Res *t;
			if(r1==0) r1=mpdb->chn->res;
			if(r2==0) r2=mpdb->chn->lastres();
	
			//check boundaries
			if(r1->id0<rr->id0) r1=rr;
			if(r2->id0>nn) r2=mpdb->chn->findsmallres(nn);
			if(r2==0) continue;

			int n1=calcinsert(r1,r2->id0);
			int n2=calcdelete(r1,r2->id0);

			//no insertion and deletion case			
			if(n1==0&&n2==0&&sharp==2) {
				for(t=r1;t;t=t->next) {
					if(t->id0>r2->id0) break;
					status[t->id0]=0;
				}	
				continue;
			}

			if(n1==0&&n2==0&&r2->id0-r1->id0+1<4) {
				status[r1->id0]=0;
				for(t=r1;t;t=t->next) {				
					if(t->id0>r2->id0) break;					
					else status[t->id0]=0;
				}
				continue;
			}	

			if(r1->last&&r2->next&&r2->id0-r1->id0+1<3) {
				status[r1->id0]=0;
				for(t=r1;t;t=t->next) {				
					if(t->id0>r2->id0) break;					
					else status[t->id0]=0;
				}
				continue;				
			}

			Res *t1,*t2;
				
			//find the most unconserved region
			t=findminscore(status,r1,r2->id0); 	
	
			if(t==0) {
				status[r1->id0]=0;
				for(t=r1;t;t=t->next) {				
					if(t->id0>r2->id0) break;					
					else status[t->id0]=0;
				}
				continue;
			}
 
			t1=r1;
			t2=r2;
			int nt1=r1->id0;
			int nt2=r2->id0;
			int nt=0;int ntx=0;

			//minimize the region of length wsize
			//for(nt=nt1;nt<=nt2;nt+=(wsize/2+0.7)) {
			//
			wsize=5;
			cerr<<"minimize regions between: "<<r1->name<<r1->id0<<"---"<<r2->name<<r2->id0<<endl;
			while(1) {  
				int ns=0;
				int y1=0;
				for(y1=r1->id0;y1<r2->id0;y1++)  if(status[y1]==1) ns++;
				if(ns==0&&wsize<r2->id0-r1->id0&&wsize<7) {
					wsize+=4;
					for(y1=r1->id0;y1<r2->id0;y1++) status[y1]=1;
				}
				
				t=findminscore(status,r1,r2->id0);
				if(t==0) break;
				if(wsize+1<r2->id0-r1->id0+1) {
					nt=t->id0-(wsize+1)/2;
					if(nt<nt1) nt=nt1;
					ntx=t->id0+(wsize+1)/2;
					if(ntx>nt2) ntx=nt2;	
				}
				else {
					nt=r1->id0;
					ntx=r2->id0;
				}
				Res *tsave=t;										
			//  			
				t1=mpdb->chn->isres0(nt);	
				/*			
				int ntx=nt+wsize;
				if(ntx>=nt2) {
					ntx=nt2;
					nt=nt2;
				}
				*/				
				t2=mpdb->chn->isres0(ntx);

				if(t1==0) t1=r1;
				if(t2==0) t2=r2;

				int n1;n1=calcinsert(t1,t2->id0);
				int n2;n2=calcdelete(t1,t2->id0);
 
				int no;no=0;

				for(t=t1;t;t=t->next) {
					if(t->id0>t2->id0) break;
					if(status[t->id0]==1) {
						no=1;break;
					}
				}
				
				//all have been done!
				if(no==0) goto re200;
				if(mpdb->chn->isalllinked(t1->id0-1,t2->id0+1)==0) {
					cerr<<"chain breakers found in the segment: "<<t1->name<<t1->id0<<"---"<<t2->name<<t2->id0<<endl;
					cerr<<"ignore this segment..."<<endl;
					goto re200;
				}
				float coeff;coeff=1;

				if(fapr==5) 	 coeff=0.25;
				else if(fapr==4) coeff=0.5;
				else if(fapr==3) coeff=1;
				else if(fapr==2) coeff=2;
				else if(fapr==1) coeff=4;
				else if(fapr==0) coeff=8;
				
				//
				coeff=2*coeff;
				//
				if(t2->id0-t1->id0+1>=10) {					
					segen->arbt=(int)(200*coeff);
				}
				else if(t2->id0-t1->id0+1>=8) {
					segen->arbt=(int)(100*coeff);
				}
				else if(t2->id0-t1->id0+1>=6) {
					segen->arbt=(int)(50*coeff);
				}
				else if(t2->id0-t1->id0+1>=4) {
					segen->arbt=(int)(30*coeff);
				}
				else{
					segen->arbt=(int)(20*coeff);
				} 			
				if(n1==0) segen->arbt=segen->arbt/2;
				
				if(test) segen->arbt=5; 

				//if((t2->id0-t1->id0+1-n1)>6&&aloop==0) goto re200;
  
				int *stem;stem=new int[2];
				stem[0]=compare[t1->id0];
				stem[1]=compare[t2->id0];
				float *xyzorg;xyzorg=mpdb->chn->gettransfer(t1,t2->id0);
				if(xyzself) delete [] xyzself;xyzself=0;
				xyzself=getpdboldtransfer(t1,t2);
				 
				tote++;tota++;
				segen->updatexyzsave(stem);
				updatepdbcopy();				 
				//delsidechain(stem);
				mpdb->chn->dihedral();
				if(TRES.logg>3) printsegment(stem);	
 	
				segen->readyminfix(t1,t2->id0);	

				//minimize sidechain only
				segen->onlysidechain=1;
				float d;d=segen->predtmin();
				//mpdb->chn->transfer(xyzorg,t1,t2->id0,0);

 				//minimize backbone
				segen->onlysidechain=0;
				segen->predtmin();
				float *xyzout;xyzout=mpdb->chn->gettransfer(t1,t2->id0);
				 
				if(stem) delete [] stem;stem=0;		
				if(xyzout==0){					 
					segen->hooksidechain(); 
					mpdb->configure();					
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1); 
					segen->setallnear();					
				}						
				if(xyzout) {	
					/*
					Res *ti=mpdb->chn->isres0(segen->start-3);
					if(ti==0) ti=mpdb->chn->res;	
					int endo=segen->end+3;	 
					mpdb->chn->transfer(xyzout,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d1=segen->disc->clash(ti,endo);
					//
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d2=segen->disc->clash(ti,endo);
						
					if(d1<d2) mpdb->chn->transfer(xyzout,t1,t2->id0,1);
					*/
					mpdb->chn->transfer(xyzout,t1,t2->id0,1);
				}
				else {
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);				 
				}

				if(xyzout) delete [] xyzout;xyzout=0;
				if(xyzorg) delete [] xyzorg;xyzorg=0;

				re200:
				if(tsave) status[tsave->id0]=0;
			 	//int nmid=(t1->id0+t2->id0)/2;				
				for(t=t1;t&&tota;t=t->next) {				
					if(t->id0>t2->id0) break;
					//else if(fabs(t->id0-nmid)<=(wsize/2-1)) status[t->id0]=0;
					//else if(t->id-t1->id0<2)  continue;
					//else if(t2->id0-t->id0<2) continue;
					else status[t->id0]=0;
				}	
				 
			}	 							
		} 

		if(tote==0) break;
	}

 

	if(status) delete [] status;status=0;
}


float *Mutate::getpdboldtransfer(Res *t1, Res *t2){

	if(pdbold==0) return 0;

	Res *r1=pdbold->chn->isres0(t1->id0);
	Res *r2=pdbold->chn->isres0(t2->id0);

	if(r1==0||r2==0||r1->name!=t1->name||r2->name!=t2->name) {
		cerr<<"pdb copy is not identical to the original"<<endl;
		return 0;
	}	

	return pdbold->chn->gettransfer(r1,t2->id0);

}
void Mutate::minwinrefine(Res *rr,int nn) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
        segen->flex=1;
        segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
        segen->randcoil=0;
	//TRES.smoothclash=1;
        segen->smoothclash=1;
        segen->part=0;
	segen->randcoil=2;
 
	segen->cutoff=8.;
	 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	mpdb->setlastresorder();
	mpdb->chn->dihedral((FILE*)0);
	int nlen=strlen(sqnto);
        int *status=new int[nlen];
		
	setnewunchanged(); //2, old, 1 move, 0 insert
	 
	Res *r;//,*t;
	int i;
	for(i=0;i<nlen;i++) status[i]=1;
	 
	refine=1000; //refinement
	int wsize=3;
	while(1) {
		int tote=0; int tota=0;
		wsize++;
		if(wsize>6) break;
		for(i=0;i<nlen;i++) status[i]=1;
		for(r=rr;r;r=r->next) {
			if(r->id0>nn) break;
			tota=0;
			//if(status[r->id0]==0) continue;
			if(rely[r->id0]==2) continue;
			//if(r->sec!='-') continue;
			//if(rely[r->id0]==1) continue;
			Res *r1,*r2;
			for(r1=r;r1;r1=r1->last){
				//if(r1->sec!='-') break;
				//if(status[r1->id0]==0) break;
				if(rely[r1->id0]>1) break;			 
			}
			if(r1==0) r1=mpdb->chn->res;
			//else	  r1=r1->next;
			for(r2=r;r2;r2=r2->next){
				//if(r2->sec!='-') break;
				//if(status[r2->id0]==0) break;	
				if(rely[r2->id0]>1) break;						 
			}
			if(r2==0) r2=mpdb->chn->lastres();
			//else      r2=r2->last;
			Res *t;
			if(r1==0) r1=mpdb->chn->res;
			if(r2==0) r2=mpdb->chn->lastres();		
			int n1=calcinsert(r1,r2->id0);
			int n2=calcdelete(r1,r2->id0);			
			if(n1==0&&n2==0) {
				for(t=r1;t;t=t->next) {
					if(t->id0>r2->id0) break;
					status[t->id0]=0;
				}	
				continue;
			}
		
			
			Res *t1,*t2;
				
			t=findminscore(status,r1,r2->id0); 		
			if(t==0) goto re200;
			
			t1=r1;
			t2=r2;
			int nt1;nt1=r1->id0;
			int nt2;nt2=r2->id0;
			int nt;nt=0;
			for(nt=nt1;nt<=nt2;nt++) { 	
				 		
				t1=mpdb->chn->isres0(nt);				
				int ntx; ntx=nt+wsize;
				if(ntx>=nt2) {
					ntx=nt2;
					nt=nt2;
				}				
				t2=mpdb->chn->isres0(ntx);

				if(t1==0) t1=r1;
				if(t2==0) t2=r2;

				int n1;n1=calcinsert(t1,t2->id0);
				int n2;n2=calcdelete(t1,t2->id0);
			
				//if(n1==0&&n2==0) goto re200;			
			 
				//if(t2->id0-t1->id0+1<5&&n1==0) {				
				//	if(t2->next) t2=t2->next;
				//	if(t1->last) t1=t1->last;				
				//}

				//if(t2->id0-t1->id0+1<5&&n2==0) {
				//	if(t2->next) t2=t2->next;
				//	if(t1->last) t1=t1->last;				
				//}			 

				//if(t2->id0-t1->id0+1<5) goto re200;	 
			
				int no;no=0;

				for(t=t1;t;t=t->next) {
					if(t->id0>t2->id0) break;
					if(status[t->id0]==1) {
						no=1;break;
					}
				}

				if(no==0) goto re200;

				if(t2->id0-t1->id0+1>=10) {
					segen->arbt=200;
				}
				else if(t2->id0-t1->id0+1>=8) {
					segen->arbt=100;
				}
				else if(t2->id0-t1->id0+1>=6) {
					segen->arbt=50;
				}
				else if(t2->id0-t1->id0+1>=4) {
					segen->arbt=30;
				}
				else{
					segen->arbt=20;
				}
 			
				if(n1==0) segen->arbt=segen->arbt/2;
				if(test) segen->arbt=5; 
				//if((t2->id0-t1->id0+1-n1)>6&&aloop==0) goto re200;
  
				int *stem;stem=new int[2];
				stem[0]=compare[t1->id0];
				stem[1]=compare[t2->id0];
				float *xyzorg;xyzorg=mpdb->chn->gettransfer(t1,t2->id0);
				if(xyzself) delete [] xyzself;xyzself=0;
				xyzself=mpdb->chn->gettransfer(t1,t2->id0);
				tote++;tota++;
				mpdb->chn->dihedral();
				if(TRES.logg) printsegment(stem);
				
				float *xyzout;xyzout=segen->myfixsegment(stem);		
				if(xyzout) {				 
					mpdb->chn->transfer(xyzout,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d1=segen->disc->clash(mpdb->chn->res,100000);
					//
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d2=segen->disc->clash(mpdb->chn->res,100000);	
					if(d1<d2) mpdb->chn->transfer(xyzout,t1,t2->id0,1);
				}
				else {
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);				 
				}

				if(xyzout) delete [] xyzout;xyzout=0;
				if(xyzorg) delete [] xyzorg;xyzorg=0;

				re200:
			 
				for(t=t1;t&&tota;t=t->next) {
					if(t->id0>t2->id0) break;
					status[t->id0]=0;
				}	
			}	 							
		} 

		if(tote==0) break;
	}

	 
	segen->mutate=0;

	if(status) delete [] status;status=0;
} 

void Mutate::setrelydssp() {

	if(rely) delete [] rely;rely=0;

        //int n=mpdb->maxresid0()+100;

        int n=strlen(sqnto);

        if(rely==0) rely=new int[n];

        int i;

        for(i=0;i<n;i++) rely[i]=0;

        Res *r;

        for(r=mpdb->chn->res;r;r=r->next) {		 
		//int n=compare[r->id0];
		if(r->sec!='-') rely[r->id0]=2;
		//else if(dssp[n]!='-')rely[r->id0]=1;
		else  rely[r->id0]=0;
        }

}

int* Mutate::getdsspmark() {

	int *rely0;
	 
        int n=strlen(sqnto);

        rely0=new int[n];

        int i;

        for(i=0;i<n;i++) rely0[i]=0;

        Res *r;

        for(r=mpdb->chn->res;r;r=r->next) {
		int n=compare[r->id0];
		if(dssp[n]=='h')rely0[r->id0]=2;
		else  if(dssp[n]=='e')rely0[r->id0]=1;
		else  if(dssp[n]=='-')rely0[r->id0]=0;
        }

	return rely0;
}

void Mutate::refinesec(Res *rr,int nn) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
        segen->flex=1;
        segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
        segen->randcoil=0;
        segen->smoothclash=0;
        segen->part=0;
	segen->randcoil=2;
	segen->cutoff=6.;
	segen->onlyenergy=0;
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	mpdb->header();
	refine=1000; //secondary refinement 

	mpdb->setlastresorder();
	mpdb->chn->dihedral((FILE*)0);
	int nlen=strlen(sqnto);
        int *status=new int[nlen];
 
	//set rely value, consistent secondary is 2 otherwise 0;
	setrelydssp();
 	int *rely0=getdsspmark();

	Res *r; 

	//status indicates if minimized.
	int i;
	for(i=0;i<nlen;i++) status[i]=1;
 	
	int ih;
        for(ih=0;ih<1;ih++) {				
		int tote=0; int tota=0;
		
		for(i=0;i<nlen;i++) status[i]=1;
		for(r=rr;r;r=r->next) {
			if(r->id0>nn) break;			
			if(rely0[r->id0]==0) continue; //not secondary		
			if(status[r->id0]==0) continue; //done before
			tota=0;
		
			Res *r1,*r2,*v1,*v2;
			for(r1=r;r1;r1=r1->last){
				if(rely0[r1->id0]!=rely0[r->id0]) break;			
				//if(r1->sec!=r->sec) break;			 
			}
			if(r1==0) r1=mpdb->chn->res;
			if(rr->id0-r1->id0>5) continue;	//check boundaries
			v1=r1;
			//find boundary of v1			
			if(rely0[v1->id0]==0) {
				while(v1->last&&rely0[v1->last->id0]==0&&r1->id0-v1->id0<3) v1=v1->last;
			}   
 			
			
			for(r2=r;r2;r2=r2->next){
				if(rely0[r2->id0]!=rely0[r->id0]) break;				 								 
			}
			if(r2==0) r2=mpdb->chn->lastres();		 
			if(r2->id0-nn>5) continue; //check boundaries	
			v2=r2;
			//find boundary of v2			
			if(rely0[v2->id0]==0) {
				while(v2->next&&rely0[v2->next->id0]==0&&v2->id-r2->id0<3) v2=v2->next;
			}
 
			Res *t;
 
			//define refinement type
			if(rely0[r->id0]==2) refine=1000; 
			else if(rely0[r->id0]==1) refine=2000;  
			else {
				refine=2000;
				status[r->id0]=0;
				continue;
			}
  
			//set fixed torsions
			if(segen->rotatm) delete [] segen->rotatm;segen->rotatm=0;
			/*
			int ntot=mpdb->manyatm()+100;
			if(segen->rotatm) delete [] segen->rotatm;segen->rotatm=0;
			segen->rotatm=new Atm*[ntot];
			
			int nn=0;
			for(nn=0;nn<ntot;nn++) segen->rotatm[nn]=0;
			nn=0;
			for(t=r1;t;t=t->next) {
				if(t->id0>r2->id0) break;
				if(t->sec=='-') continue;
				if(t->sec!=r->sec) continue;
				segen->rotatm[nn++]=t->isatmid(1);
				segen->rotatm[nn++]=t->isatmid(2);		
			} 			 			

			int i=0;
			for(nn=0;nn<ntot;nn++) {
				if(segen->rotatm[nn]==0) continue;
				segen->rotatm[i++]=segen->rotatm[nn];
			}
			segen->rotatm[i]=0;
			*/

 			//minimize secondary now
			cerr<<"refine regions between: "<<v1->name<<v1->id0<<"---"<<v2->name<<v2->id0<<endl;
			Res *t1,*t2;
			int nt1,nt2;
			for(nt1=v1->id0;nt1<=r1->id0;nt1+=20) 	
			for(nt2=v2->id0;nt2>=r2->id0;nt2-=20) { 	
			 			
				t1=mpdb->chn->isres0(nt1);								
				t2=mpdb->chn->isres0(nt2);

				if(t1==0) t1=v1;
				if(t2==0) t2=v2;

				int n1=calcinsert(t1,t2->id0);
				int n2=calcdelete(t1,t2->id0);
 
				int no=0;

				for(t=t1;t;t=t->next) {
					if(t->id0>t2->id0) break;
					if(status[t->id0]==1) {
						no=1;break;
					}
				}

				if(no==0) continue;
				if(mpdb->chn->isalllinked(t1->id0-1,t2->id0+1)==0) {
					cerr<<"chain breakers found in the segment: "<<t1->name<<t1->id0<<"---"<<t2->name<<t2->id0<<endl;
					cerr<<"ignore this segment..."<<endl;
					continue;
				}
				float coeff=1;

                                if(fapr==5)      coeff=0.25;
                                else if(fapr==4) coeff=0.5;
                                else if(fapr==3) coeff=1;
                                else if(fapr==2) coeff=2;
                                else if(fapr==1) coeff=4;
                                else if(fapr==0) coeff=8;

				//
                                coeff=2*coeff;
                                //
                                if(t2->id0-t1->id0+1>=10) {
                                        segen->arbt=(int)(200*coeff);
                                }
                                else if(t2->id0-t1->id0+1>=8) {
                                        segen->arbt=(int)(100*coeff);
                                }
                                else if(t2->id0-t1->id0+1>=6) {
                                        segen->arbt=(int)(50*coeff);
                                }
                                else if(t2->id0-t1->id0+1>=4) {
                                        segen->arbt=(int)(30*coeff);
                                }
                                else{
                                        segen->arbt=(int)(20*coeff);
                                }
                                if(n1==0) segen->arbt=segen->arbt/2;

                                if(test) segen->arbt=5;
							 
  
				int *stem=new int[2];
				stem[0]=compare[t1->id0];
				stem[1]=compare[t2->id0];
				float *xyzorg=mpdb->chn->gettransfer(t1,t2->id0);
				if(xyzself) delete [] xyzself;xyzself=0;
				xyzself=getpdboldtransfer(t1,t2);
				 
				tote++;tota++;
				segen->updatexyzsave(stem);
				updatepdbcopy();
				delsidechain(stem);
				mpdb->chn->dihedral();
				if(TRES.logg>3) printsegment(stem);
				
				float *xyzout=segen->myfixsegment(stem);
				
				if(stem) delete [] stem;stem=0;	
				if(xyzout==0){
                                        segen->hooksidechain();
                                        mpdb->configure();
                                        mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
                                        segen->setallnear();
                                }
	
				if(xyzout) {				 
					Res *ti=mpdb->chn->isres0(segen->start-3);
                                        if(ti==0) ti=mpdb->chn->res;
                                        int endo=segen->end+3;
                                        mpdb->chn->transfer(xyzout,t1,t2->id0,1);
                                        segen->disc->setupall(mpdb,10,1);
                                        strcpy(segen->disc->force,"uDd");
                                        float d1=segen->disc->clash(ti,endo);
                                        //
                                        mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
                                        segen->disc->setupall(mpdb,10,1);
                                        strcpy(segen->disc->force,"uDd");
                                        float d2=segen->disc->clash(ti,endo);
                                        if(d1<d2) {
						mpdb->chn->transfer(xyzout,t1,t2->id0,1);
						mpdb->chn->header(t1->id0-3,t2->id0+3);
					}
				}
				
				if(xyzout) delete [] xyzout;xyzout=0;
				if(xyzorg) delete [] xyzorg;xyzorg=0;
				
			}
			status[r->id0]=0;
			for(t=v1;t&&tota;t=t->next) {
				if(t->id0>v2->id0) break;				 
				status[t->id0]=0;
			}	
			//goto re400; 	 							
		} 

		if(tote==0) break;
	}

 
	if(rely0) delete [] rely0;rely0=0;
	if(status) delete [] status;status=0;
} 


void Mutate::refineminsec(Res *rr,int nn) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
        segen->flex=1;
        segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
        segen->randcoil=0;
        segen->smoothclash=0;
        segen->part=0;
	segen->randcoil=2;
	segen->cutoff=6.;
	segen->onlyenergy=0;
 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	refine=1000; //secondary refinement 

	mpdb->setlastresorder();
	mpdb->chn->dihedral((FILE*)0);
	int nlen=strlen(sqnto);
        int *status=new int[nlen];
 
	//set rely value, consistent secondary is 2 otherwise 0;
	setrelydssp();
 	int *rely0=getdsspmark();

	Res *r; 

	//status indicates if minimized.
	int i;
	for(i=0;i<nlen;i++) status[i]=1;
 	
	int ih;
        for(ih=0;ih<1;ih++) {				
		int tote=0; int tota=0;
		
		for(i=0;i<nlen;i++) status[i]=1;
		for(r=rr;r;r=r->next) {
			if(r->id0>nn) break;			
			if(rely0[r->id0]==0) continue; //not secondary		
			if(status[r->id0]==0) continue; //done before
			tota=0;
		
			Res *r1,*r2,*v1,*v2;
			for(r1=r;r1;r1=r1->last){
				if(rely0[r1->id0]!=rely0[r->id0]) break;			
				//if(r1->sec!=r->sec) break;			 
			}
			if(r1==0) r1=mpdb->chn->res;
			if(rr->id0-r1->id0>5) continue; //check boundaries

			v1=r1;
			//find boundary of v1			
			if(rely0[v1->id0]==0) {
				while(v1->last&&rely0[v1->last->id0]==0&&r1->id0-v1->id0<3) v1=v1->last;
			}   
 			
			
			for(r2=r;r2;r2=r2->next){
				if(rely0[r2->id0]!=rely0[r->id0]) break;				 								 
			}
			if(r2==0) r2=mpdb->chn->lastres();		 
			if(r2->id0-nn>5) continue; //check boundaries

			v2=r2;
			//find boundary of v2			
			if(rely0[v2->id0]==0) {
				while(v2->next&&rely0[v2->next->id0]==0&&v2->id-r2->id0<3) v2=v2->next;
			}
 
			Res *t;
 
			//define refinement type
			if(rely0[r->id0]==2) refine=1000; 
			else if(rely0[r->id0]==1) refine=2000;  
			else {
				refine=2000;
				status[r->id0]=0;
				continue;
			}
  
			//set fixed torsions
			if(segen->rotatm) delete [] segen->rotatm;segen->rotatm=0;
			/*
			int ntot=mpdb->manyatm()+100;
			if(segen->rotatm) delete [] segen->rotatm;segen->rotatm=0;
			segen->rotatm=new Atm*[ntot];
			
			int nn=0;
			for(nn=0;nn<ntot;nn++) segen->rotatm[nn]=0;
			nn=0;
			for(t=r1;t;t=t->next) {
				if(t->id0>r2->id0) break;
				if(t->sec=='-') continue;
				if(t->sec!=r->sec) continue;
				segen->rotatm[nn++]=t->isatmid(1);
				segen->rotatm[nn++]=t->isatmid(2);		
			} 			 			

			int i=0;
			for(nn=0;nn<ntot;nn++) {
				if(segen->rotatm[nn]==0) continue;
				segen->rotatm[i++]=segen->rotatm[nn];
			}
			segen->rotatm[i]=0;
			*/
			cerr<<"minimize regions between: "<<v1->name<<v1->id0<<"---"<<v2->name<<v2->id0<<endl;
 			//minimize secondary now
			Res *t1,*t2;
			int nt1,nt2;
			for(nt1=v1->id0;nt1<=r1->id0;nt1+=2) 	
			for(nt2=v2->id0;nt2>=r2->id0;nt2-=2) { 	
			 			
				t1=mpdb->chn->isres0(nt1);								
				t2=mpdb->chn->isres0(nt2);

				if(t1==0) t1=v1;
				if(t2==0) t2=v2;

				int n1=calcinsert(t1,t2->id0);
				int n2=calcdelete(t1,t2->id0);
 
				int no=0;

				for(t=t1;t;t=t->next) {
					if(t->id0>t2->id0) break;
					if(status[t->id0]==1) {
						no=1;break;
					}
				}

				if(no==0) continue;
				if(mpdb->chn->isalllinked(t1->id0-1,t2->id0+1)==0) {
					cerr<<"chain breakers found: "<<t1->name<<t1->id0<<"---"<<t2->name<<t2->id0<<endl;
					cerr<<"ignore this segment..."<<endl;
					continue;
				}
				float coeff=1;

                                if(fapr==5)      coeff=0.25;
                                else if(fapr==4) coeff=0.5;
                                else if(fapr==3) coeff=1;
                                else if(fapr==2) coeff=2;
                                else if(fapr==1) coeff=4;
                                else if(fapr==0) coeff=8;

				//
                                coeff=2*coeff;
                                //
                                if(t2->id0-t1->id0+1>=10) {
                                        segen->arbt=(int)(200*coeff);
                                }
                                else if(t2->id0-t1->id0+1>=8) {
                                        segen->arbt=(int)(100*coeff);
                                }
                                else if(t2->id0-t1->id0+1>=6) {
                                        segen->arbt=(int)(50*coeff);
                                }
                                else if(t2->id0-t1->id0+1>=4) {
                                        segen->arbt=(int)(30*coeff);
                                }
                                else{
                                        segen->arbt=(int)(20*coeff);
                                }
                                if(n1==0) segen->arbt=segen->arbt/2;

                                if(test) segen->arbt=5;
							 
  
				int *stem=new int[2];
				stem[0]=compare[t1->id0];
				stem[1]=compare[t2->id0];
				float *xyzorg=mpdb->chn->gettransfer(t1,t2->id0);
				if(xyzself) delete [] xyzself;xyzself=0;
				xyzself=getpdboldtransfer(t1,t2);
				 
				tote++;tota++;
				segen->updatexyzsave(stem);
				updatepdbcopy();
				//delsidechain(stem);
				mpdb->chn->dihedral();
				if(TRES.logg>3) printsegment(stem);
				
				segen->readyminfix(t1,t2->id0);

                                //minimize sidechain only
                                segen->onlysidechain=1;
                                float d=segen->predtmin();
                                //mpdb->chn->transfer(xyzorg,t1,t2->id0,0);

                                //minimize backbone
                                segen->onlysidechain=0;
                                segen->predtmin();
                                float *xyzout=mpdb->chn->gettransfer(t1,t2->id0);
				
				if(stem) delete [] stem;stem=0;	
				if(xyzout==0){
                                        segen->hooksidechain();
                                        mpdb->configure();
                                        mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
                                        segen->setallnear();
                                }
	
				if(xyzout) {				 
					Res *ti=mpdb->chn->isres0(segen->start-3);
                                        if(ti==0) ti=mpdb->chn->res;
                                        int endo=segen->end+3;
                                        mpdb->chn->transfer(xyzout,t1,t2->id0,1);
                                        segen->disc->setupall(mpdb,10,1);
                                        strcpy(segen->disc->force,"uDd");
                                        float d1=segen->disc->clash(ti,endo);
                                        //
                                        mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
                                        segen->disc->setupall(mpdb,10,1);
                                        strcpy(segen->disc->force,"uDd");
                                        float d2=segen->disc->clash(ti,endo);
                                        if(d1<d2) mpdb->chn->transfer(xyzout,t1,t2->id0,1);
				}
				
				if(xyzout) delete [] xyzout;xyzout=0;
				if(xyzorg) delete [] xyzorg;xyzorg=0;
				
			}
			status[r->id0]=0;
			for(t=v1;t&&tota;t=t->next) {
				if(t->id0>v2->id0) break;				 
				status[t->id0]=0;
			}	
			//goto re400; 	 							
		} 

		if(tote==0) break;
	}

 
	if(rely0) delete [] rely0;rely0=0;
	if(status) delete [] status;status=0;
} 

void Mutate::refinesmoothsec(Res *rr,int nn) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
        segen->flex=1;
        segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
        segen->randcoil=0;
        segen->smoothclash=TRES.smoothclash;
        segen->part=0;
	segen->randcoil=2;
	segen->cutoff=6.;
	segen->onlyenergy=0;	 

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	refine=-1000; //secondary refinement 

	mpdb->setlastresorder();
	mpdb->chn->dihedral((FILE*)0);
	int nlen=strlen(sqnto);
        int *status=new int[nlen];
 
	//set rely value, consistent secondary is 2 otherwise 0;
	setrelydssp();
 	int *rely0=getdsspmark();

	Res *r; 

	//status indicates if minimized.
	int i;
	for(i=0;i<nlen;i++) status[i]=1;
 	
	int ih;
        for(ih=0;ih<1;ih++) {				
		int tote=0; int tota=0;
		
		for(i=0;i<nlen;i++) status[i]=1;
		for(r=rr;r;r=r->next) {
			if(r->id0>nn) break;			
			if(rely0[r->id0]==0) continue; //not secondary		
			if(status[r->id0]==0) continue; //done before
			tota=0;
		
			Res *r1,*r2,*v1,*v2;
			for(r1=r;r1;r1=r1->last){
				if(rely0[r1->id0]!=rely0[r->id0]) break;			
				//if(r1->sec!=r->sec) break;			 
			}
			if(r1==0) r1=mpdb->chn->res;
			if(rr->id0-r1->id0>5) continue; //check boundaries
			
			v1=r1;
			//find boundary of v1			
			if(rely0[v1->id0]==0) {
				while(v1->last&&rely0[v1->last->id0]==0&&r1->id0-v1->id0<3) v1=v1->last;
			}   
 			
			
			for(r2=r;r2;r2=r2->next){
				if(rely0[r2->id0]!=rely0[r->id0]) break;				 								 
			}
			if(r2==0) r2=mpdb->chn->lastres();		 
			if(r2->id0-nn>5) continue; //check boundaries
			
			v2=r2;
			//find boundary of v2			
			if(rely0[v2->id0]==0) {
				while(v2->next&&rely0[v2->next->id0]==0&&v2->id-r2->id0<3) v2=v2->next;
			}
 
			Res *t;
 
			//define refinement type
			if(rely0[r->id0]==2) refine=-1000; 
			else if(rely0[r->id0]==1) refine=-2000;  
			else {
				refine=-2000;
				status[r->id0]=0;
				continue;
			}
  
			//set fixed torsions
			if(segen->rotatm) delete [] segen->rotatm;segen->rotatm=0;
			/*
			int ntot=mpdb->manyatm()+100;
			if(segen->rotatm) delete [] segen->rotatm;segen->rotatm=0;
			segen->rotatm=new Atm*[ntot];
			
			int nn=0;
			for(nn=0;nn<ntot;nn++) segen->rotatm[nn]=0;
			nn=0;
			for(t=r1;t;t=t->next) {
				if(t->id0>r2->id0) break;
				if(t->sec=='-') continue;
				if(t->sec!=r->sec) continue;
				segen->rotatm[nn++]=t->isatmid(1);
				segen->rotatm[nn++]=t->isatmid(2);		
			} 			 			

			int i=0;
			for(nn=0;nn<ntot;nn++) {
				if(segen->rotatm[nn]==0) continue;
				segen->rotatm[i++]=segen->rotatm[nn];
			}
			segen->rotatm[i]=0;
			*/
			cerr<<"minimize regions between: "<<v1->name<<v1->id0<<"---"<<v2->name<<v2->id0<<endl;
 			//minimize secondary now
			Res *t1,*t2;
			int nt1,nt2;
			for(nt1=v1->id0;nt1<=r1->id0;nt1+=2) 	
			for(nt2=v2->id0;nt2>=r2->id0;nt2-=2) { 	
			 			
				t1=mpdb->chn->isres0(nt1);								
				t2=mpdb->chn->isres0(nt2);

				if(t1==0) t1=v1;
				if(t2==0) t2=v2;

				int n1=calcinsert(t1,t2->id0);
				int n2=calcdelete(t1,t2->id0);
 
				int no=0;

				for(t=t1;t;t=t->next) {
					if(t->id0>t2->id0) break;
					if(status[t->id0]==1) {
						no=1;break;
					}
				}

				if(no==0) continue;
				if(mpdb->chn->isalllinked(t1->id0-1,t2->id0+1)==0) {
					cerr<<"chain breakers found: "<<t1->name<<t1->id0<<"---"<<t2->name<<t2->id0<<endl;
					cerr<<"ignore this segment..."<<endl;
					continue;
				}
				float coeff=1;

                                if(fapr==5)      coeff=0.25;
                                else if(fapr==4) coeff=0.5;
                                else if(fapr==3) coeff=1;
                                else if(fapr==2) coeff=2;
                                else if(fapr==1) coeff=4;
                                else if(fapr==0) coeff=8;

				//
                                coeff=2*coeff;
                                //
                                if(t2->id0-t1->id0+1>=10) {
                                        segen->arbt=(int)(200*coeff);
                                }
                                else if(t2->id0-t1->id0+1>=8) {
                                        segen->arbt=(int)(100*coeff);
                                }
                                else if(t2->id0-t1->id0+1>=6) {
                                        segen->arbt=(int)(50*coeff);
                                }
                                else if(t2->id0-t1->id0+1>=4) {
                                        segen->arbt=(int)(30*coeff);
                                }
                                else{
                                        segen->arbt=(int)(20*coeff);
                                }
                                if(n1==0) segen->arbt=segen->arbt/2;

                                if(test) segen->arbt=5;
							 
  
				int *stem=new int[2];
				stem[0]=compare[t1->id0];
				stem[1]=compare[t2->id0];
				float *xyzorg=mpdb->chn->gettransfer(t1,t2->id0);
				if(xyzself) delete [] xyzself;xyzself=0;
				xyzself=getpdboldtransfer(t1,t2);
				 
				tote++;tota++;
				segen->updatexyzsave(stem);
				updatepdbcopy();
				//delsidechain(stem);
				mpdb->chn->dihedral();
				if(TRES.logg>3) printsegment(stem);
				
				segen->readyminfix(t1,t2->id0);

                                //minimize sidechain only
                                segen->onlysidechain=1;
                                float d=segen->predtmin();
                                //mpdb->chn->transfer(xyzorg,t1,t2->id0,0);

                                //minimize backbone
                                segen->onlysidechain=0;
                                segen->predtmin();
                                float *xyzout=mpdb->chn->gettransfer(t1,t2->id0);
				
				if(stem) delete [] stem;stem=0;	
				if(xyzout==0){
                                        segen->hooksidechain();
                                        mpdb->configure();
                                        mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
                                        segen->setallnear();
                                }
	
				if(xyzout) {	
					/*			 
					Res *ti=mpdb->chn->isres0(segen->start-3);
                                        if(ti==0) ti=mpdb->chn->res;
                                        int endo=segen->end+3;
                                        mpdb->chn->transfer(xyzout,t1,t2->id0,1);
                                        segen->disc->setupall(mpdb,10,1);
                                        strcpy(segen->disc->force,"uDd");
                                        float d1=segen->disc->clash(ti,endo);
                                        //
                                        mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
                                        segen->disc->setupall(mpdb,10,1);
                                        strcpy(segen->disc->force,"uDd");
                                        float d2=segen->disc->clash(ti,endo);
                                        if(d1<d2) mpdb->chn->transfer(xyzout,t1,t2->id0,1);
					*/
					mpdb->chn->transfer(xyzout,t1,t2->id0,1);
				}
				
				if(xyzout) delete [] xyzout;xyzout=0;
				if(xyzorg) delete [] xyzorg;xyzorg=0;
				
			}
			status[r->id0]=0;
			for(t=v1;t&&tota;t=t->next) {
				if(t->id0>v2->id0) break;				 
				status[t->id0]=0;
			}	
			//goto re400; 	 							
		} 

		if(tote==0) break;
	}

 
	if(rely0) delete [] rely0;rely0=0;
	if(status) delete [] status;status=0;
} 

void Mutate::minwinrefine0(Res *rr,int nn) {

	if(segen) delete segen;segen=0;
        segen=new Segen();
        segen->flex=1;
        segen->arbt=100;
        segen->pdb=mpdb;
        segen->cid=mpdb->chn->id;
        segen->mutate=this;
        segen->randcoil=0;
	//TRES.smoothclash=1;
        segen->smoothclash=1;
        segen->part=0;
	segen->randcoil=2;
 
	segen->cutoff=8.;
	 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;

	mpdb->setlastresorder();
	mpdb->chn->dihedral((FILE*)0);
	int nlen=strlen(sqnto);
        int *status=new int[nlen];
	//TRES.setnewengcoeff(3,0.5);	
	setnewunchanged();
	 
	Res *r;//,*t;
	int i;
	for(i=0;i<nlen;i++) status[i]=1;
 
	refine=1000; //refinement
	int wsize=3;
	while(1) {
		sethbonddssp(4);
		setrelydssp();
		int tote=0; int tota=0;
		
		if(wsize>6) break;
		for(i=0;i<nlen;i++) status[i]=1;
		for(r=rr;r;r=r->next) {
			if(r->id0>nn) break;
			if(status[r->id0]==0) continue;
			tota=0;
			//if(status[r->id0]==0) continue;
			//if(rely[r->id0]==2) continue;
			//if(r->sec!='-') continue;
			//if(rely[r->id0]==1) continue;
			Res *r1,*r2;
			for(r1=r;r1;r1=r1->last){
				//if(r1->sec!='-') break;
				//if(status[r1->id0]==0) break;
				//if(rely[r1->id0]>1) break;
				if(r1->sec!=r->sec) break;			 
			}
			if(r1==0) r1=mpdb->chn->res;
			//else	  r1=r1->next;
			for(r2=r;r2;r2=r2->next){
				//if(r2->sec!='-') break;
				//if(status[r2->id0]==0) break;	
				//if(rely[r2->id0]>1) break;	
				if(r2->sec!=r->sec) break;					 
			}
			if(r2==0) r2=mpdb->chn->lastres();
			//else      r2=r2->last;
			
			if(r->sec=='-') {
				if(r1->last) {
					if(r1->sec!='-') r1=r1->next;
				}

				if(r2->next) {
					if(r2->sec!='-') r2=r2->last;
				}
			}
 			else {
				goto re300;
			}

			Res *t;
			if(r1==0) r1=mpdb->chn->res;
			if(r2==0) r2=mpdb->chn->lastres();		
			int n1;n1=calcinsert(r1,r2->id0);
			int n2;n2=calcdelete(r1,r2->id0);			
			if(n1==0&&n2==0) {
				for(t=r1;t;t=t->next) {
					if(t->id0>r2->id0) break;
					status[t->id0]=0;
				}	
				continue;
			}
		
			
			Res *t1,*t2;
				
			t=findminscore(status,r1,r2->id0); 		
			if(t==0) goto re200;
			
			t1=r1;
			t2=r2;
			int nt1;nt1=r1->id0;
			int nt2;nt2=r2->id0;
			int nt;nt=0;
			int wsize;wsize=3;
			re100:

			refine=100; //think as loopy
			wsize+=2;
			if(wsize==7) wsize=6;
			if(wsize>6) continue;
			for(nt=nt1;nt<=nt2-wsize;nt++) { 	
				 		
				t1=mpdb->chn->isres0(nt);				
				int ntx;ntx=nt+wsize;
				if(ntx>=nt2) {
					ntx=nt2;
					nt=nt2;
				}				
				t2=mpdb->chn->isres0(ntx);

				if(t1==0) t1=r1;
				if(t2==0) t2=r2;

				int n1;n1=calcinsert(t1,t2->id0);
				int n2;n2=calcdelete(t1,t2->id0);
			
				//if(n1==0&&n2==0) goto re200;			
			 
				//if(t2->id0-t1->id0+1<5&&n1==0) {				
				//	if(t2->next) t2=t2->next;
				//	if(t1->last) t1=t1->last;				
				//}

				//if(t2->id0-t1->id0+1<5&&n2==0) {
				//	if(t2->next) t2=t2->next;
				//	if(t1->last) t1=t1->last;				
				//}			 

				//if(t2->id0-t1->id0+1<5) goto re200;	 
			
				int no;no=0;

				for(t=t1;t;t=t->next) {
					if(t->id0>t2->id0) break;
					if(status[t->id0]==1) {
						no=1;break;
					}
				}

				if(no==0) goto re200;

				if(t2->id0-t1->id0+1>=10) {
					segen->arbt=200;
				}
				else if(t2->id0-t1->id0+1>=8) {
					segen->arbt=100;
				}
				else if(t2->id0-t1->id0+1>=6) {
					segen->arbt=50;
				}
				else if(t2->id0-t1->id0+1>=4) {
					segen->arbt=30;
				}
				else{
					segen->arbt=20;
				}
 			
				if(n1==0) segen->arbt=segen->arbt/2;
				if(test) segen->arbt=5; 
				//if((t2->id0-t1->id0+1-n1)>6&&aloop==0) goto re200;
  
				int *stem;stem=new int[2];
				stem[0]=compare[t1->id0];
				stem[1]=compare[t2->id0];
				float *xyzorg;xyzorg=mpdb->chn->gettransfer(t1,t2->id0);
				if(xyzself) delete [] xyzself;xyzself=0;
				xyzself=mpdb->chn->gettransfer(t1,t2->id0);
				tote++;tota++;
				mpdb->chn->dihedral();
				if(TRES.logg) printsegment(stem);
				
				float *xyzout;xyzout=segen->myfixsegment(stem);		
				if(xyzout) {				 
					mpdb->chn->transfer(xyzout,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d1=segen->disc->clash(mpdb->chn->res,100000);
					//
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d2=segen->disc->clash(mpdb->chn->res,100000);	
					if(d1<d2) mpdb->chn->transfer(xyzout,t1,t2->id0,1);
				}
				else {
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);				 
				}

				if(xyzout) delete [] xyzout;xyzout=0;
				if(xyzorg) delete [] xyzorg;xyzorg=0;

				re200:
			 
				for(t=t1;t&&tota;t=t->next) {
					if(t->id0>t2->id0) break;
					status[t->id0]=0;
				}	
			}
			 
			goto re100;

			re300:

			if(r->sec=='h') refine=1000; //treated as refinement
			else if(r->sec=='e') refine=2000;
			else refine=2000;

			Res *v1,*v2;
			for(v1=r1;v1;v1=v1->last) {			
				if(v1->sec!='-') break;
			}

			for(v2=r2;v2;v2=v2->next) {			
				if(v2->sec!='-') break;
			}

			if(v1==0) v1=mpdb->chn->res;
			if(v2==0) v2=mpdb->chn->lastres();

			t1=r1;
			t2=r2;
			nt1=v1->id0;
			nt2=v2->id0;
			nt=0;
			wsize=8;
			//re400:
			//wsize++;
			//if(wsize>16) continue;
			
			if(segen->rotatm) delete [] segen->rotatm;segen->rotatm=0;
			segen->rotatm=new Atm*[3000];
			
			int nn=0;
			for(nn=0;nn<3000;nn++) segen->rotatm[nn]=0;
			nn=0;
			for(t=r1;t;t=t->next) {
				if(t->id0>r2->id0) break;
				if(t->sec=='-') continue;
				segen->rotatm[nn++]=t->isatmid(1);
				segen->rotatm[nn++]=t->isatmid(2);		
			} 			 			

			for(nt1=v1->id0;nt1<=r1->id0;nt1+=2) 
			for(nt2=r2->id0;nt2<=v2->id0;nt2+=2) { 	
			 			
				t1=mpdb->chn->isres0(nt1);				
				//int ntx=nt+wsize;
				//if(ntx>=nt2) {
				//	ntx=nt2;
				//	nt=nt2;
				//}				
				t2=mpdb->chn->isres0(nt2);

				if(t1==0) t1=r1;
				if(t2==0) t2=r2;

				int n1=calcinsert(t1,t2->id0);
				int n2=calcdelete(t1,t2->id0);
			
				//if(n1==0&&n2==0) goto re200;			
			 
				//if(t2->id0-t1->id0+1<5&&n1==0) {				
				//	if(t2->next) t2=t2->next;
				//	if(t1->last) t1=t1->last;				
				//}

				//if(t2->id0-t1->id0+1<5&&n2==0) {
				//	if(t2->next) t2=t2->next;
				//	if(t1->last) t1=t1->last;				
				//}			 

				//if(t2->id0-t1->id0+1<5) goto re2000;	 
			
				int no=0;

				for(t=t1;t;t=t->next) {
					if(t->id0>t2->id0) break;
					if(status[t->id0]==1) {
						no=1;break;
					}
				}

				if(no==0) goto re2000;

				if(t2->id0-t1->id0+1>=10) {
					segen->arbt=200;
				}
				else if(t2->id0-t1->id0+1>=8) {
					segen->arbt=100;
				}
				else if(t2->id0-t1->id0+1>=6) {
					segen->arbt=50;
				}
				else if(t2->id0-t1->id0+1>=4) {
					segen->arbt=30;
				}
				else{
					segen->arbt=20;
				}
 				
				if(n1==0) segen->arbt=segen->arbt/2;
				if(test) segen->arbt=10; 
				//if((t2->id0-t1->id0+1-n1)>6&&aloop==0) goto re2000;
  
				int *stem;stem=new int[2];
				stem[0]=compare[t1->id0];
				stem[1]=compare[t2->id0];
				float *xyzorg;xyzorg=mpdb->chn->gettransfer(t1,t2->id0);
				if(xyzself) delete [] xyzself;xyzself=0;
				xyzself=mpdb->chn->gettransfer(t1,t2->id0);
				tote++;tota++;
				mpdb->chn->dihedral();
				if(TRES.logg) printsegment(stem);
				
				float *xyzout;xyzout=segen->myfixsegment(stem);
				//segen->readyminfix(t1,t2->id0);
				//float d=segen->predtmin();
				//float *xyzout=mpdb->chn->gettransfer(t1,t2->id0);
				if(stem) delete [] stem;stem=0;		
				if(xyzout) {				 
					mpdb->chn->transfer(xyzout,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d1=segen->disc->clash(mpdb->chn->res,100000);
					//
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);
					segen->disc->setupall(mpdb,10,1);
					strcpy(segen->disc->force,"uDd");
					float d2=segen->disc->clash(mpdb->chn->res,100000);	
					if(d1<d2) mpdb->chn->transfer(xyzout,t1,t2->id0,1);
				}
				else {
					mpdb->chn->transfer(xyzorg,t1,t2->id0,1);				 
				}

				if(xyzout) delete [] xyzout;xyzout=0;
				if(xyzorg) delete [] xyzorg;xyzorg=0;

				re2000:
			 
				for(t=t1;t&&tota;t=t->next) {
					if(t->id0>t2->id0) break;
					if(t->id0>=r2->id0) continue;
					if(t->id0<=r1->id0) continue;
					status[t->id0]=0;
				}	
			}
			//goto re400; 	 							
		} 

		if(tote==0) break;
	}

 

	if(status) delete [] status;status=0;
} 

void Mutate::minsidechain(Res *rr,int nn) {

	mpdb->setflgr(-99999);
	
	Res *r;
	for(r=rr;r;r=r->next) {
		if(r->id0>nn) break;
		if(strchr("PGA",r->name)) continue;
		r->flag=10000;
	}
	 	
	Scap scprd; 
	//08/16/2002
	strcpy(scprd.force,"124");
	TRES.setdonaldcharge();
	scprd.dielectric=15;
	//
	scprd.singletorsion=1;
	scprd.colonyline=1;
	scprd.colony=2;
	scprd.ncolony=1;
	scprd.nncolony=1;
	scprd.nummore=100;
	scprd.bmax=2;
	scprd.tormax=2;
	scprd.ring=1;
	scprd.pdb=mpdb;
	scprd.includeself=1;		
	scprd.resultout=0;

	if(fapr==5) scprd.nummore=10;
	else if(fapr==4) scprd.nummore=20;
	else if(fapr==3) scprd.nummore=35;
	else if(fapr==2) scprd.nummore=50;	
	else if(fapr==1) scprd.nummore=70;
	else if(fapr==0) scprd.nummore=100;

	setrotamerorder(0); 
 		
	scprd.scpred(1);
	scprd.pdb=0;
	for(Chn *c=mpdb->chn;c;c=c->next) c->setfreenextonmore();

	setrotamerorder(1);
	//
	TRES.switchcharge(1);
	//
	if(TRES.logg>3) mpdb->write("side");
}

void Mutate::minbadsidechain(Res *rr,int nn,int ff) {
//ff =0, using mix rotamer
//ff =1, using small rotamer
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();

	mpdb->setflgr(-99999);
	
	Res *r;

	cerr<<endl<<"find conserved residues..."<<endl;
	for(r=rr;r;r=r->next) {
                if(r->id0>nn) break;
                Res *rr=getres(r);
                if(rr&&r->name==rr->name) {
                        cerr<<"conserved residues ignored in scap: "<<r->id+se->start<<","<<r->name<<endl;
                }
        }

	cerr<<endl<<"find unconserved residues..."<<endl;
	for(r=rr;r;r=r->next) {
		if(r->id0>nn) break;
		if(strchr("PGA",r->name)) continue;
		Res *rr=getres(r);
		if(rr&&r->name==rr->name) continue;
		r->flag=10000;
		cerr<<"unconserved residues will be minimized with scap: "<<r->id+se->start<<","<<r->name<<endl;
	}
	 	
	
	Scap scprd; 
	//08/16/2002
	strcpy(scprd.force,"124");
	TRES.setdonaldcharge();
	scprd.dielectric=15;
	//end
	scprd.singletorsion=1;
	scprd.colonyline=1;
	scprd.colony=2;
	scprd.ncolony=1;
	scprd.nncolony=1;
	scprd.nummore=100;
	scprd.bmax=2;
	scprd.tormax=2;
	scprd.ring=1;
	scprd.pdb=mpdb;
	scprd.includeself=1;		
	scprd.resultout=0;

	if(fapr==5) scprd.nummore=10;
	else if(fapr==4) scprd.nummore=20;
	else if(fapr==3) scprd.nummore=35;
	else if(fapr==2) scprd.nummore=50;	
	else if(fapr==1) scprd.nummore=70;
	else if(fapr==0) scprd.nummore=100;
	
	setrotamerorder(ff); 
 		
	scprd.scpred(1);
	scprd.pdb=0;
	for(Chn *c=mpdb->chn;c;c=c->next) c->setfreenextonmore();

	setrotamerorder(1);
	TRES.switchcharge(1);
	if(TRES.logg>3) mpdb->write("side");
}
void Mutate::setrotamerorder(int n) {

	Pdb *s1=TRES.findrotamername("side_small_rotamer"); 
	Pdb *s2=TRES.findrotamername("side_mix_rotamer");
	if(s2==0) s2=TRES.findrotamername("side_large_rotamer");
	if(n) {
		if(s1->token) delete [] s1->token;
		s1->token=strdup("sidechain");
		if(s2->token) delete [] s2->token;
		s2->token=strdup("sidechain0");
	}
	else {
		if(s1->token) delete [] s1->token;
		s1->token=strdup("sidechain0");
		if(s2->token) delete [] s2->token;
		s2->token=strdup("sidechain");		
	}
}

void Mutate::minsidechain0() {

	setcontact(mpdb->chn->res,100000);
	mpdb->chn->setresenergy(mpdb->chn->res,100000,4);
		
	int nres=mpdb->chn->lastres()->id0+100;
	
	Res **ren=new Res*[nres];

	float *temp=new float[nres];
	int   *order=new int[nres];

	int n=0;Res *r;
	for(r=mpdb->chn->res;r;r=r->next) {
		if(strchr("PGA",r->name)) continue;
		ren[n]=r;
		temp[n]=-r->energy;
		n++;
	}

	Qsort cc;

	cc.sort(temp,n,order);	

	int m;
	for(m=0;m<n;m++) {
		if(-temp[m]<0) break;
	}	
	
	if(m==0) {
		delete [] temp;temp=0;
		delete [] order;order=0;
		delete [] ren;ren=0;
		return;
	}

	int i;
	mpdb->setflgr(-99999);
	for(i=0;i<m/2;i++) {		
		int j=order[i];
		r=ren[j];
		r->flag=10000;
	}
	Scap scprd; 
	//08/16/2002
	strcpy(scprd.force,"124");
	TRES.setdonaldcharge();
	scprd.dielectric=15;
	//end
	scprd.singletorsion=1;
	scprd.colonyline=1;
	scprd.colony=2;
	scprd.ncolony=1;
	scprd.nncolony=1;
	scprd.nummore=100;
	scprd.bmax=2;
	scprd.tormax=2;
	scprd.ring=1;
	scprd.pdb=mpdb;
	scprd.includeself=1;		
	scprd.resultout=0;
	scprd.scpred(1);
	scprd.pdb=0;
	for(r=mpdb->chn->res;r;r=r->next) r->setfreenextonmore();
 	TRES.switchcharge(1);
	delete [] temp;temp=0;
	delete [] order;order=0;
	delete [] ren;ren=0;
	if(TRES.logg)mpdb->write("side");
}

int Mutate::getlength(int mm) {

	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 1;

	//only insertion case
	if(se->seqngap[mm]=='-'||sqnto[mm]!='-') return 1;

	int n;

	int i;

	n=0;
	for(i=mm;i>=0;i--) {
		if(se->seqngap[i]!='-'&&owner->seqngap[i]=='-') n++;
		else break;
	}
	
	int nlen=strlen(sqnto);
	for(i=mm+1;i<nlen;i--) {
		if(se->seqngap[i]!='-'&&owner->seqngap[i]=='-') n++;
		else break;
	}
	
	if(n<10) return 100;
	else return 100;
}

void Mutate::minsegment(Res *r1,Res *r2){
	
	setcontact(r1,r2->id0);
	
	int n1,n2;
	
	Atm *a;
	while(1) {
 
		int m1=compare[r1->id0];
		int m2=compare[r2->id0];
 
		n1=0;
		n2=0;
			
		int i;	
		for(i=m1;i<=m2;i++) {
			if(sqnto[i]!='-') n1++;
			if(owner->seqngap[i]!='-') n2++;
		}	
 
		if(n1<n2) break;
		
		a=findclash(r1,r2);
		
		if(a->contact->total<0) break; break;				
	}
}




 
Atm *Mutate::findclash(Res *r1,Res *r2) {

	Res *r;
	Atm *a,*b;
	float d=-999999;
	b=0;
	for(r=r1;r;r=r->next) {
		if(r->id0>r2->id0) break;
		for(a=r->atm;a;a=a->next) {
			if(a->contact==0) continue;
			if(a->contact->total>d) {
				d=a->contact->total;
				b=a;
			}
		}	
	}
	return b;
}
 
void Mutate::setcontact(Res *r1,int nid) {

	mpdb->deletecontact();
	//
	Lattice *lat=new Lattice();
	lat->putoff();
  	lat->flag=1;
  	lat->grdsiz=2.;
  	lat->radall=15.;
  	lat->ready(mpdb);
	lat->puton(mpdb);
	
	Res *t;
	Atm *a,*all[50],*b;
	int m,j,n,i,k;
 
  	for(t=r1;t;t=t->next) {
		if(t->id0>nid) break;
		t->energy=0;
		for(a=t->atm;a;a=a->next) {
			
			j=0;n=0;
			a->energy=0;
  			while(a->near[j]){
    				if(lat->exist(a->near[j])) { 
       					all[n++]=a->near[j];
       					lat->putoff(a->near[j]);
    				}
    				j++; 
  			}
			
			lat->getcell(a,6);
			lat->cutoff(a,6);
			a->contact=new Contact();
			m=lat->nget;
			a->contact->setatmarray(m);
			a->contact->setengarray(m);
			a->contact->total=0;
			k=0;
			for(m=0;m<lat->nget;m++) {
				b=lat->obtain[m]->atm;						 
				a->contact->atm[k]=b;	
				a->contact->energy[k]=clash(a,b);
				a->energy+=a->contact->energy[k];				
				k++;
			}
			a->contact->num=k;
			a->contact->total=a->energy;
			t->energy+=a->contact->total;
			if(0) {
				//for(i=0;i<n;i++)cerr<<i<<" "<<a->name<<" "<<a->res->id0<<" "<<all[i]->name<<" "<<all[i]->res->id0<<endl;	
			}
			for(i=0;i<n;i++)lat->puton(all[i]);
		}
		if(TRES.logg>3) cerr<<"energy: "<<t->name<<t->id0<<" "<<t->energy<<endl;
	}

	delete lat;lat=0;
}



float Mutate::clash(Atm *aa0,Atm *aa1) {

	float ta,tb,tc,tn;
	float dr,dn,e,d,elon,xr[3],ds,tt,de;
  	int   j,isp;  
  	int   mr; 

  	if(segen->smoothclash) {
        	ta=TRES.engcoeff->coeff[0];
        	tb=TRES.engcoeff->coeff[1];
        	tc=TRES.engcoeff->coeff[2];
        	tt=TRES.engcoeff->coeff[3];
        	tn=TRES.engcoeff->coeff[4];
  	}
  	else {
        	ta=TRES.smooth[TRES.smt*4+0];
        	tb=TRES.smooth[TRES.smt*4+1];
        	tc=TRES.smooth[TRES.smt*4+2];
        	tn=TRES.smooth[TRES.smt*4+3];
        	tt=0;
  	}
 
	mr=0; e=0;
     	d=aa0->tatm->eng->radius+aa1->tatm->eng->radius;
     	j=aa0->tatm->hbond*aa1->tatm->hbond;
     	if(j==2||j==3||j==6||j==9){d=3.0;mr=1;}
     	if(aa0->tatm->hbond==5) d-=aa0->tatm->eng->radius;
     	if(aa1->tatm->hbond==5) d-=aa1->tatm->eng->radius;
     	if(d<=0) return 0;
     	if(aa0->res->name=='C'&&aa1->res->name=='C')
     	{
        	if(aa0->tatm->id==5&&aa1->tatm->id==5) {d=2.0;mr=2;}
        	else if(aa0->tatm->id==5&&aa1->tatm->id==4) d=3.;
        	else if(aa0->tatm->id==4&&aa1->tatm->id==5) d=3.;
     	}
     	if(aa0->tatm->name[1]=='H'&&aa1->tatm->name[1]=='H') d=1.0;
     	for(j=0;j<3;j++) xr[j]=(aa0->xyz[j]-aa1->xyz[j])/d;
     	dr=TRES.distsqr(xr);
 	ds=dr;
     	if(dr<0.1) dr=0.1; 
     	else if(dr>4) return 0;
     	if(aa0->res->name=='C'&&aa1->res->name=='C')
     	{
        	if(aa0->tatm->id==5&&aa1->tatm->id==5)e+=-10.*exp(-(dr-1)*(dr-1)); 
     	}
     	elon=sqrt(aa0->tatm->eng->epslon*aa1->tatm->eng->epslon);
     	elon=ta*elon*exp(-dr*tb);
     	dr=1/(tt+dr);
      		
     	dn=pow(dr,tn);  
    	dn=elon*dn*(dn-tc);   
     	 
     	isp=aa0->tatm->ispolar*aa1->tatm->ispolar;
     	if(dn<0)
     	{
      		if(mr==1)       e+= dn*1.50;
      		else if(mr==2)  e+= dn*5;
      		else if(isp==1) e+= dn*1.50;  
      		else if(isp==2) e+= dn*1.25; 
      		else if(isp==3) e+= dn*0.75; 
      		else            e+= dn;
     	}
     	else
     	{
      		if(mr==1)       e+= dn*0.50;
      		else if(mr==2)  e+= dn*0.20;
      		else if(isp==1) e+= dn*0.50; 
      		else if(isp==2) e+= dn*0.75;
      		else if(isp==3) e+= dn*1.25;
      		else            e+= dn;
     	}

	if(aa0->tatm->eng->charge!=0&&aa1->tatm->eng->charge!=0)
     	{
            de=1/(0.5+ds*d*d);
	    //de=sqrt(de);
            dn=aa0->tatm->eng->charge*aa1->tatm->eng->charge*de*332.3/20;
	    //if(dn>0) dn=0;
            e+=dn;
     	}
	return e;
}



void Mutate::setreson(Res **tmp) {
	Res **reson=0;

	if(reson==0) {
		int n=strlen(sqnto);
		reson=new Res*[n+10];
		reson[0]=0;
	}
	if(tmp==0||tmp[0]==0) reson[0]=0;
	else {
		int n=0;
		while(tmp[n]) {
			reson[n]=tmp[n];
			n++;
		}
		reson[n]=0;
	}
}

void Mutate::setnearreliable() {

	int nlen=strlen(sqnto);
 
	if(rely) delete [] rely;
	rely=new int[nlen];
	
	int i,j;
	for(i=0;i<nlen;i++)  rely[i]=1000;
 
	Res *r,*r0;
	for(r=mpdb->chn->res;r;r=r->next) {
		i=compare[r->id0];
		j=owner->match[i];
		r0=owner->resn[j];
		if(r0==0) continue;
		rely[r->id0]=(int)r->directrmsdanyway(r0,0,3);		
	}
}

Res *Mutate::findhelixend(Res *s,int f) {

	int n;
	if(f==0) {
		Res *t0=0;
		for(Res *t=s;t;t=t->last){
			n=compare[t->id0];
			if(dssp[n]!='h') break;
			t0=t;	
		}
		return t0;
	}
	else {
		Res *t0=0;
		for(Res *t=s;t;t=t->next){
			n=compare[t->id0];
			if(dssp[n]!='h') break;
			t0=t;	
		}
		return t0;
	}
}

int * Mutate::actallhelix(int *stem) {

	int n1=stem[0];
	int n2=stem[1];
 
	//all helix?
	int i;
	for(i=n1;i<=n2;i++) {
		if(dssp[i]!='h') return stem;
	}
 	
	Res *r1=resn[match[n1]];
	Res *r2=resn[match[n2]];
		 
	StrFmt *parent=owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return 0;

	if(r1->last==0||r2->next==0) return stem;

	Res *rr;
	 
	Res *m1=findhelixend(r1,0);
	Res *m2=findhelixend(r2,1);
		
	if(m1==0||m1->last==0) {
		int *out=new int[4];
		out[0]=compare[m1->id0];
		out[1]=stem[1];
		out[2]=stem[2];
		out[3]=0;
		delete [] stem;stem=0;
		return out;
	}
	else if(m2==0||m2->next==0){
		int *out=new int[4];
		out[0]=stem[0];
		out[1]=compare[m2->id0];		
		out[2]=stem[2];
		out[3]=0;
		delete [] stem;stem=0;
		return out;
	}

	float a1=0,a2=0;
	
	for(rr=r1;rr;rr=rr->last) {
		if(rr->id0>=m1->id0) break;
		i=compare[rr->id0];
		a1+=parent->sitescore[i];			 		
	} 

	for(rr=r2;rr;rr=rr->next) {
		if(rr->id0>=m2->id0) break;
		i=compare[rr->id0];
		a2+=parent->sitescore[i];			 
	} 
	
	if(a1>a2) {
		int *out=new int[4];
		out[0]=stem[0];
		out[1]=compare[m2->id0];		
		out[2]=stem[2];
		out[3]=0;
		delete [] stem;stem=0;
		return out; 
	}
	else {
		int *out=new int[4];
		out[0]=stem[0];
		out[1]=compare[m2->id0];		
		out[2]=stem[2];
		out[3]=0;
		delete [] stem;stem=0;
		return out; 
	}
}

int *Mutate::resetstemcomposite(int *stem) {
	
	int n1,n2;

	n1=stem[1];n2=stem[0];

	int i;
	for(i=stem[0];i<=stem[1];i++) {
		if(sqnto[i]=='-') continue;
		n1=min(n1,i);
		n2=max(n2,i);
	}
	stem[0]=n1;
	stem[1]=n2;
	return stem;
} 
