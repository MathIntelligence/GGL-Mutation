#include "source.h"

Chiangle::Chiangle() {
	rotdir=1;
	name=0;
	angle=0;
	numa=0;
	next=0;
	number=0;
	pdb=0;
	filename=0;
	omega=0;
	chain=0;
	oid=0;
	endres=100000;
	wavedir=1;
}

Chiangle::~Chiangle(){

	Strhandler cc;
	name=cc.strdel(name);
	numa=cc.intdel(numa);
	angle=cc.floatdel(angle);
	chain=cc.strdel(chain);
	oid=cc.intdel(oid);
	if(next) delete next;
	if(pdb) delete pdb;
	if(filename) delete [] filename;filename=0;
}

void Chiangle::setuprecord(int tid,float ang) {

	
	Atm *a=pdb->getatmbyoid(tid);
	
	if(a==0) {
		cerr<<"no atom with id:"<<tid<<endl;
		exit(0);
	}
	
	name=new char[10];
	name[0]=a->res->name;
 	numa=new int[10];
	angle=new float*[10];
	int i;
	for(i=0;i<10;i++) angle[i]=0;
	angle[0]=new float[10];	
	oid=new int[10];
	oid[0]=a->res->oid;
	chain=new char[10];
	chain[0]=a->res->chn->id;
	Atm *b;int m=0;int gt=0;
	for(b=a->res->atm;b;b=b->next) {
 
		if(b->tatm->id==0) {
			m++;
		}
		else if(b->tatm->id==1) {				 
			m++;
		}
		else if(b->tatm->id==2) {				 
			m++;
		}
		else if(b->tatm->rotate) {				 
			m++;
		}
		else {
			continue;
		} 
		if(b==a) {gt=1;angle[0][m-1]=ang;}
		else angle[0][m-1]=999;	
	}
	numa[0]=m;
	number=1;
	Atm *a0=a->bond[0];
	if(gt==0&&a0) {
		if(a->res->chn->id!=' '&&a->res->chn->id!='-') cerr<<"the bond in chain: "<<a->res->chn->id<<endl;
		cerr<<"the bond is: "<<a0->res->name<<a0->res->oid<<" "<<a0->name<<" -- "<<a->res->name<<a->res->oid<<" "<<a->name<<endl;
		cerr<<"the above bond can not be rotated"<<endl;		 
		number=0;
		return;	 
	}
	else if(a0) {

		if(a->res->chn->id!=' '&&a->res->chn->id!='-') cerr<<"the bond in chain: "<<a->res->chn->id<<endl;
		cerr<<"the bond is: "<<a0->res->name<<a0->res->oid<<" "<<a0->name<<" -- "<<a->res->name<<a->res->oid<<" "<<a->name<<endl;
		 
	}
	else if(a0==0) { 
		cerr<<"the bond at: "<<a->res->name<<a->res->oid<<" "<<a->name<<endl;
		cerr<<"does not have parent atoms"<<endl;
		cerr<<"ignored.."<<endl;
		number=0;
		return;
	}
}
void Chiangle::read(char *filename) {

	FILE *fp;
	int i;
	char *mlib,line[200];

	fp=0;
	fp=fopen(filename,"r");
	
	mlib=RCS["library"];
	i=0;	
	sprintf(line,"%s/%s",mlib,filename);
	i=strlen(mlib);
	while(fp==0&&(fp=fopen(line,"r"))==NULL)
  	{
   		mlib=mlib+i+1;
   		i=strlen(mlib);
   		if(*mlib=='\0')
   		{
			cerr<<"warning:could not open: "<<filename<<endl;
			return;
		}
   		sprintf(line,"%s/%s",mlib,filename);
  	}

	Strhandler cc;
	int m=strlen("#start");

	Chiangle *tt=0;		

	char myname[2000];
	int  mynuma[2000];
	float *myangle[2000];
	int  mynumber;
	while(fgets(line,200,fp)!=NULL) {

		char *s=strdup(line);

		s=cc.clearendchar(s,"\n\t ");

		if(s==0) continue;
		if(strlen(s)<=1) {
			delete [] s;s=0;continue;
		}
		if(s[0]=='!') {
			delete [] s;s=0;
			continue;
		}
		if(strncmp(s,"#start",m)==0) {
			if(tt==0)  tt=this;
			else	{
				tt->name=cc.opnstr(myname);
				tt->numa=cc.opnint(mynuma,mynumber);
				tt->angle=new float*[mynumber];
				for(int ii=0;ii<mynumber;ii++) tt->angle[ii]=myangle[ii];
				tt->number=mynumber;
				tt->next=new Chiangle();
				tt=tt->next;	
			}  
			myname[0]='\0';     
			mynuma[0]=0;
			for(int ii=0;ii<2000;ii++) myangle[0]=0;
			mynumber=0;
		}
		if(s[0]=='#') {
			delete [] s;
			s=0;
			continue;
		}
		
		myname[mynumber]=s[0];
		myangle[mynumber]=parseangle(s);
		mynuma[mynumber]=parseanglenum(s);
		mynumber++;
	}
	fclose(fp);
	if(mynumber>0) {
		tt->name=cc.opnstr(myname);
		tt->numa=cc.opnint(mynuma,mynumber);
		tt->angle=new float*[mynumber];
		for(int ii=0;ii<mynumber;ii++) tt->angle[ii]=myangle[ii];
		tt->number=mynumber;
	}

	mynumber=0;
	for(tt=this;tt;tt=tt->next) mynumber++;
	cerr<<"the total number of conformations of protein segment read:"<<mynumber<<endl;
	
}

void Chiangle::readone(char *filename0) {

	FILE *fp;
 
	char line[200];

	if(filename0==0) {
		cerr<<"the file name is not given"<<endl;
		exit(0);
	}
	filename=strdup(filename0);
	fp=fopen(filename,"r");
	if(fp==0) {
		cerr<<"the file:"<<filename<<" does not exist"<<endl;
		exit(0);
	}
	
	Strhandler cc;
	int m=strlen("#start");

	Chiangle *tt=0;		

	int myoid[2000];
	char mychain[2000];
	char myname[2000];
	int  mynuma[2000];
	float *myangle[2000];
	int  mynumber;
	while(fgets(line,200,fp)!=NULL) {

		char *s=strdup(line);

		s=cc.clearendchar(s,"\n\t ");

		if(s==0) continue;
		if(strlen(s)<=1||s[0]=='!') {
			delete [] s;s=0;continue;
		}
		if(strncmp(s,"#start",m)==0) {
			if(tt==0)  tt=this;
			else	{
				mychain[mynumber]='\0';
				tt->chain=cc.opnstr(mychain);
				tt->oid=cc.opnint(myoid,mynumber);
				myname[mynumber]='\0';
				tt->name=cc.opnstr(myname);
				tt->numa=cc.opnint(mynuma,mynumber);
				tt->angle=new float*[2*mynumber];
				int ii;
				for(ii=0;ii<2*mynumber;ii++) tt->angle[ii]=0;
				for(ii=0;ii<mynumber;ii++) tt->angle[ii]=myangle[ii];
				tt->number=mynumber;
				tt->next=new Chiangle();
				tt=tt->next;	
			}  
			myname[0]='\0';     
			mynuma[0]=0;
			for(int ii=0;ii<2000;ii++) myangle[0]=0;
			mynumber=0;
		}
		if(s[0]=='#') {
			delete [] s; s=0;
			continue;
		}
		Tres *tr=TRES[s[0]];
		if(tr==0) {
			cerr<<"Warning! ignored: "<<s<<endl;
			delete [] s;s=0;
			continue;
		}
		myname[mynumber]=s[0];
		myoid[mynumber]=atoi(s+1);
		char tmp[100];
		sprintf(tmp,"%i",myoid[mynumber]); 
		char *ss=strstr(s,tmp);
		if(ss==0) {
			delete [] s;s=0;
			continue;
		}
		ss+=strlen(tmp);
		if(ss==0) {
			delete [] s;s=0;
			continue;
		}
		int nss=strlen(ss);
		int ii;
		char jj='-';
		for(ii=0;ii<nss;ii++) {
			if(ss[ii]==' ') continue;			
			if((ss[ii]>='A'&&ss[ii]<='Z')||(ss[ii]>='a'&&ss[ii]<='z')) {
				jj=ss[ii];
				ss[ii]=' ';
			} 
			else break;
		} 
		mychain[mynumber]=jj;
		myangle[mynumber]=parseangle(s);
		mynuma[mynumber]=parseanglenum(s);		
		
		mynumber++;
		if(s) delete [] s; s=0;
	}
	fclose(fp);
	if(tt==0) {
		tt=this;
	}
	if(mynumber>0) {
		mychain[mynumber]='\0';
		tt->chain=cc.opnstr(mychain);
		tt->oid=cc.opnint(myoid,mynumber);
		myname[mynumber]='\0';
		tt->name=cc.opnstr(myname);
		tt->numa=cc.opnint(mynuma,mynumber);
		tt->angle=new float*[mynumber*2];
		int ii;
		for(ii=0;ii<mynumber*2;ii++) tt->angle[ii]=0;
		for(ii=0;ii<mynumber;ii++) tt->angle[ii]=myangle[ii];
		tt->number=mynumber;
	}
	for(tt=this;tt;tt=tt->next) {
		cerr<<"the total number of residues specified:"<<tt->number<<endl;
	}
}


void Chiangle::buildallstructure() {

	Chiangle *tt;
	for(tt=this;tt;tt=tt->next) tt->buildstructure();
}
 
void Chiangle::buildstructure() {

	Rotate rot;

	if(name==0) return;
	int n=strlen(name);
	if(n==0) return;
	
	//
	if(!pdb) {
		char *name0=new char[n+100];
		sprintf(name0,"A%sA",name);
		
		pdb =new Pdb;
		pdb->chn=new Chn;
		pdb->chn->create(name0);
		pdb->configure();
		pdb->chn->dihedral();
		if(name0) delete [] name0;name0=0;
		Res *r;
		int n=0;
		for(r=pdb->chn->res->next;r;r=r->next) { 
			if(r->next==0) continue;
			r->oid=oid[n];n++;	
		}
		pdb->chn->res->oid=-100000;
		pdb->chn->lastres()->oid=100000;
	}
	else {
		Chn *c;
		for(c=pdb->chn;c;c=c->next) {
			c->header();
			Tres *t=TRES['A'];
			Res *a=new Res(t); 
			 
			a->nemp=9;
			a->temp=new float[9];
			int ii;
			for(ii=0;ii<9;ii++) a->temp[ii]=t->head[ii];
			 
			a->chn=c;
			a->id=0;
			c->increaseid(c->res,100000,1);			
			a->next=c->res;
			c->res=a;			 			
			c->configure();
			rot.link(a,a->next,0,0);
		 	if(TRES.logg>3) c->write("ss1");	
			a=new Res(t); 
			 
			a->nemp=9;
			a->temp=new float[9];
		 
			for(ii=0;ii<9;ii++) a->temp[ii]=t->head[ii];
			 
			a->chn=c;
			Res *r=c->lastres();
			r->next=a;a->id=r->id+1;
			c->configure();
			rot.link(r,a,a->id0,1);	
			c->res->oid=-100000;
			c->lastres()->oid=100000;
			if(TRES.logg>3) c->write("ss2");
			c->header();	
		}
	}

	pdb->setlastresorder();

	int endres0=endres;

	n=0;
	Res *r;
	Atm *a;
	Chn *c;
	pdb->setflgr(-99999);
	int ii=0;
	for(ii=0;ii<number;ii++) {
		if(chain&&chain[ii]=='-') c=pdb->chn;
		else c=pdb->ischain(chain[ii]); 
		if(c==0) {
			cerr<<"no chain with id:"<<chain[ii]<<endl;
			cerr<<"ignored.."<<endl;
			continue;
		}
		endres=endres0;
		Res *eres0=0;
		if(c->lastres()->last==0||c->res->next==0) {
			cerr<<"no residue in the pdb just read"<<endl;
			continue;
		}
		if(endres>c->lastres()->last->oid) eres0=c->lastres();
		else if(endres<c->res->next->oid)   eres0=c->res;	
		else eres0=c->isresoid(endres);
		endres=eres0->oid;

		r=c->isresoid(oid[ii]);
		if(r==0) {
			if(chain[ii]=='-')  
			cerr<<"Warning! no residue found for residue with id:"<<oid[ii]<<endl;
			else 
			cerr<<"Warning! no residue found for residue with id:"<<oid[ii]<<" chain id:"<<c->id<<endl;
			continue;
		}
		if(r->name!=name[ii]) {
			cerr<<"Warning! residue name does not match for res:"<<name[ii]<<" "<<r->name<<endl;
			continue;
		}
	  	if(wavedir==1&&r->oid>=endres) {
			cerr<<"warning! residue id "<<r->name<<r->oid<<" trespasses the end id of the segment: "<<endres;
			cerr<<"-->Ignored!"<<endl;
			continue;
		}
		else if(wavedir==0&&r->oid<=endres) {
			cerr<<"warning! residue id "<<r->name<<r->oid<<" trespasses the end id of the segment: "<<endres;
			cerr<<"-->Ignored!"<<endl;
			continue;
		}
		int m=0;
		for(a=r->atm;a;a=a->next) {
			
			float *f=angle[ii];
			 
			if(f==0) continue;
			float ang=a->chi;
			if(a->tatm->id==0&&omega==0) {
				m++;
			}
			else if(a->tatm->id==1&&(omega==0||omega==1||omega==2)) {				 
				m++;
			}
			else if(a->tatm->id==2&&(omega==0||omega==1||omega==2)) {				 
				m++;
			}
			else if(a->tatm->rotate&&a->tatm->id>=4&&(omega==0||omega==1||omega==3)) {				 
				m++;
			}
			else {
				continue;
			}
			if(m>numa[ii]) continue;
			ang=f[m-1];
			if(ang>360) continue;
			float rt=a->gettorsionangle();
			float ee=ang-rt;
			if(wavedir==0) ee=-ee; //modified to take account of reverse rotation
			if(rotdir==0&&wavedir==1) ee=ang;
			else if(rotdir==0&&wavedir==0) ee=-ang;//modified
			if(ee==0) continue;
			
			if(a->tatm->id==1||a->tatm->id==2) {
				//rot.rotate(a,ee,1);								
				//cerr<<a->res->id0<<" "<<eres0->id0<<wavedir<<endl;
				rot.rotate(a,ee,eres0->id0,wavedir);
			}
			
			else if(a->tatm->id==0&&wavedir==1) {
				//rot.rotateomega(a,ee,100000,1);
				//cerr<<a->res->id0<<" "<<eres0->id0<<" "<<wavedir<<endl;
				rot.rotateomega(a,ee,eres0->id0,wavedir);
			}
			
			else {
				rot.rotate(a,ee);
			}
		}

		if(r->name=='P') {
			pdb->setflgr('P',10000);
		}
	}

 	endres=endres0;

	pdb->transform(5,-99999);
	Scap scprd;
	scprd.pdb=pdb;
	scprd.hookside();
	scprd.pdb=0;
		

	for(c=pdb->chn;c;c=c->next) {
	 
		r=c->res;
		c->res=r->next;
		r->next=0;
		delete r;
		r=c->lastres();
		Res *r0=c->isres0(r->id0-1);
		r0->next=0;
		r->next=0;
		delete r;
		pdb->configure();
	}
}


int Chiangle::parseanglenum(char *line) {
	Strhandler cc;
        char *s=strdup(line);
        s=cleardoublespace(s);
        s=cc.clearendchar(s,"\n\t ");
        if(s==0) return 0;

        char **tt=cc.pairbytoken(s," ");

        int n=cc.gettotnum(tt);
	if(n<=2) {
                delete [] s;s=0;
                return 0;
        }
	else {
		tt=cc.strdel(tt);
		s=cc.strdel(s);
		return n-2;
	}
}

float *Chiangle::parseangle(char *line) {

	Strhandler cc;
	char *s=strdup(line);
	s=cleardoublespace(s);
	s=cc.clearendchar(s,"\n\t ");
	if(s==0) return 0; 

	char **tt=cc.pairbytoken(s," ");

	int n=cc.gettotnum(tt);

	if(n<=2) {
		delete [] s;s=0;
		return 0;
	}
	float *d=new float[n-2];
	int j=0;
	for(int i=2;i<n;i++) {
		d[j++]=atof(tt[i]);
	}
	tt=cc.strdel(tt);
	s=cc.strdel(s);
	return d;
}

char *Chiangle::cleardoublespace(char *s) {


	while(1) {
		int n=strlen(s);
		int i=0,j=0,m=0;
		for(i=0;i<n;i++) {
			if(i==n-1) {s[j++]=s[i];continue;}
			if(s[i]==' '&&s[i+1]==' ') {m++;continue;}
			s[j++]=s[i];			
		}
		s[j]='\0';
		if(m==0) return s;
	}
}

int Chiangle::ifexist(Chiangle *s){

	if(s==0) return 1;

	Chiangle *t;

	for(t=this;t;t=t->next) {
		if(t->number!=s->number) continue;
		int i; int nn=1;
		for(i=0;i<t->number;i++) {
			if(t->numa[i]!=s->numa[i]) {nn=0;break;}
			int j=t->numa[i];
			int k;
			for(k=0;k<j;k++) {
				if(t->angle[i][k]!=s->angle[i][k]) {nn=0;break;}
			}
			if(nn==0) break;
		}
		if(nn==1) return 1;
	}
	return 0;
}

void Chiangle::add(Chiangle *s) {

	if(s==0) return;
	Chiangle *t;
	if(ifexist(s)==1) {
		delete s;s=0;
		return;
	}
	for(t=this;t->next;t=t->next);
	t->next=s;
}

void Chiangle::addanyway(Chiangle *s) {

        if(s==0) return;
        Chiangle *t;
        for(t=this;t->next;t=t->next);
        t->next=s;
}

void Chiangle::size(int n){
 
	Strhandler cc;
	angle=cc.floatdel(angle);
	numa=cc.intdel(numa);
	
	angle=new float*[n];
	numa=new int[n];
	int i;
	for(i=0;i<n;i++) {
		angle[i]=0;
		numa[i]=0;
	}
	 
}


void Chiangle::resize(int n){
 
	angle=(float **)realloc(angle,n*sizeof(float));
	numa=(int *)realloc(numa,n*sizeof(int));	
	 
}
 
Chiangle *Chiangle::get(int n) {

	int nn=0;

	Chiangle *s;

	for(s=this;s;s=s->next) {

		if(nn==n) return s;
		nn++;
	}
	return 0;
}

int Chiangle::getsize(){

	int nn=0;

	Chiangle *s;

	for(s=this;s;s=s->next) {

		nn++;

	}
	return nn;
}

float Chiangle::getangle(int nn,int i) {

	if(angle==0||numa==0) return 999;

	if(nn>=number) return 999;
	
	if(i>=numa[nn]) return 999;
	
	if(angle[nn]==0) return 999;

	return angle[nn][i];
}	
