#include "source.h"

StrFmt::StrFmt() {
       initial();
}
void StrFmt::initial() {
	profix=0;
	alnback=1;
        tune=0;
	fileName=0;
	code=0;
	start=1;
	end=1;
	cid='0';
	seqn=0;
	seqngap=0;
    	match=0;     
	compare=0;
	source=' ';
	pdb=0;
	next=0;
	more=0;
	up=0;
	last=0;
	token=0;
	flag=0;
	resn=0;
	mutate=0;
	sitescore=0;
	oldgap=0;
	matrix=strdup("blosum65mt");
	gap=findmatrixminvalue()*4;
	zipcode='\0';
	addh=0;
	other=0;
	minfast=0;
	onlyrefine=0;
	//shiftcut=0.0;
	iscom=0;
	numaln=0;
	refstart=-100000;
	refend=100000;
}
StrFmt::~StrFmt() {

	Strhandler cc;
	fileName=cc.strdel(fileName);
	code=cc.strdel(code);
	seqn=cc.strdel(seqn);
	seqngap=cc.strdel(seqngap);
    	if(match) delete [] match;
    	if(mutate) delete mutate;
	if(compare) delete [] compare;
	if(pdb&&(flag&TRES.constant->pdbundeletable)==0) delete pdb;
	if(next) delete next;
	if(more) delete more;
	if(sitescore) delete [] sitescore;
	if(oldgap) delete [] oldgap;
	if(matrix) delete [] matrix;
	last=0;up=0;
}

StrFmt::StrFmt(StrFmt *se) {
	initial();
	tune=se->tune;
	fileName=0;
	if(se->fileName) fileName=strdup(se->fileName);
	code=0;
	if(se->code) code=strdup(se->code);
	start=se->start;
	end=se->end;
	cid=se->cid;
	seqn=0;
	if(se->seqn) seqn=strdup(se->seqn);
	seqngap=0;
	if(se->seqngap) seqngap=strdup(se->seqngap);
    	match=0;     
	compare=0;
	if(se->match&&se->compare) {
		int nlen=strlen(seqngap);
		match=new int[nlen];
		compare=new int[nlen];
		int i;
		for(i=0;i<nlen;i++) {
			match[i]=se->match[i];
			compare[i]=se->compare[i];
		}
	}
	source=se->source;
	pdb=0;
	if(se->pdb) {
		pdb=new Pdb(se->pdb);
		pdb->configure();
	}
	next=0;
	more=0;
	up=0;
	last=0;
	token=0;
	if(se->token) token=strdup(se->token);
	flag=se->flag;
	resn=0;
	if(se->resn) {
		int nlen=strlen(seqngap);
		resn=new Res*[nlen];
		int i;
		for(i=0;i<nlen;i++) {
			resn[i]=0;
		}

		int n=strlen(seqn);

		for(i=0;i<n;i++) {
			if(se->resn[i]==0) continue;
			resn[i]=pdb->chn->isres(se->resn[i]->id);
			if(resn[i]==0) {
				cerr<<se->resn[i]->name<<se->resn[i]->id<<"  has error in composite mapping!"<<endl;
				continue;
			}
		}
	}
	mutate=0;
	sitescore=0;
	if(se->sitescore) {
		
		int nlen=strlen(seqngap);
		sitescore=new float[nlen];
		int i;
		for(i=0;i<nlen;i++) sitescore[i]=se->sitescore[i];
	}
	oldgap=0;
	if(se->oldgap) oldgap=strdup(se->oldgap);
	matrix=0;
	if(se->matrix) matrix=strdup(se->matrix);
	gap=se->gap;
	zipcode=se->zipcode;
	addh=se->addh;
	other=0;
	//shiftcut=se->shiftcut;
	iscom=0;
	numaln=se->numaln;
	int i;
	for(i=0;i<numaln;i++) {
		resaln[i]=se->resaln[i];
	}
	refstart=se->refstart;
	refend=se->refend;
	profix=se->profix;
}

void StrFmt::setgapvalue(){

}

void StrFmt::setmatrixname(char *s) {

	StrFmt *t,*t0;

	for(t=this;t;t=t->next)
        for(t0=t;t0;t0=t0->more) {
		if(t0->matrix) delete [] t0->matrix;
		t0->matrix=strdup(s);
	}
}

void StrFmt::setmutateresn() {
	
	Pdb *mutatepdb=0;

	StrFmt *t,*t0;

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

int  StrFmt::isbreaker(){

	Res *r;
	pdb->chn->header();
	//check id sequence
	for(r=pdb->chn->res;r;r=r->next) {
		if(r->next==0) continue;
		if(r->next->id-r->id>1&&pdb->chn->islinked(r,r->next)==1) {	
			cerr<<"warning!"<<endl;
			cerr<<"the residue "<<r->next->name<<r->next->oid<<" "<<r->name<<r->oid<<" is linked!"<<endl;
			cerr<<"they should have sequential id numbers"<<endl;
			int n=r->next->id-r->id-1;
			Res *t;
			for(t=r->next;t;t=t->next) t->id-=n;	
		}		 
	}

	for(r=pdb->chn->res;r;r=r->next) {
		if(r->next==0) continue;
		if(r->next->id-r->id>1)  return r->id;
	}
	return -1;
}

void StrFmt::setresn() {

	if(seqn==0) return;
	StrFmt *root=getrootStrFmt();
	if(root->onlyrefine==0) {
		cerr<<endl;
		cerr<<"find corresponding residues between pdb and pir for the structure pir:"<<code<<endl;
		cerr<<endl;
	}
	if(resn==0) delete [] resn;resn=0;

	int  n=strlen(seqngap);
	int  nlen=n;

	resn=new Res *[n+100];
	
	int i;

	for(i=0;i<n+100;i++)   resn[i]=0;

	int error=0;
	//if the start and end position specified
	if(start!=-1&&end!=-1) {
		 				 
		Res *rr;

		//check if the start and end position right
        	int nn=0;
		int mm=strlen(seqn);
		for(i=start;i<=end;i++) {
			rr=pdb->chn->isresoid(i);
			if(rr==0) {
				continue;
			}
			else if(rr==0) {
				cerr<<endl<<endl;
				cerr<<"**********warning!***********"<<endl;
				cerr<<"error in specifying the start and end residue ids in the pir alignment:"<<endl;
				printoutonly(stderr); 
				cerr<<endl;
				cerr<<"the program fixs it anyway by doing alignments between residues in pdb and in pir...."<<endl;
				cerr<<endl<<endl;
				error=1;
				break;
			}
			else if(nn>=mm||rr->name!=seqn[nn]) {
				cerr<<endl<<endl;
				cerr<<"**********warning!***********"<<endl;
				cerr<<"error in specifying the start and end residue ids in the pir alignment:"<<endl;
				printoutonly(stderr); 
				cerr<<endl;
				cerr<<"the program fixs it anyway by doing alignments between residues in pdb and in pir...."<<endl;
				cerr<<endl<<endl;
				error=1;
				break;
			}
			nn++;	
		} 
         	if(nn!=mm&&error==0) {
			cerr<<endl<<endl;
			cerr<<"**********warning!***********"<<endl;
                        cerr<<"error in specifying the start and end residue ids in the pir alignment:"<<endl;
                        printoutonly(stderr);
                        cerr<<"the program fixs it anyway by doing alignments between residues in pdb and in pir...."<<endl;
                        cerr<<endl<<endl;
                        error=1;
		}
        	if(error==0) {
			nn=0;
			for(i=start;i<=end;i++) { 
                        	rr=pdb->chn->isresoid(i);
				if(rr==0) continue;
                        	else if(rr&&rr->name==seqn[nn]) {
					resn[nn]=rr;
					nn++;
				}
				else {
					error=1;
					break;
				}							 
                	}
			resn[nn]=0;
			seqn[nn]='\0';
		}
	}
	else {
		error=1;
	}
	Algn algn;
	if(error==0) {
		setmatch();
		return;
		//goto re200;
	}
	else {
		for(i=0;i<nlen;i++) resn[i]=0;
		setmatch();
	}
	cerr<<"performing Needleman-Wunsch alignment to find the corresponding residue"<<endl;
	cerr<<"in pdb file and in pir file..."<<endl;
	
	if(pdb->chn->seqcard==0) {
		pdb->chn->setseqcard();
		if(pdb->chn->seqcard==0) pdb->chn->seqcard=strdup(" "); 
	}

	//setmatch();
	
	algn.setsequence(pdb->chn->seqcard,seqn);
	algn.defalgnid();
	 
	//char ** result = cmain(pdb->chn->seqcard,seqn);
	
	int n1,n2;
	n1=0;
	n2=0;
	//int slen= strlen(result[0]);
	int slen=algn.getroutelength();
	char **result=algn.output(stderr);
	nlen=strlen(seqngap);
	for(i=0;i<nlen;i++) seqngap[i]='-';
	int ne=strlen(seqn);
	for(i=0;i<ne;i++) seqn[i]='-';
	for(i=0;i<slen;i++) {
		
		if(result[0][i]!='-'&&result[0][i]!='X') n1++;
		if(result[1][i]!='-'&&result[1][i]!='X') n2++;
		//identical
		if(result[0][i]==result[1][i]&&result[0][i]!='-'&&result[0][i]!='X') {
			resn[n2-1]=pdb->chn->isres0(n1-1);
			Res *rr=resn[n2-1];
			seqn[n2-1]=rr->name;
			int jj=compare[n2-1];
                        seqngap[jj]=rr->name;
		}
		//different name
		else if(result[0][i]!='-'&&result[0][i]!='X'&&result[1][i]!='-'&&result[1][i]!='X') {			
			Res *rr=pdb->chn->isres0(n1-1);
			resn[n2-1]=rr;			 
			seqn[n2-1]=rr->name;
			int jj=compare[n2-1];
			seqngap[jj]=rr->name;
			cerr<<endl;
			cerr<<"warning! the residue names in pdb and alignment are different:"<<rr->name<<rr->id0<<"--"<<result[1][i]<<endl;
			cerr<<"changing the residue name in alignment to that in pdb..."<<endl;
			cerr<<endl;
			
		}
		//delete case, not exist in pdb
		else if((result[0][i]=='-'||result[0][i]=='X')&&result[1][i]!='-'&&result[1][i]!='X') { 
			cerr<<endl;
			cerr<<"warning! residue: "<<result[1][i]<<" at position: "<<i<<" does not exist in PDB"<<endl;
			cerr<<"deleting this residue from the alignment..."<<endl;
			cerr<<endl;			 
		}	
		//insert case
		else if((result[1][i]=='-'||result[1][i]=='X')&&result[0][i]!='-'&&result[0][i]!='X') { 
			Res *rr=pdb->chn->isres0(n1-1);
			int j;
			int jt=0;		
			for(j=0;j<i;j++) {
				if(result[1][j]!='-'&&result[1][j]!='X') {
					jt++;
					break;
				}
			}

			for(j=i+1;j<slen;j++) {
				if(result[1][j]!='-'&&result[1][j]!='X') {
					jt++;
					break;
				}
			}
			if(jt==2) {
				cerr<<endl;
				cerr<<"warning! residue: "<<rr->name<<rr->id0<<" does not exist in pir alignment"<<endl;			
				cerr<<"trying to insert this residue to the alignment..."<<endl;	
				cerr<<endl;			
			}				 
		}				
	}
	
	seqn[n2]='\0';
	resn[n2]=0;

	//re200:

	int m=0;	
	for(i=0;i<n2;i++) {
		if(resn[i]==0) continue;
		resn[m]=resn[i];
		m++;
	}
	resn[m]=0;

	re200:
	setmatch();
	
	int kep=strlen(seqn);
	for(int ii=0;ii<kep;ii++) {
			if(TRES.logg>3) cerr<<ii<<" "<<resn[ii]->name<<" "<<seqn[ii]<<" "<<seqngap[compare[ii]]<<endl;
	}

	checkcont();
	
	Strhandler cc;
	cc.strdel(result);
	
	//treatbreaker();
}

void StrFmt::checkresn() {
	int n=strlen(seqn);
	int i;
	for(i=0;i<n;i++) {
		if(resn[i]==0) {
			int j=compare[i];
			seqngap[j]='-';
		}	
	}
}

void StrFmt::setsize(int nn) {
	
	if(seqngap==0||seqn==0) return;

	int ngap=strlen(seqngap);
	int mlen=strlen(seqn);

	char *s1;
	int i;

	if(seqngap) {
		s1=new char[nn];
		strcpy(s1,seqngap);
		delete []  seqngap; seqngap=s1;
	}

	if(seqn) {
		s1=new char[nn];
		strcpy(s1,seqn);
		delete []  seqn; seqn=s1;
	}

	Res **tt;
	if(resn) {
		tt=new Res*[nn];
		for(i=0;i<mlen;i++) tt[i]=resn[i];
		delete [] resn; resn=tt;
		resn[mlen]=0;
	}

	int *mt;
	if(match) {
		mt=new int[nn];
		for(i=0;i<ngap;i++) mt[i]=match[i];
		delete [] match;match=mt;
	}

	int *nt;
	if(compare) {
		nt=new int[nn];
		for(i=0;i<ngap;i++) nt[i]=compare[i];
		delete [] compare;compare=nt;
	}
	
}

void StrFmt::checkcont() {
//insert residues

	StrFmt *parent=getparentStrFmt();
	StrFmt *se=parent->findsequencefmt();	
	if(se==0) return;

	//set size
	int nt=strlen(seqngap);
	StrFmt *t;
	for(t=parent;t;t=t->more) t->setsize(nt+nt+1000);

 	//find need to be inserted case

	re200:
	int nfd=0;
	setmatch();
	int i;
	int nlen=strlen(seqn);
	for(i=0;i<nlen;i++) {
		if(i==nlen-1) continue;
		Res *r1=resn[i];
		Res *r2=resn[i+1];
		if(r1->next!=r2) {
			Res *tt=r1->next;
			int j=compare[i];
			if(seqngap[j+1]=='-') {
				seqngap[j+1]=tt->name;					
			} 
			else {
				int jj;	
				int j=compare[i];
				int nt=strlen(seqngap);
				
				for(jj=nt-1;jj>=j+1;jj--){
					for(t=parent;t;t=t->more) {
						t->seqngap[jj+1]=t->seqngap[jj];
					}
					//seqngap[jj+1]=seqngap[jj];
					//se->seqngap[jj+1]=se->seqngap[jj];
				}
				for(t=parent;t;t=t->more) {
					t->seqngap[nt+1]='\0';
				}
				//seqngap[nt+1]='\0';
				//se->seqngap[nt+1]='\0';
				seqngap[j+1]=tt->name;
				se->seqngap[j+1]='-';
				for(t=parent;t;t=t->more) {
					if(strcmp(t->token,"composite")!=0) continue;
					if(t->seqngap[j]==t->seqngap[j+2]) {
						t->seqngap[j+2]=t->seqngap[j]; 
					}   
					else if(t->seqngap[j]==se->zipcode) {
						t->seqngap[j+2]=t->seqngap[j]; 
					}
					else if(t->seqngap[j+1]==se->zipcode) {
						t->seqngap[j+2]=t->seqngap[j+1]; 
					}
					else {
						t->seqngap[j+2]=t->seqngap[j]; 
					}
				}				
			}
			int jj;	
			for(jj=nlen-1;jj>=i+1;jj--){
				resn[jj+1]=resn[jj];
			}
			resn[i+1]=tt;
			//for(jj=i+1;jj<nlen;jj++) {
			for(jj=nlen-1;jj>=i+1;jj--){
				seqn[jj+1]=seqn[jj];
			}
			seqn[i+1]=tt->name;
			seqn[nlen+1]='\0';
			nfd++;
			goto re200;
		}
	}

	if(nfd==0) return;

}

int StrFmt::findpost(Res *t) {

	int nlen=strlen(seqn);

	int i;

	for(i=0;i<nlen;i++) {
		if(t==resn[i]) {
			return compare[i];
		}
	}
	return -1;
}

void StrFmt::treatbreaker0() {

	StrFmt *parent=getparentStrFmt();
        StrFmt *se=parent->findsequencefmt();
 
	Res *r,*r1;
		 
	//int kep=0;while(resn[kep])kep++;
	int kep=strlen(seqn);
	int nlen=strlen(seqngap);
	
	int i;//int ii;
	char pp1[100],pp2[100],pp3[100];
	
	 
	setmatch();se->setmatch();
	//char *dssp=setnewdssp(); 
	for(i=0;i<kep-1;i++) {
		
		r=resn[i];
		r1=resn[i+1];
		if(r==0||r1==0) {
			cerr<<"warning! no residue found for:"<<seqn[i]<<i<<" in the alignment"<<endl;
			continue;	
		}	
		if(r1->id-r->id==1)  continue;

		int n1=compare[i];
		int n2=compare[i+1];	 	
		if(n2-n1>50) continue;
		int ii;
		for(ii=0;ii<kep;ii++) {
			if(TRES.logg) cerr<<ii<<" "<<resn[ii]->name<<" "<<seqn[ii]<<" "<<seqngap[compare[ii]]<<endl;
		}
		//int t1=r1->id-r->id-1;
		//int t2=n2-n1-1;
 		 
		int p1=0;
		int p2=0;
		int d1=0;
		int d2=0;
		int dd1=0;
		int dd2=0; 
		for(ii=n1-5;ii<=n2+5;ii++) {
			if(ii<0||ii>=nlen) continue;
			if(seqngap[ii]!='-') pp1[p1++]=seqngap[ii];
			if(se->seqngap[ii]!='-') pp2[p2++]=se->seqngap[ii];
			if(ii==n1) {d1=p1;dd1=p2;}
			if(ii==n2) {d2=p1;dd2=p2;} 
		}			
		if(d2==d1) d2=d1+1;
		if(dd2==dd1) dd2=dd1+1;	 
		pp1[p1]='\0';pp2[p2]='\0';
		if(p1==p2) continue;
		else if(p1<p2) {
			for(ii=n1-5;ii<=n2+5;ii++) {
				if(ii>=0&&ii<nlen) { 
					seqngap[ii]='-';
					se->seqngap[ii]='-';				 
				}
			}
 
			for(ii=0;ii<p2;ii++) pp3[ii]='-';
			pp3[p2]='\0';
			for(ii=0;ii<d1;ii++) pp3[ii]=pp1[ii];
			for(ii=d2-1;ii<p1;ii++) {
				if(ii<0) continue;
				pp3[ii+p2-p1]=pp1[ii];
			}
 			
			int tt=0;
			for(ii=n1-5;ii<=n2+5;ii++) {
				if(ii<0||ii>=nlen) continue;
				if(tt==p2) break;
				se->seqngap[ii]=pp2[tt];	
			 	seqngap[ii]=pp3[tt];
				tt++; 
			}	
			setmatch();	
			se->setmatch();
		}
		else if(p1>p2) {
			
			for(ii=n1-5;ii<=n2+5;ii++) {
				if(ii>=0&&ii<nlen) { 
					seqngap[ii]='-';
					se->seqngap[ii]='-';				 
				}
			}
 
			for(ii=0;ii<p1;ii++) pp3[ii]='-';
			pp1[p1]='\0';
			for(ii=0;ii<dd1;ii++) pp3[ii]=pp2[ii];
			for(ii=dd2-1;ii<p2;ii++) {
				if(ii<0) continue;
				pp3[ii+p1-p2]=pp2[ii];
			}
 
			int tt=0;
			
			for(ii=n1-5;ii<=n2+5;ii++) {
				if(ii<0||ii>=nlen) continue; 
				if(tt==p1) break;
			 	se->seqngap[ii]=pp3[tt];	
			 	seqngap[ii]=pp1[tt];
				tt++;
			} 
			setmatch();	
			se->setmatch();
		}		
	}

	
}


void StrFmt::treatbreaker() {

	StrFmt *parent=getparentStrFmt();
        StrFmt *se=parent->findsequencefmt();
 
	Res *r,*r1;
		 
	int kep=strlen(seqn);
	int nlen=strlen(seqngap);
	
	int i; 
	//char pp1[100],pp2[100],pp3[100];
	
	Rotate rot; 
	setmatch();se->setmatch();

	//char *dssp=setnewdssp(); 
	reallocatmatch(nlen*2);se->reallocatmatch(nlen*2);
 
	for(i=0;i<kep-1;i++) {
		
		r=resn[i];
		r1=resn[i+1];
		if(r==0||r1==0) {
			cerr<<"warning! no residue found for:"<<seqn[i]<<i<<" in the alignment"<<endl;
			continue;	
		}	
		if(r1->id-r->id==1)  continue;
		
		 
		int n1=compare[i];
		int n2=compare[i+1];	

		int ii;
		int jj=0;
		
		for(ii=n1;ii<=n2;ii++) {
			if((se->seqngap[ii]!='-'&&seqngap[ii]=='-')||(se->seqngap[ii]=='-'&&seqngap[ii]!='-')) {
				jj=1;
				break;
			}  
		}	
		if(jj==1) continue;
		cerr<<"chain breaker at:"<<r->name<<r->oid<<" "<<r1->name<<r1->oid<<endl;
		cerr<<"since there is no insertion/deletion in the chain breaker region"<<endl;
		cerr<<"adjust the local alignment to account for the chain breaker..."<<endl;
		
		int nn=strlen(seqngap);
		if(parent->sitescore[n1]<=parent->sitescore[n2]) {		
			for(ii=nn-1;ii>=n1+1;ii--) {		
				seqngap[ii+1]=seqngap[ii];
				se->seqngap[ii+1]=se->seqngap[ii];
			}
			seqngap[n1+1]='-'; seqngap[nn+1]='\0';
			se->seqngap[n1+1]=se->seqngap[n1];
			se->seqngap[n1]='-';
			se->seqngap[nn+1]='\0';			
			setmatch();se->setmatch();setnewsitescore();		
			kep=strlen(seqn);
			//i++;	
		} 
		else {
			for(ii=nn-1;ii>=n2+1;ii--) {		
				seqngap[ii+1]=seqngap[ii];
				se->seqngap[ii+1]=se->seqngap[ii];
			}
			seqngap[n2+1]='-'; seqngap[nn+1]='\0';
			se->seqngap[n2+1]=se->seqngap[n2];
			se->seqngap[n2]='-';
			se->seqngap[nn+1]='\0';			
			setmatch();se->setmatch();setnewsitescore();		
			kep=strlen(seqn);
			//i++;	
		}				 
	}
	
}

void StrFmt::readModelAlgnFmt(char *f) {

	cerr<<endl;
	cerr<<"reading pir alignment blocks..."<<endl;
	cerr<<endl;

	Strhandler cc;

	if(f==0) return;

	fileName=strdup(f);

	//char **lines = cc.opnfilebylinesimple(f);

	//if(lines==0) return;

	//int n= cc.gettotnum(lines);
	char line[1000],*lines;

	//int i=0;
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

	StrFmt *tmp,*xmp;

	tmp=0;xmp=0;

	int num=0;

	int nl=strlen("#start");

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
		if(strncmp(lines,"#start",nl)==0) {
			if(xmp==0) {
				xmp=this;
			}
			else {
				xmp->next=new StrFmt();
				xmp->next->up=xmp;
				xmp=xmp->next;
			}
			tmp=0;
			continue;
		}
		if(lines[0]=='!'||lines[0]=='#') continue;
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
				tmp->more=new StrFmt();
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
			else if(strncasecmp(lines,"structure:",strlen("structure:"))==0) {
				tmp->source='U';
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
			else if(strncasecmp(lines,"hetero:",strlen("hetero:"))==0) {
				tmp->source='X';
				tmp->token=strdup("hetero");
				cerr<<"do not know the meaning:"<<lines<<endl;
				exit(0);
			}
			else if(strncasecmp(lines,"hetero:",strlen("hetero:"))==0) {
				tmp->source='X';
				tmp->token=strdup("hetero");
				cerr<<"do not know the meaning:"<<lines<<endl;
				exit(0);
			}			
			else if(strncasecmp(lines,"complex:",strlen("complex:"))==0) {
				tmp->source='U';
				tmp->token=strdup("complex");
				cerr<<"do not know the meaning:"<<lines<<endl;
				exit(0);
			}
			else if(strncasecmp(lines,"composite:",strlen("composite:"))==0) {
				tmp->source='U';
				tmp->token=strdup("composite");	
				tmp->numaln=0;
			}
			else if(strncasecmp(lines,"secondary:",strlen("secondary:"))==0) {
				tmp->source='U';
				tmp->token=strdup("secondary");	
				tmp->numaln=0;
			}
			else if(strncasecmp(lines,"sequence:",strlen("sequence:"))==0) {
				tmp->token=strdup("sequence");
			}
			else if(strncasecmp(lines,"sequence-",strlen("sequence-"))==0) {                                
                                tmp->token=strdup("sequence");
				char *s1=strdup(lines);
				char *s2=strchr(s1,':');
				if(s2) *s2='\0';
				char *s3=strchr(s1,'-');
				if(s3) {
					s3++;
					char *s4=strdup(s3);	
					s4=cc.clearchar(s4,"\r\n\t* ");
					if(s4) tmp->zipcode=s4[0];
					if(s4) delete [] s4;
					if(tmp->zipcode=='*'||tmp->zipcode=='-'||tmp->zipcode=='?'||tmp->zipcode=='\0'||tmp->zipcode=='+') {
						cerr<<"wrong: "<<lines<<endl<<endl;
						cerr<<"sequence id could not -,? or non-character"<<endl;
						cerr<<"using a-Z or 0-9"<<endl;
						exit(0);
					}
				}
				if(s1) delete [] s1;
				if(tmp->zipcode=='*'||tmp->zipcode=='-'||tmp->zipcode=='+') {
					cerr<<"\n\n\nthe '*-+' could not be used as id for the sequence pir"<<endl;
					cerr<<lines<<endl;
					exit(0);
				}
                        }
			
			else {
				cerr<<"\n\n\ndoes not know the meaning of:"<<endl;
				cerr<<lines<<endl;
				cerr<<"\n\n\nrejected!!\nit is not either structure or sequence token!"<<endl;
				
				exit(0);
			}

			char **tt=cc.pairbytokensimple(lines,":");

			int to=cc.gettotnum(tt);
			int j;
			int isc=0;
			if(strcmp(tmp->token,"composite")==0) isc=1; 
			 
			for(j=1;j<to&&isc;j++) {
				char *s=strdup(tt[j]);
				s=cc.clearendchar(s," \n\t\r"); 	
				if(s==0||strlen(s)==0) {
					if(s) delete [] s;s=0;
					continue;	
				}	 
				if(j==1) {
                                        tmp->code=strdup(tt[j]);
                                        tmp->code=cc.clearendchar(tmp->code,"\r\t\n ");
                                }
				else {
					
					tmp->resaln[tmp->numaln]=atoi(tt[j]);
					tmp->numaln++;
				}         
				if(s) delete [] s;s=0;                       
			}

			for(j=0;j<to&&isc==0;j++) {
				if(strlen(tt[j])==0) continue;
				if(j==1) {
                                        tmp->code=strdup(tt[j]);
                                        tmp->code=cc.clearendchar(tmp->code,"\r\t\n ");
                                }
                                else if(j==2) {
					tmp->start=atoi(tt[j]);
					char st[1000];
					sprintf(st,"%i",tmp->start);
					if(strstr(tt[j],st)==0) tmp->start=1;
                                }
                                else if(j==3) {
                                        if(tt[j]&&strlen(tt[j])>0) {
                                                //tmp->cid=tt[j][0];
						int ni=strlen(tt[j]);
						int ii;
						for(ii=0;ii<ni;ii++) {
							if(tt[j][ii]==' ') continue;
							tmp->cid=tt[j][ii];
							break;
						}
                                        }
                                }
                                else if(j==4) {
                                        tmp->end=atoi(tt[j]);
					char st[1000];
                                        sprintf(st,"%i",tmp->end);
                                        if(strstr(tt[j],st)==0) tmp->end=1;
                                }
                                else if(j==5) {
                                        if(tt[j]&&strlen(tt[j])>0) {
						
						//char c=tt[j][0];
						char c='0';
						int ni=strlen(tt[j]);
                                                int ii;
                                                for(ii=0;ii<ni;ii++) {
                                                        if(tt[j][ii]==' ') continue;
                                                        c=tt[j][ii];
                                                        break;
                                                }
						if(c!=tmp->cid) {
							cerr<<"the chain id should be the same:"<<endl;
							cerr<<lines<<endl;
							exit(0);
						}
						tmp->cid=c;
                                                //tmp->cid=tt[j][0];
                                        }
                                }

			}
			tt=cc.strdel(tt);
		}
		else if(tmp){
			tmp->seqngap=cc.straddup(tmp->seqngap,lines);
		}
		if(lines) delete [] lines;lines=0;
	}

	
	StrFmt *c,*c0;
  
	//
	int ii=0;
	for(c=this;c;c=c->next) ii++;
	cerr<<endl;
	cerr<<"the number of alignment blocks read: "<<ii<<endl;
	cerr<<endl;

	//uppercase
	cerr<<endl<<"change residue name into uppercase ..."<<endl<<endl;

	Strhandler ccc;
	for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
		if(c0->token==0) {
			cerr<<"no token found:"<<endl;
			c0->printoutonly(stderr);
			exit(0);
		}
		if(c0->token&&strcmp(c0->token,"composite")==0) continue;
		if(c0->token&&strcmp(c0->token,"secondary")==0) continue;
		if(c0->seqngap==0) continue;
		if(c0->token&&strcmp(c0->token,"sequence")==0) { 
			c0->seqngap=ccc.uppercase(c0->seqngap);
		}
		else if(c0->token&&strcmp(c0->token,"structure")==0) { 
			c0->seqngap=ccc.uppercase(c0->seqngap);
		}
	}

	//remove bad residue name
        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
	    if(strcmp(c0->token,"sequence")==0)  {
		//if(c0->cid=='0') c0->cid=' ';
		c0->seqngap=cc.clearchar(c0->seqngap," *\n\t\r");
		if(c0->seqngap&&c0->seqngap[0]=='\0') {
			delete [] c0->seqngap; c0->seqngap=0;
	    	}
	    	if(c0->seqngap) if(TRES.logg>3) cerr<<c0->seqngap<<endl;
	    }
	    else if(strcmp(c0->token,"structure")==0)  {
	    	//if(c0->cid=='0') c0->cid='0';
	    	c0->seqngap=cc.clearchar(c0->seqngap," *\n\t\r");
	    	if(c0->seqngap&&c0->seqngap[0]=='\0') {
			delete [] c0->seqngap; c0->seqngap=0;
	    	}
	    	if(c0->seqngap) if(TRES.logg>3) cerr<<c0->seqngap<<endl;
	    }
	    else if(strcmp(c0->token,"composite")==0)  {
	    	//if(c0->cid=='0') c0->cid='0';
	    	c0->seqngap=cc.clearchar(c0->seqngap," *\n\t\r");
	    	if(c0->seqngap&&c0->seqngap[0]=='\0') {
			delete [] c0->seqngap; c0->seqngap=0;
	    	}
	    	if(c0->seqngap) if(TRES.logg>3) cerr<<c0->seqngap<<endl;
	    }
	    else if(strcmp(c0->token,"secondary")==0)  {
	    	//if(c0->cid=='0') c0->cid='0';
	    	c0->seqngap=cc.clearchar(c0->seqngap," *\n\t\r");
	    	if(c0->seqngap&&c0->seqngap[0]=='\0') {
			delete [] c0->seqngap; c0->seqngap=0;
	    	}
	    	if(c0->seqngap) if(TRES.logg>3) cerr<<c0->seqngap<<endl;
	    }
	    else {
		c0->seqngap=cc.clearchar(c0->seqngap," *\n\t\r");
	    	if(c0->seqngap&&c0->seqngap[0]=='\0') {
			delete [] c0->seqngap; c0->seqngap=0;
	    	}
	    	if(c0->seqngap) if(TRES.logg>3) cerr<<c0->seqngap<<endl;
	    }
        }
	cerr<<endl<<"remove non-standard residue name, treated as gap..."<<endl<<endl;
	removeallnonstandard();
	clearemptyseq();	 
	char *ss=ccc.lastindexof(fileName,"/");
	if(ss&&ss[0]!='\0') {
		char *tt=strdup(ss);
		delete [] fileName;
		fileName=strdup(tt);
		delete [] tt;
	}
}

void StrFmt::printoutallnoseqmesg()
{
	StrFmt *c;
        for(c=this;c;c=c->next) {
		if(strcmp(c->token,"sequence")!=0&&strcmp(c->token,"structure")!=0) continue;
		c->printoutnoseqmesg();
	}
}

void StrFmt::exec()
{
        StrFmt *c;	 
	StrFmt *root=getrootStrFmt();
        for(c=this;c;c=c->next) {
		if(c->mutate==0) {
			cerr<<"ignore the set of alignments:"<<endl;
			c->printout(stderr);
			continue;
		}
		c->mutate->actmutate0();
		if(c->mutate->sharp>=1) {
			if(root->onlyrefine==0) {
			cerr<<endl;
			cerr<<"optimizing initial model of the "<<c->getblockids()<<"th alignment blocks..."<<endl;
			cerr<<endl;
			}
			c->mutate->optimize();			 
		}
	}				 
}
void StrFmt::optimizecom(){

	StrFmt *c;
 
        for(c=this;c;c=c->next) {	
		if(c->iscom==0) continue;	
		if(c->mutate==0) {
			cerr<<"some problem occurs in the set of alignments:"<<endl;
			c->printout(stderr);
			continue;
		}
		if(c->mutate->code&&c->mutate->sharp>=2) {
			cerr<<endl;
                        cerr<<"refinement models... "<<endl;
			cerr<<"the model with code name:"<<c->mutate->code<<endl;
                        cerr<<endl;
			c->mutate->refinement();	
		}
	}
}
void StrFmt::optimize(){

	StrFmt *c;
	StrFmt *root=getrootStrFmt(); 
        for(c=this;c;c=c->next) {		
		if(c->mutate==0) {
			cerr<<"some problem occurs in the set of alignments:"<<endl;
			c->printout(stderr);
			continue;
		}
		
		if(c->mutate->code&&c->mutate->sharp>=2) {
			if(root->onlyrefine==0) {
			cerr<<endl;
                        cerr<<"refinement models of the "<<c->getblockids()<<"th from alignment blocks..."<<endl;
			cerr<<"the model with code name:"<<c->mutate->code<<endl;
                        cerr<<endl;
			}
			if(minfast==0) c->mutate->refinement();	
			else	       c->mutate->smoothmin();
		}
	}
}
void StrFmt::writefinal() {
	StrFmt *c;
	cerr<<endl<<"write out the final models..."<<endl<<endl;
	for(c=this;c;c=c->next) {		
		if(c->mutate==0) {
			cerr<<"ignore the set of alignments:"<<endl;
			c->printout(stderr);
			continue;
		}
		if(c->mutate->code) c->mutate->writefinalout();	
	}
}

void StrFmt::beforewritefinal() {
        StrFmt *c;
	setupatmid0();
        cerr<<endl<<"check out the final models..."<<endl<<endl;
        for(c=this;c;c=c->next) {
                if(c->mutate==0) {
                        cerr<<"ignore the set of alignments:"<<endl;
                        c->printout(stderr);
                        continue;
                }
                if(c->mutate->code) c->mutate->beforewritefinalout();
        }
}

void StrFmt::setupatmid0() {
        StrFmt *c;
        for(c=this;c;c=c->next) {
                if(c->mutate->code) c->mutate->mpdb->configureatmid0();
        }
}

void StrFmt::setshift(float dcut,int nn) {

	StrFmt *c;
        for(c=this;c;c=c->next) {		
		StrFmt *parent=c;
		Mutate *cc=c->mutate;
		if(cc->sharp<5) continue;
		if(parent->other) delete parent->other;parent->other=0;
		parent->setotherfmt(1,dcut,nn);
		parent->setotherfmt(2,dcut,nn);
		if(parent->other==0) continue;
		if(TRES.logg>3) parent->other->printoutall(stderr);
		parent->other->prepare();
		 
		StrFmt *tmp;
		for(tmp=parent->other;tmp;tmp=tmp->next) {
			if(tmp->mutate==0)  continue;				
			tmp->mutate->seglen=cc->seglen;
			tmp->mutate->asloop=cc->asloop;
			tmp->mutate->flag=c->mutate->flag;
			tmp->mutate->window=cc->window;
			tmp->mutate->refineid=cc->refineid;
			tmp->mutate->sharp=cc->sharp;
			tmp->mutate->fapr=cc->fapr;
			tmp->mutate->rmsd=cc->rmsd;				
			tmp->mutate->actmutate0();
			if(tmp->mutate->sharp>=1)tmp->mutate->optimize();
		} 			 
	}
}

void StrFmt::createComStrFmt() {
 
	StrFmt *c;
	cerr<<endl<<"handling blocks for composite pirs..."<<endl<<endl;
	for(c=this;c;c=c->next) {
		if(c->iscom) continue;
		if(c->mutate==0) continue;
		StrFmt *t;
		StrFmt *se=c->getsequenceStrFmt(); 
		for(t=c;t;t=t->more) {
			if(strcmp(t->token,"composite")) continue;
			StrFmt *c0=new StrFmt(c->mutate->owner);
			StrFmt *se0=new StrFmt(se);
			StrFmt *com0=new StrFmt(t);
			//
			StrFmt *gg=this; while(gg->next) gg=gg->next;	
			c0->up=gg;
			gg->next=c0;
			//
			c0->more=se0;
			se0->more=com0;
			com0->more=0;
			com0->last=se0;
			se0->last=c0;
			//
			if(c0->sitescore==0) {
				if(se0->sitescore) {
					c0->sitescore=se0->sitescore;
					se0->sitescore=0;
				}
				else if(com0->sitescore) {
					c0->sitescore=com0->sitescore;
					com0->sitescore=0;
				}
			}
			//
			c0->mutate=new Mutate(c->mutate);
			c0->mutate->owner=c0;
			c0->mutate->autotmp=0;
			if(c0->mutate->code) delete [] c0->mutate->code;
			c0->mutate->code=0;
			c0->mutate->code=strdup(t->code);
			c0->iscom=1;
			//			 
		}
	}
}

 

void StrFmt::combine()
{
        StrFmt *c;
	//Mutate *cc[10000]; 
	StrFmt *root=getrootStrFmt();
        for(c=this;c;c=c->next) {
		if(c->iscom==0) continue;
		StrFmt *se = c->getsequenceStrFmt();		
		if(c->mutate==0) {
			cerr<<"ignore the set of alignments:"<<endl;
			c->printout(stderr);
		}
		if(se==0) {
			cerr<<"\n\n\nno sequence pir found:"<<endl<<endl;
			c->printout(stderr);
		}
		if(se->iscomposite()==0) continue;		
		c->mutate->combine();
		if(c->mutate->segen) c->mutate->segen->mutate=0;
		if(c->mutate->sharp>=1&&0) { //modified on 07/22/2002 to take out automatic refinement
                        if(root->onlyrefine==0) {
                        cerr<<endl;
                        cerr<<"optimizing initial composite model"<<endl;
                        cerr<<endl;
                        }
			int cc=c->mutate->sharp;
			c->mutate->sharp=4;
                        c->mutate->refineminment();
			c->mutate->sharp=cc;
			if(c->mutate->out>0&&c->mutate->code) {
                		char line[100];
                		sprintf(line,"%s_composite_min.%i.pdb",c->mutate->code,c->mutate->refineid++);
                		cerr<<endl<<"output the minimized composite structure: "<<line<<endl<<endl;
                		c->mutate->mpdb->write(line);
        		}
                }
		
	}
}


void StrFmt::setmutate(char *ss,int nn)
{
        StrFmt *c;
        for(c=this;c;c=c->next) {
		//if(strcmp(ss,"shiftcut")==0) c->shiftcut=nn/10;
		if(c->mutate==0) {
			cerr<<"some problem occurs in the set of alignments:"<<endl;
			c->printout(stderr);
		}
		if(strcmp(ss,"fapr")==0) c->mutate->fapr=nn;		 		 
		//if(strcmp(ss,"seed")==0) c->mutate->ranseed=nn;
		if(strcmp(ss,"test")==0) c->mutate->test=nn;
		if(strcmp(ss,"sharp")==0) c->mutate->sharp=nn;
		if(strcmp(ss,"out")==0) c->mutate->out=nn;
		if(strcmp(ss,"asloop")==0) c->mutate->asloop=nn;
		if(strcmp(ss,"rmsd")==0) c->mutate->rmsd=nn/100;
		if(strcmp(ss,"seglen")==0) c->mutate->seglen=nn;
		if(strcmp(ss,"restraint")==0) c->mutate->restraint=nn;		
	}
}



void StrFmt::treatcomposite() {
	StrFmt *t=getStrFmt("composite");
	if(t==0) return;
	StrFmt *se=getsequenceStrFmt();
	if(se==0) return;
	if(t->seqngap) return;
	
	int nlen=strlen(se->seqngap);
	t->seqngap=new char[nlen+10];
	int i;
	for(i=0;i<nlen;i++) {
		t->seqngap[i]='-';
	}
	t->seqngap[nlen]='\0';
}

int StrFmt::getblockids(){

	StrFmt *root=getrootStrFmt();
	StrFmt *t;
	StrFmt *parent=getparentStrFmt();
	int i=0;
	for(t=root;t;t=t->next) {
		if(t==parent) return i;
		i++;
	}
	return -1;
}

void StrFmt::checknonesegmentcomposite() {

	StrFmt *root=getrootStrFmt();
	StrFmt *t,*t0;
	for(t=root;t;t=t->next) 
	for(t0=t;t0;t0=t0->more) {
		if(strcmp(t0->token,"composite")!=0) continue;
		StrFmt *se=t->getStrFmt("sequence");	
		if(se==0) continue;
		int i=0;
		while(1) {
			i++;
			int *stem=t0->getresalncomposite(i);
			if(stem==0) break;
			int j;
			int m=0;
			for(j=stem[0];j<=stem[1];j++) {
				if(se->seqngap[j]!='-') m++;
			}
			
			if(m==0) {
				cerr<<endl;
				cerr<<"the "<<i<<"th segment in the composite does not correspond to any subsequence"<<endl;
				t0->printoutonly(stderr);
				cerr<<endl;
				cerr<<"try to ignore this composite segment"<<endl;
				for(j=stem[0];j<=stem[1];j++) {	
					t0->seqngap[j]=se->zipcode;
				}
				int mm=0;
				for(j=0;j<t0->numaln;j++) {
					if(j==2*i-2||j==2*i-1) continue;
					t0->resaln[mm++]=t0->resaln[j];
				}
				t0->numaln=mm;
			}
			delete [] stem;stem=0;								 
		} 	
	}
}
void StrFmt::reindexcomposite() {

	StrFmt *root=getrootStrFmt();
	StrFmt *t,*t0;
	for(t=root;t;t=t->next) 
	for(t0=t;t0;t0=t0->more) {
		if(strcmp(t0->token,"composite")) continue;
		StrFmt *te=t0; 		
		StrFmt *se=t->getsequenceStrFmt();
		if(se==0) continue;
		int nlen=strlen(te->seqngap);
		char *tmp=strdup(te->seqngap);

		int i;
		for(i=0;i<nlen;i++)  tmp[i]=se->zipcode;
	 
		int n=0;
		for(i=0;i<nlen;i++) {
			if(se->oldgap[i]=='-') continue;
			int m=se->compare[n];
			tmp[m]=te->seqngap[i];
			n++;
			if(tmp[m]==' '||tmp[m]=='-') {				 
				continue;
			}
			else {
				StrFmt *fd=getStrFmtWithId(tmp[m]);
				if(fd==0) {
					cerr<<"\n\n\nthere is no sequence with id:"<<tmp[m]<<endl;
					cerr<<"\ncheck: "<<endl;
					te->printoutonly(stderr);
					exit(0);
				}					
			}		
		}
		if(te->seqngap) delete [] te->seqngap;
		te->seqngap=tmp;
		if(TRES.logg>3) cerr<<te->seqngap<<endl;
	}
	
}



int StrFmt::iscomposite() {
	
	//only one alignment
	StrFmt *root=getrootStrFmt();
	StrFmt *t;
	int j=0;
	for(t=root;t;t=t->next) {
		j++;
	}
	if(j<2) return 0;
	
 
	t=getStrFmt("composite");
	StrFmt *se=getsequenceStrFmt();
	if(t==0) return 1;
	int n=strlen(se->seqngap);

	int i;
	
	for(i=0;i<n;i++){
		if(se->seqngap[i]=='-') continue;
		if(t->seqngap[i]!=se->zipcode) return 1; 
	}
	return 0;
}

void StrFmt::printoutnoseqmesg()
{
        StrFmt *c;

        for(c=this;c;c=c->more) {
		if(strcmp(c->token,"sequence")==0) return;
	}

	cerr<<"Warning! some alignments does not have query"<<endl;
}


void StrFmt::clearemptyseq() {

        StrFmt *c,*c0;

        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
		if(strcmp(c0->token,"sequence")!=0&&strcmp(c0->token,"structure")!=0) continue;
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


void StrFmt::deletegap() {

	if(seqngap==0) return;

	if(seqn) delete [] seqn;seqn=0;	
	if(match) delete [] match;match=0;
	if(compare) delete [] compare;compare=0;

	int n=strlen(seqngap);

	seqn=new char[n+1];
	match=new int[n+1];
	compare=new int[n+1];

	int i;

	for(i=0;i<n+1;i++) match[i]=-1;
	for(i=0;i<n+1;i++) seqn[i]='\0';
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

void StrFmt::deleteallgap() {

	StrFmt *t;
	StrFmt *tt;

	for(tt=this;tt;tt=tt->next)
	for(t=tt;t;t=t->more) {
		if(strcmp(t->token,"sequence")==0) {
			t->deletegap();
		}
		else if(strcmp(t->token,"structure")==0)  {
			t->deletegap();
		}
	}
}

void StrFmt::takeofzero() {

        StrFmt *t;
        StrFmt *tt;

        for(tt=this;tt;tt=tt->next)
        for(t=tt;t;t=t->more) {
		if(t->token==0) token=strdup("UNKNOWN");
        }
}

StrFmt *StrFmt::findfirstsequencefmt() {

        StrFmt *t;
        StrFmt *tt;

        for(tt=this;tt;tt=tt->next)
        for(t=tt;t;t=t->more) {
                if(strcmp(t->token,"sequence")==0) return t;
        }

	return 0;
}

void StrFmt::removeallnonstandard() {

	StrFmt *c,*c0;

        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
		if(strcmp(c0->token,"sequence")==0||strcmp(c0->token,"structure")==0) {
        		c0->removenonstandard();
		}
        }
}

void StrFmt::removenonstandard() {

	//if(s==0) return;
	if(token==0) return;
	if(strcmp(token,"sequence")!=0&&strcmp(token,"structure")!=0) return;

	if(seqngap==0) return;

	int n=strlen(seqngap);

	Tres *tt;

	for(int i=0;i<n;i++) {
		if(seqngap[i]=='-') continue;
		tt=TRES[seqngap[i]];
		if(tt==0)   {
			cerr<<"the residue :"<<i<<" "<<seqngap[i]<<" is not standard! treated as dash!"<<endl;
			seqngap[i]='-';
		}
	}
	return;
}

void StrFmt::createstructure(char *s) {

        pdb=new Pdb();
        pdb->chn=new Chn();
        pdb->chn->create(s);
        pdb->configure();
	pdb->chn->header();
	pdb->chn->start=start;
	if(cid!='0'&&cid!='1'&&cid!='\0') pdb->chn->id=cid;
}

void StrFmt::createstructure() {

	createstructure(seqn);
}
 
void StrFmt::setseqnstructure() {

	StrFmt *c,*c0,*c1;

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

StrFmt *StrFmt::getlongestsequence() {

	StrFmt *c,*c0,*c1;

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


void StrFmt::setseqnstructure(Pdb *p) {

	StrFmt *c,*c0;

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

void StrFmt::setalldefaultseq() {

	StrFmt *c,*c0;

        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
                if(strcmp(c0->token,"structure")==0) c0->setdefaultstructseq();
		else if(strcmp(c0->token,"sequence")==0) c0->setdefaultqueryseq();
        }
}

void StrFmt::setdefaultqueryseq() {

	StrFmt *c=getrootStrFmt();

	int n=-1000;

        StrFmt *t,*t0,*g=0;

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

void StrFmt::setdefaultstructseq() {

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
void StrFmt::getstructure() {

	if(pdb==0) {
		pdb=new Pdb();		
		pdb->read(code,cid);
	}
	if(pdb->chn==0) {
		cerr<<"the protein structure "<<code<<" "<<cid<<" does not exist"<<endl;
		exit(0);
	}
	
	//assign pdbs
	StrFmt *root=getrootStrFmt();
	StrFmt *t,*c;
	for(t=root;t;t=t->next)
	for(c=t;c;c=c->more){
		if(strcmp(t->token,"structure")!=0) continue;
		if(c->code==0) continue;
		if(c->pdb) continue;
		if(strcmp(code,c->code)!=0) continue;	
		if(cid!=c->cid) continue;	
		c->pdb=new Pdb(pdb);
		c->pdb->configure();
	}
		

	//check default sequence and deletegaps 
	if(strcmp(token,"structure")==0&&seqngap==0) {
		Chn *c=pdb->ischain(cid);
		if(c==0) c=pdb->chn;
		seqngap=c->getseqn();
		StrFmt *s=getStrFmt("sequence");
		if(s&&s->seqngap==0) s->seqngap=pdb->chn->getseqn();
		s=getrootStrFmt();
		s->deleteallgap();		 
	}

	//delete all other chains	 
	deletechain();

	//complex or hetero
	if(strcmp(token,"structure")!=0) return;	
	
	root=getrootStrFmt();
	/*
	if(root->addh) { 
		pdb->chn->header();
		pdb->chn->addhatoms(1);
		pdb->configure();
	}
	*/
	//set up structure
	if(pdb) {
		pdb->checkpdb(0);
		PdbFix pdbfix;
		int n=pdbfix.checkpdb(pdb);
		if(n) {
			cerr<<"the template structure has missing atoms..."<<endl;
			cerr<<"try profix program to fix it..."<<endl;
			pdbfix.pdb=pdb;
			pdbfix.onlybackbone=0;
			pdbfix.pdb->configure();
			pdbfix.ready();
			pdbfix.myfixnow();
			pdb=pdbfix.pdb;
			pdbfix.pdb=0;
			if(TRES.logg>3) pdb->write("fix.pdb");
		}
	}
	 if(root->addh) {
                pdb->chn->header();
                pdb->chn->addhatoms(1);
                pdb->configure();
        }

	setresn();
}


void StrFmt::getallstructure() {

	StrFmt *c,*c0;

	for(c=this;c;c=c->next)
	for(c0=c;c0;c0=c0->more) {
		if(strcmp(c0->token,"structure")==0) {
			c0->getstructure();
		}		
		else if(strcmp(c0->token,"complex")==0) {
			c0->getstructure();
		}	
		else if(strcmp(c0->token,"hetero")==0) {
			c0->getstructure();
		}
	}
}


void StrFmt::deleteallchain() {

	StrFmt *c,*c0;

        for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {

                if(strcmp(c0->token,"structure")==0) {
                        c0->deletechain();
                }
        }
}

void StrFmt::deletechain() {

	if(strcmp(token,"hetero")!=0) {
		if(pdb==0||pdb->chn==0||pdb->chn->next==0) return;
	}
	else {
		return; //change later on
	}
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
		if(c->seqcard==0) c->seqcard=strdup(" ");
		Algn aa;
		aa.setsequence(c->seqcard,seqn);
		aa.defalgn();
		float d=aa.calctotalscore();
		if(c->id==cid&&profix==1) d=100000000;
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

StrFmt * StrFmt ::findsequencefmt() {

	StrFmt *c;

	for(c=this;c;c=c->more) {

		if(c->token&&strcmp(c->token,"sequence")==0) return c;
		
	}

	return 0;
}

StrFmt * StrFmt ::findstructurefmt() {

        StrFmt *c;

        for(c=this;c;c=c->more) {

                if(c->token&&strcmp(c->token,"structure")==0) return c;

        }

        return 0;
}

void StrFmt::setmutatesec() {

	Pdb *mutatepdb=0;
	StrFmt *c;

        for(c=this;c;c=c->next) {
                if(strcmp(c->token,"structure")) continue;
                Chn *cc;
		//for(cc=c->mutatepdb->chn;cc;cc=cc->next) {
		for(cc=mutatepdb->chn;cc;cc=cc->next) {
			cc->header();
			cc->buildhbond();
			cc->setdsspstr();
			
			cc->write("sss");
			Res *r;
                	for(r=cc->res;r;r=r->next) if(TRES.logg) cerr<<r->id<<r->name<<" "<<r->sec<<endl;
                	for(r=cc->res;r;r=r->next) {
                        	for(HBondList *h=r->hbond;h;h=h->next) {
                                	if(TRES.logg) cerr<<h->donor->res->name<<" "<<h->donor->res->id<<" "<<h->donor->name<<" ";
                                	if(TRES.logg) cerr<<h->acceptor->res->name<<" "<<h->acceptor->res->id<<" "<<h->acceptor->name<<endl;
				}
                        }
		}
        }
}

void StrFmt::setmutatesidepdb() {

	StrFmt *c;

	for(c=this;c;c=c->next) {
		if(strcmp(c->token,"structure")) continue;
		c->setmutatesidechainpdb();
	}
}

StrFmt *StrFmt::getparentStrFmt() {

	StrFmt *c;

        for(c=this;c->last;c=c->last);

	return c;
}

StrFmt *StrFmt::getrootStrFmt() {

        StrFmt *c,*c0;
        for(c=this;c->last;c=c->last);
	for(c0=c;c0->up;c0=c0->up);
        return c0;
}


StrFmt *StrFmt::getsequenceStrFmt() {

        StrFmt *c,*a;

	c=getparentStrFmt();

	for(a=c;a;a=a->more) {
		if(strcmp(a->token,"sequence")==0) return a;
	}
        return 0;
}
StrFmt *StrFmt::getStrFmtWithId(char nm) {

        StrFmt *c,*a;
	StrFmt *root=getrootStrFmt();
 
	for(a=root;a;a=a->next) 
	for(c=a;c;c=c->more) {
		if(c->zipcode==nm) return c;
	}
        return 0;
}

StrFmt *StrFmt::getStrFmt(char *nm) {

        StrFmt *c,*a;

	c=getparentStrFmt();

	for(a=c;a;a=a->more) {
		if(strcmp(a->token,nm)==0) return a;
	}
        return 0;
}

int StrFmt::iscodeexist(char *nm) {

        StrFmt *c,*a;
	if(nm==0) return 0;
	c=getparentStrFmt();
	
	for(a=c;a;a=a->more) {
		if(a->code&&strcmp(a->code,nm)==0) return 1;
	}
        return 0;
}

void StrFmt::checklength() {

	StrFmt *c;

        for(c=this;c;c=c->next) {
                c->checkalignmentlength();
        }
}



void StrFmt::checkalignmentlength() {

	StrFmt *c;

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


void StrFmt::setmutatesidechainpdb() {
	Pdb *mutatepdb=0;
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

	StrFmt *parent=getparentStrFmt();
	StrFmt *se = parent->getsequenceStrFmt();
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
	//mutatepdb->write("mutated.out1");
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

void StrFmt::treatallbreaker(){
	StrFmt *t,*t0;

	cerr<<endl<<"find out chain breakers and resolve..."<<endl<<endl;
	for(t=this;t;t=t->next) {
		for(t0=t;t0;t0=t0->more) {
			if(strcmp(t0->token,"structure")==0) {
				t0->setmatch();
				t0->pdb->chn->header();
				t0->treatbreaker();			
			} 
		}
	}	
}

void StrFmt::addhatoms() {

	StrFmt *t,*t0;

	for(t=this;t;t=t->next) {
		for(t0=t;t0;t0=t0->more) {
			if(strcmp(t0->token,"structure")==0) {
				if(t0->pdb==0||t0->pdb->chn==0) continue;
				Chn *c;
				for(c=t0->pdb->chn;c;c=c->next)
				c->addhatoms(1);
				t0->pdb->configure();
			} 
		}
	}	
}

void StrFmt::reorderstrftm(){
//obsolete
	StrFmt *t,*t0,*r,*r0;
	
	re100:
	r=0; 
	int n=0; 	
	for(t=this;t;t=t->next) {
		r0=0;			
		for(t0=t;t0;t0=t0->more) {
			if(strcmp(t0->token,"structure")!=0) {
				r0=t0;
				continue;
			}
			if(r0&&r) {
				r->next=t0;
				t0->next=t->next;
				t->next=0;
				r0->more=t0->more;
				t0->more=t;	
				n++;		
			} 			
			if(n) goto re200;
			r0=t0;
		}
		r=t;
	}

	re200:

	if(n) goto re100;

}

void StrFmt::preparesuperimpose() {
	//reorder so that structure first
	 
	checklegalsuperimpose();
	deleteallgap();
	getallstructure();
	checklegalsuperimpose();
	StrFmt *t,*t0;
	int n=0;
	for(t=this;t;t=t->next) {
		for(t0=t;t0;t0=t0->more) {
			n++;
			if(strcmp(t0->token,"structure")==0) t0->setmatch();
		}
	}
	if(n!=2) {
		cerr<<"only two pir alignments allowed!"<<endl;
		exit(0);
	}
	if(pdb==0||pdb->chn==0||pdb->chn->res==0) {
		cerr<<"no residue found for the first alignmetn structure"<<endl;
		exit(0);
	}
	if(more->pdb==0||more->pdb->chn==0||more->pdb->chn->res==0) {
		cerr<<"no residue found for the second alignmetn structure"<<endl;
		exit(0);
	}
}
void StrFmt::writeseqxyz() {	
	checklegalseqxyz();
	checkseqname();
	deleteallgap();
	
	StrFmt *t,*t0;
	for(t=this;t;t=t->next) {
		for(t0=t;t0;t0=t0->more) {
			if(t0->code==0) continue;
			t0->createstructure();
			char line[1000];
			if(t0->pdb==0) continue;
			char *s=pdb->getname(t0->code);			 
			if(s==0) continue;
			sprintf(line,"%s.pdb",s);
			if(s) delete [] s;s=0;
			t0->pdb->write(line);
		}
	}
	
}

 

void StrFmt::structuresuperimpose(char *line,int flg) {

	Stralg algn;
	if(line==0) return;	
	int n=max(pdb->manyres(),more->pdb->manyres())+1000;
	algn.flag=1;
	algn.alga=new Res*[n];
	algn.algb=new Res*[n];

	int i=0;
	for(i=0;i<n;i++) {algn.alga[i]=0; algn.algb[i]=0;}
	 
	int nlen=strlen(seqngap);
	int m=0;
	for(i=0;i<nlen;i++) {
		if(seqngap[i]=='-'||more->seqngap[i]=='-') continue;
		int j1=match[i];
		int j2=more->match[i];
		if(j1==-1||j2==-1) continue;
		if(resn[j1]==0||more->resn[j2]==0) continue;
		algn.alga[m]=resn[j1];
		algn.algb[m]=more->resn[j2];
		m++;
	}

	algn.superimpose(flg);
	Strhandler cc;
	//pdb->chn->id='A';
	pdb->info=cc.strdel(pdb->info);
	pdb->endinfo=cc.strdel(pdb->endinfo);
	//more->pdb->chn->id='B';
	more->pdb->info=cc.strdel(more->pdb->info);
        more->pdb->endinfo=cc.strdel(more->pdb->endinfo);
	pdb->chn->next=more->pdb->chn;
	pdb->writeold(line);
	pdb->chn->next=0;
}

void StrFmt::setdefaultseqn(){

	StrFmt *t;
	//set default structure/sequence
	for(t=this;t;t=t->next) {
		StrFmt *seq=t->getStrFmt("sequence");
		StrFmt *str=t->getStrFmt("structure");
		if(str==0||seq==0) continue;
		if(seq->seqngap==0&&str->seqngap) {
			seq->seqngap=strdup(str->seqngap);
			cerr<<"assign the sequence to the structure pir from the sequence pir"<<endl;
		}
		else if(seq->seqngap&&str->seqngap==0) {
			str->seqngap=strdup(seq->seqngap);
			cerr<<"assign the sequence to the sequnece pir from the structure pir"<<endl;
		}
	}
}

int *StrFmt::getresalncomposite(int ii) {
	
	StrFmt *t,*s;
	
	for(t=this;t;t=t->next) {
		for(s=t;s;s=s->more) {
			if(strcmp(s->token,"composite")!=0) continue;
			if(s->seqngap==0) continue;
			int n=strlen(s->seqngap);
			StrFmt *se=t->getStrFmt("sequence");
			if(se==0) continue;
			int i;
			char c='\0';
			int j=0;
			for(i=0;i<n;i++) {
				char v=s->seqngap[i];
				if(v==' '||v==se->zipcode) {
					c=v;
					continue;
				}
				if(v==c) continue;
				c=v;
				j++;
				if(j==ii) {
					int i1;
					int *stem=new int[2];
					stem[0]=i;
					for(i1=i;i1<n;i1++) {
						v=s->seqngap[i1];
						if(v==c) stem[1]=i1;
						else break;
					}
					return stem;
				}
			}		
			
		}
	}
	return 0;
}
void StrFmt::checkresalncomposite() {

	StrFmt *t,*s;
	
	for(t=this;t;t=t->next) {
		for(s=t;s;s=s->more) {
			if(strcmp(s->token,"composite")!=0) continue;
			if(s->seqngap==0) continue;
			int n=strlen(s->seqngap);
			StrFmt *se=t->getStrFmt("sequence");
			if(se==0) continue;
			int i;
			char c='\0';
			int j=0;
			for(i=0;i<n;i++) {
				char v=s->seqngap[i];
				if(v==' '||v==se->zipcode) {
					c=v;
					continue;
				}
				if(v==c) continue;
				c=v;
				j++;
			}	
			cerr<<endl;	
			cerr<<"the number of segments possibly to be replaced: "<<j<<endl;
			if(j*2!=s->numaln) {
				cerr<<endl;	
				cerr<<"error in specifying conformations alignment boundaries"<<endl;
				cerr<<"require "<<j*2<<" numbers in token line for all the segments"<<endl;
				cerr<<"however forund "<<s->numaln<<" numbers in token line"<<endl;
				cerr<<"each segment possibly to be replaced needs two numbers (start and end)"<<endl;
				s->printoutonly(stderr); 			
				exit(0);
			}
			cerr<<endl;		
		}
	}

}
void StrFmt::prepare() {
	if(onlyrefine==0) {
		cerr<<endl<<"check error in alignment blocks..."<<endl<<endl; 
	}
	checkresalncomposite();
	checklegal();
	setdefaultseqn(); 
	checklegal();
	checkseqcomname();
	deleteallgap();
	//get structure and check
	getallstructure();	
	checklegal();
	//
	setsecstruct();
	setalloldgap();
	setallseqnresn();
	//find zero segments
	if(onlyrefine==0) {
		cerr<<endl<<"check error in composite pirs..."<<endl;
	}
	checknonesegmentcomposite();
	setalloldgap();	
	//tuning the sequence
	
	tuneallalign();
	setnewsitescore();//to handle treatbreaker 
	treatallbreaker();
	setnewsitescore();
	setnest();
	//if tuned, then need to reasses
	//reindexcomposite();

	if(onlyrefine>0) return;
	if(alnback==0) return;
	char line[1000];
	sprintf(line,"%s~",fileName);
	FILE *fp=fopen(line,"w");
	cerr<<endl<<"write down the alignment blocks actually used in the program to: "<<line<<endl<<endl;
	StrFmt *t,*t0;
	for(t=this;t;t=t->next) {
		fprintf(fp,"#start\n");
		for(t0=t;t0;t0=t0->more) {
			fprintf(fp,">P1;%s\n",t0->code);
			fprintf(fp,"%s:%s:%c:%i:%c:%i\n",t0->token,t0->code,t0->cid,t0->start,t0->cid,t0->end);	
			fprintf(fp,"%s\n",t0->seqngap);
		}	
		StrFmt *se = t->findstructurefmt();
		if(se==0) continue;
		char *dssp=se->setnewdssp();
		if(dssp) {
			fprintf(fp,">P1;%s; secondary structure\n",se->code);
			fprintf(fp,"%s:%s:%c:%i:%c:%i\n",se->token,se->code,se->cid,se->start,se->cid,se->end);
			fprintf(fp,"%s\n",dssp);
			if(dssp) delete [] dssp;dssp=0;
		}
		
		fprintf(fp,"#end\n");
	}
	fclose(fp);
	
	for(t=this;t;t=t->next) {
		StrFmt *se = t->findstructurefmt();
		if(se==0) continue;
		int n=strlen(se->seqngap);
		for(int i=0;i<n;i++) {
			int j=se->match[i];
			if(j==-1) continue;
			if(TRES.logg) cerr<<i<<" "<<j<<" "<<se->resn[j]->name<<" "<<se->seqngap[i]<<endl;
		}
	}
}

void StrFmt::setseqnlimit(int lent) {

	StrFmt *s;
	for(s=this;s;s=s->next) {
		StrFmt *parent=getparentStrFmt();
		StrFmt *se=parent->findsequencefmt(); 
		StrFmt *tr=parent->findstructurefmt();
		if(se==0||tr==0) continue;
		int n=strlen(se->seqngap);
		int i,j;
		for(i=0;i<n;i++) {
			if(se->seqngap[i]=='-'||tr->seqngap[i]!='-') continue;
			int nn=0;
			for(j=i;j<n;j++) {
				if(se->seqngap[j]=='-'&&tr->seqngap[j]=='-') continue;
				if(se->seqngap[j]!='-'&&tr->seqngap[j]=='-') {
					nn++;
				}
				else {
					break;
				}				
			}
			if(nn<=lent) {
				i=j;
				continue;
			}
			else {
				int ii;
				for(ii=i;ii<j;ii++) {
					se->seqngap[i]=tr->seqngap[i];					
				}
				i=j;
			}
		}
		se->setmatch();
		tr->setmatch();
	}
}

void StrFmt::prepareprofix(int lent) {
	if(onlyrefine==0) {
		cerr<<endl<<"check error in alignment blocks..."<<endl<<endl; 
	}
	checkresalncomposite();
	checklegal();
	setdefaultseqn(); 
	checklegal();
	checkseqcomname();
	deleteallgap();
	//get structure and check
	getallstructure();
	setseqnlimit(lent);	
	checklegal();
	//
	setsecstruct();
	setalloldgap();
	setallseqnresn();
	//find zero segments
	if(onlyrefine==0) {
		cerr<<endl<<"check error in composite pirs..."<<endl;
	}
	checknonesegmentcomposite();
	setalloldgap();	
	//tuning the sequence
	
	tuneallalign();
	setnewsitescore();//to handle treatbreaker 
	treatallbreaker();
	setnewsitescore();
	setnest();
	//if tuned, then need to reasses
	//reindexcomposite();

	if(onlyrefine>0) return;
	if(alnback==0) return;
	char line[1000];
	sprintf(line,"%s~",fileName);
	FILE *fp=fopen(line,"w");
	cerr<<endl<<"write down the alignment blocks actually used in the program to: "<<line<<endl<<endl;
	StrFmt *t,*t0;
	for(t=this;t;t=t->next) {
		fprintf(fp,"#start\n");
		for(t0=t;t0;t0=t0->more) {
			fprintf(fp,">P1;%s\n",t0->code);
			fprintf(fp,"%s:%s:%c:%i:%c:%i\n",t0->token,t0->code,t0->cid,t0->start,t0->cid,t0->end);	
			fprintf(fp,"%s\n",t0->seqngap);
		}	
		StrFmt *se = t->findstructurefmt();
		if(se==0) continue;
		char *dssp=se->setnewdssp();
		if(dssp) {
			fprintf(fp,">P1;%s; secondary structure\n",se->code);
			fprintf(fp,"%s:%s:%c:%i:%c:%i\n",se->token,se->code,se->cid,se->start,se->cid,se->end);
			fprintf(fp,"%s\n",dssp);
			if(dssp) delete [] dssp;dssp=0;
		}
		
		fprintf(fp,"#end\n");
	}
	fclose(fp);
	
	for(t=this;t;t=t->next) {
		StrFmt *se = t->findstructurefmt();
		if(se==0) continue;
		int n=strlen(se->seqngap);
		for(int i=0;i<n;i++) {
			int j=se->match[i];
			if(j==-1) continue;
			if(TRES.logg) cerr<<i<<" "<<j<<" "<<se->resn[j]->name<<" "<<se->seqngap[i]<<endl;
		}
	}
}

void StrFmt::setallseqnresn() {

	StrFmt *t,*t0;

	for(t=this;t;t=t->next) {
		for(t0=t;t0;t0=t0->more) {
			if(strcmp(t0->token,"structure")==0) {
				t0->setseqnresn();		
				t0->setmatch();
			}
			/*
			if(strcmp(t0->token,"sequence")==0) {
				t0->setseqnresn();		
				t0->setmatch();
			}
			*/
		}
	}	
}

void StrFmt::printoutall(FILE *fp) {

	StrFmt *t0;
	for(t0=this;t0;t0=t0->next) {
		if(t0->code) {
			t0->printout(fp);
		}
	}
}
void StrFmt::printout(FILE *fp) {

	StrFmt *t,*t0;
	fprintf(fp,"#start\n");
	t=this;
        for(t0=t;t0;t0=t0->more) {
              fprintf(fp,">P1;%s\n",t0->code);
              if(strcmp(t0->token,"composite")!=0) {
			fprintf(fp,"%s:%s:%c:%i:%c:%i\n",t0->token,t0->code,t0->cid,t0->start,t0->cid,t0->end);
	      }
	      else {
			
			fprintf(fp,"%s:%s",t0->token,t0->code);
			int n;
			for(n=0;n<t0->numaln;n++) {
				fprintf(fp,":%i",t0->resaln[n]);
			}
			fprintf(fp,"\n");
	      }
              fprintf(fp,"%s\n",t0->seqngap);
        }       
        fprintf(fp,"#end\n"); 
}
 
void StrFmt::printoutonly(FILE *fp) {

	StrFmt *t0;
	//fprintf(fp,"#start\n");
        t0=this;
              fprintf(fp,">P1;%s\n",t0->code);
	      if(t0->token==0||strcmp(t0->token,"composite")!=0) {
			fprintf(fp,"%s:%s:%c:%i:%c:%i\n",t0->token,t0->code,t0->cid,t0->start,t0->cid,t0->end);
	      }
	      else {
			
			fprintf(fp,"%s:%s",t0->token,t0->code);
			int n;
			for(n=0;n<t0->numaln;n++) {
				fprintf(fp,":%i",t0->resaln[n]);
			}
			fprintf(fp,"\n");
	      }
              fprintf(fp,"%s\n",t0->seqngap);
              //if(t0->iscom==0) fprintf(fp,"%s:%s:%c:%i:%c:%i\n",t0->token,t0->code,t0->cid,t0->start,t0->cid,t0->end);
	      //else fprintf(fp,"%s:%s:%i:%i\n",t0->token,t0->code,t0->start,t0->end);
              //fprintf(fp,"%s\n",t0->seqngap);
          
        //fprintf(fp,"#end\n"); 
}

void StrFmt::checkseqcomname() {

	StrFmt *t,*t0,*t1,*t2;
	//check the name if identical
	for(t=this;t;t=t->next)  {
		for(t1=t;t1;t1=t1->more) { 
			if(t1->code==0&&strcmp(t1->token,"structure")==0)  {
				cerr<<"name missing:" <<endl;
				t1->printoutonly(stderr);
				exit(0);
			}
			else if(t1->code==0&&strcmp(t1->token,"composite")==0)  {
				cerr<<"name missing:" <<endl;
				t1->printoutonly(stderr);
				exit(0);
			}
			else if(t1->code==0) {
				continue;
			}
			if(strcmp(t1->token,"sequence")!=0&&strcmp(t1->token,"composite")!=0) continue;
			for(t0=t;t0;t0=t0->next){
				for(t2=t0;t2;t2=t2->more)  {
					if(t2->code==0&&strcmp(t2->token,"structure")==0)  {
						cerr<<"name missing:" <<endl;
						t2->printoutonly(stderr);
						exit(0);
					}
					else if(t2->code==0&&strcmp(t2->token,"composite")==0)  {
						cerr<<"name missing:" <<endl;
						t2->printoutonly(stderr);
						exit(0);
					}
					else if(t2->code==0) {
						continue;
					}
					 
					if(strcmp(t2->token,"sequence")&&strcmp(t2->token,"composite")) continue;
					if(t1==t2) continue;
					if(t1->code&&t2->code&&strcmp(t1->code,t2->code)==0)  {
						cerr<<"sequence code are identical:" <<endl;
						t1->printoutonly(stderr);
						cerr<<"\n\n"<<endl;
						t2->printoutonly(stderr);
						cerr<<"\n\n\ngive different names to the above pir alginmetn!!" <<endl;
						exit(0);
					}
				}
			}
		}
	}
}

void StrFmt::checkseqname() {

	StrFmt *t,*t0,*t1,*t2;

	for(t=this;t;t=t->next)  {
		t1=t->getsequenceStrFmt();
		if(t1==0) continue;
		for(t0=t->next;t0;t0=t0->next) {
			t2=t0->getsequenceStrFmt();
			if(t2==0) continue;
			if(t1->code&&t2->code&&strcmp(t1->code,t2->code)==0) {
				cerr<<"sequence code are identical:" <<endl;
				t->printout(stderr);
				cerr<<"\n\n"<<endl;
				t0->printout(stderr);
				cerr<<"\n\n\ngive different names to the sequence above!!" <<endl;
				exit(0);
			}
		}
	}
}

void StrFmt::checklegal() {

	StrFmt *t,*t0,*t1,*t2;
	StrFmt *c,*c0,*c1,*c2;

	//check token
 	for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {	    
	    if(c0->token==0) {
		cerr<<"\n\n\ntoken missing:\n\n"<<endl;
		c0->printoutonly(stderr);
		exit(0);
	    }	
	    if(strcmp(c0->token,"sequence")!=0&&strcmp(c0->token,"structure")!=0)continue;
	   
	    if(strcmp(c0->token,"structure")==0&&c0->code==0) {
		cerr<<"\n\n\nplease provide the name of the structure:"<<endl;
		c0->printoutonly(stderr);
		exit(0);
	    }
	}

	//check composite zipcode
	for(c=this;c;c=c->next){
       		c0=c->getStrFmt("composite");
		if(c0==0) continue;
		t=c->getStrFmt("sequence");
		if(t==0) continue;
		if(t->zipcode=='\0') {
			cerr<<"you need to specify the id of the sequence pir:"<<endl; 
			t->printoutonly(stderr);
			exit(0);
		}
	}
	for(t=this;t;t=t->next) {
		int i,j;
		

		//one and only one sequence and structure
		i=0;j=0;	 
		for(t0=t;t0;t0=t0->more){
			if(strcmp(t0->token,"sequence")==0) i++;			 
			if(strcmp(t0->token,"structure")==0) j++;
		}
		
		if(i!=1||j>1) {
			cerr<<"the set of alignments should have one for sequence and one for structure"<<endl;
			t->printout(stderr);
			exit(0);	 
		}
 		
		//if composite, then all sequence should exist
		t1=t->getStrFmt("composite");
		for(t0=t;t0&&t1;t0=t0->more) {
			if(strcmp(t0->token,"composite")==0&&t0->seqngap==0) {
				cerr<<"pir content missing"<<endl;
				t0->printoutonly(stderr);
				exit(0);
			}
			else if(strcmp(t0->token,"sequence")==0&&t0->seqngap==0) {
				cerr<<"pir sequence missing!"<<endl;
				t0->printoutonly(stderr);
				exit(0);
			}
			else if(strcmp(t0->token,"structure")==0&&t0->seqngap==0) {
				cerr<<"pir sequence missing!"<<endl;
				t0->printoutonly(stderr);
				exit(0);
			}
			else if(t0->seqngap==0) {
				cerr<<"pir content missing!"<<endl;
				t0->printoutonly(stderr);
				exit(0);
			}
		}

		//check all sequence should be equal
		i=0;
		for(t0=t;t0;t0=t0->more) {
			if(t0->seqngap) i=max(i,strlen(t0->seqngap));
		}
		 
		for(t0=t;t0;t0=t0->more){
			if(t0->seqngap==0) continue;
			int k=strlen(t0->seqngap);
			if(k!=i) {
				cerr<<"the length of alignments in the set should be equal"<<endl;
				t->printout(stderr);
				exit(0);
			} 			
		}
 
		//check if name identical
		for(t0=t;t0;t0=t0->more){
			if(t0->code==0) continue;
			if(strcmp(t0->token,"sequence")!=0&&strcmp(t0->token,"composite")!=0)continue;	
					 
			for(t1=t->next;t1;t1=t1->next) 
			for(t2=t1;t2;t2=t2->more) {
				if(t2==t1) continue;
				if(t2->code==0) continue;
				if(strcmp(t2->token,"sequence")!=0&&strcmp(t2->token,"composite")!=0)continue;			 
				if(strcmp(t0->code,t2->code)==0) {
					cerr<<"sequnce name identical:"<<endl;
					t0->printoutonly(stderr);
					t2->printoutonly(stderr);
				}
			}
		}
	}	

	//check output requirement
	int i=0;
	for(t1=this;t1;t1=t1->next) 
	for(t2=t1;t2;t2=t2->more) {
		if(strcmp(t2->token,"sequence")!=0&&strcmp(t2->token,"composite")!=0)continue;
		if(t2->code) i++; 
	}
	if(i==0) {
		cerr<<endl<<endl;
		cerr<<"there is no output requirement in the file"<<endl;
		cerr<<"the output requirement specified as sequence,sequenceI,sequenceF"<<endl;
		cerr<<"for the token in pir file"<<endl;
		exit(0);
	}

	//check identical sequence id
	
	for(c=this;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
		if(strcmp(c0->token,"sequence")!=0) continue;
		if(c0->zipcode=='\0') continue;
		
		for(c1=c->next;c1;c1=c1->next)
        	for(c2=c1;c2;c2=c2->more) {
			if(strcmp(c2->token,"sequence")!=0) continue;
			if(c2->zipcode=='\0') continue;
			if(c2->zipcode==c0->zipcode) {
				cerr<<"\n\n\ntwo pir have identical id..\n"<<endl;
				c0->printout(stderr);
				c1->printout(stderr);
				exit(0);
			}
		}
	}
}
void StrFmt::checklegalseqxyz() {

	StrFmt *t,*t0;

	for(t=this;t;t=t->next) {
		int j;
		j=0;
		for(t0=t;t0;t0=t0->more){			
			if(strcmp(t0->token,"sequence")==0) j++;
			else {
				cerr<<"the token:"<<t0->token<<" is not allowed in seqxyz"<<endl;
				exit(0);
			}
		}
		if(j==0) {
			cerr<<"the set of alignments should have at least one sequence pir"<<endl;
			t->printout(stderr);
			exit(0);	 
		}
	}	
	
}
void StrFmt::checklegalsuperimpose() {

	StrFmt *t,*t0;

  	
	for(t=this;t;t=t->next) {
		int i,j;
		i=0;j=0;
		for(t0=t;t0;t0=t0->more){			
			i++;
		}
		if(i>2) {
			cerr<<"the set of alignments has more than 2 sequences"<<endl;
			t->printout(stderr);
			exit(0);
		}
		else if(i<2) {
			cerr<<"the set of alignments has less than 2 sequences"<<endl;
			t->printout(stderr);
			exit(0);
		}
		i=0;j=0;
		for(t0=t;t0;t0=t0->more){			
			if(strcmp(t0->token,"structure")==0) j++;
		}
		if(j!=2) {
			cerr<<"the set of alignments should have two for structure"<<endl;
			t->printout(stderr);
			exit(0);	 
		}

		i=strlen(t->seqngap);
		j=strlen(t->more->seqngap);
		if(i!=j) {
			cerr<<"the length of alignments in the set should be equal"<<endl;
			t->printout(stderr);
			exit(0);
		}
	}	
	
}
void StrFmt::setseqnresn() {

	int n=strlen(seqn);

        //remove residue does not exist in structure
        int i,j;
        for(i=0;i<n;i++) {
                if(resn[i]==0) {
			j=compare[i];
                        seqngap[j]='-';
                }
        }
	int m=0;
	for(i=0;i<n;i++) {
		if(resn[i]==0) continue;
		resn[m++]=resn[i];
	}
	resn[m]=0;
	for(i=m;i<n;i++) resn[i]=0;

	m=0;
	n=strlen(seqngap);
	for(i=0;i<n;i++) {
		if(seqngap[i]=='-') continue;
		seqn[m++]=seqngap[i];
	}
	seqn[m]='\0';
}

float *StrFmt::gethbondbounds(Atm *a,Atm *a0){

	Res *r=a->res;
	Res *r0=a0->res;

    	int n=compare[r->id0];
    	int n0=compare[r0->id0];

    	if(seqngap[n]!=r->name||seqngap[n0]!=r0->name) {
        	cerr<<"Warning! the residue id and the name does not match";
    	}

    	StrFmt *t,*t0;

    	//float *dst=0;

	float tt[1000];
	StrFmt *mtt[1000];

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

float *StrFmt::getfardistbounds(Atm *a,Atm *a0) {

	//float dp[1000];

	//int ndp=0;

	Res *r=a->res;
    	Res *r0=a0->res;

    	int n=compare[r->id0];
    	int n0=compare[r0->id0];

    	if(seqngap[n]!=r->name||seqngap[n0]!=r0->name) {
        	cerr<<"Warning! the residue id and the name does not match";
    	}

    	StrFmt *t,*t0;

    	//float *dst=0;

	float tt[1000];
	StrFmt *mtt[1000];

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


float *StrFmt::getdistbounds(Atm *a,Atm *a0) {


	Res *r=a->res;
    	Res *r0=a0->res;

    	StrFmt *t,*t0;


	float tt[1000];
	StrFmt *mtt[1000];

	//get all constraints

	int ntt=0;
    	for(t=this;t;t=t->next) {
		StrFmt *s=t->findsequencefmt();
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

void StrFmt::setsitescore() {

	StrFmt *t,*t0;
        for(t0=this;t0;t0=t0->next) {

                for(t=t0;t;t=t->more) {
			if(strcmp(t->token,"structure")==0) continue;
			
			if(t->sitescore) delete [] t->sitescore;
			int n=strlen(t->seqngap);
			t->sitescore=new float[n+1];
			int i=0;
			for(i=0;i<n+1;i++) t->sitescore[i]=0;		
					
			for(i=0;i<n;i++) t->sitescore[i]=t->getalgnweight(i);
		}
	}
}



void StrFmt::setnewsitescore() {

	StrFmt *root=getrootStrFmt();
	if(root->onlyrefine==0) {
	cerr<<endl<<"calculate sequence conserve score at each residue between sequence and structure pirs..."<<endl<<endl;
	}
	StrFmt *t;
        for(t=this;t;t=t->next) {
		//if(strcmp(t->token,"structure")!=0) continue;
		if(t->sitescore) delete [] t->sitescore;
		int n=strlen(t->seqngap);
		t->sitescore=new float[n+1];
		StrFmt *se=t->getStrFmt("sequence");
		StrFmt *st=t->getStrFmt("structure");
		
		int i=0;
		for(i=0;i<n+1;i++) t->sitescore[i]=0;	
		if(se==0||st==0) continue;
		
		for(i=0;i<n;i++)   t->sitescore[i]=t->getsimplealgnweight(i);				
		for(i=0;i<n;i++) {
			if(root->onlyrefine==0) {
			cerr<<"conserve score  "<< i<<" "<<st->seqngap[i]<<"---"<<se->seqngap[i]<<":"<<t->sitescore[i]<<endl;
			}
		}
	}

	avgscore=calcavgscore();
}

float StrFmt::calcavgscore(){

	int n=strlen(seqngap);

	int m=0;
	float a=0;
	for(int i=0;i<n;i++) {
		if(sitescore[i]==0) continue;
		a+=sitescore[i];
		m++;
	}
	if(m) return a/m;
	else  return 0;
}


float *StrFmt::getdistbounds(Atm *a,Atm *a0,float low,float high) {

	 
	Res *r=a->res;
    	Res *r0=a0->res;
 
    	StrFmt *t,*t0;
 
	float tt[200];
	StrFmt *mtt[200];

	//get all constraints

	int ntt=0;
	int men=0;
    	for(t=this;t;t=t->next) {
		StrFmt *s=t->findsequencefmt();
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

float StrFmt::getalgnweight(int n)  {

	int len=strlen(seqngap);

        StrFmt *t=getparentStrFmt();
        StrFmt *s0=t->findsequencefmt();

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
                        //if(n1<n2) a+=blosum65mt[n1*(n1+1)/2+n2]/x;
                        //else      a+=blosum65mt[n2*(n2+1)/2+n1]/x;
			if(n1>n2) a+=blosum65mt[n1*(n1+1)/2+n2]/x;
			else      a+=blosum65mt[n2*(n2+1)/2+n1]/x;
                }
        }

	if(blosum65mt) delete [] blosum65mt;
        if(am*bm==0) return 0;

        float c=a/sqrt(am*bm);

        return c;
}

float StrFmt::findmatrixminvalue() {
 
	Algn algn;

	float *blosum65mt=algn.getscoretable("blosum65mt");

	if(blosum65mt==0) return 0;

	int n=strlen(algn.aaorder);

	float e=10000;

	for(int i=0;i<n*(n-1)/2;i++) {
		if(blosum65mt[i]<e) e=blosum65mt[i];	 			
	}

	if(blosum65mt) delete [] blosum65mt;

	return e;
}


float StrFmt::getnewalgnweight(int n)  {

	int len=strlen(seqngap);

	Algn algn;

	float *blosum65mt=algn.getscoretable(matrix);
 
	if(blosum65mt==0) return 0;
	
        float a=0;
        float am=0;
        float bm=0;
        for(int i=0;i<len;i++)  {
                //if(seqngap[i]=='-'&&s0->seqngap[i]=='-') continue;
		if(fabs(i-n)>5) continue;
		for(StrFmt *t=this;t;t=t->more) 
		for(StrFmt *t0=t->more;t0;t0=t0->more) {
					
                	int n1= TRES.getresid(t->seqngap[i]);
                	int n2= TRES.getresid(t0->seqngap[i]);
                	if(n1==-1&&n2==-1) continue;
			float x=1./(1+abs(i-n));
			x=pow(x,0.25);
                	if(n1==-1) {
                        	//bm+=blosum65mt[n2*(n2+1)/2+n2]*x;
				a+=gap*x;
                	}
                	else if(n2==-1) {
                        	//am+=blosum65mt[n1*(n1+1)/2+n1]*x;
				a+=gap*x;
                	}
                	else {
                        	//am+=blosum65mt[n1*(n1+1)/2+n1]*x;
                        	//bm+=blosum65mt[n2*(n2+1)/2+n2]*x;
                        	//if(n1<n2) a+=blosum65mt[n1*(n1+1)/2+n2]*x;
                        	//else      a+=blosum65mt[n2*(n2+1)/2+n1]*x;
				if(n1>n2) a+=blosum65mt[n1*(n1+1)/2+n2]*x;
                        	else      a+=blosum65mt[n2*(n2+1)/2+n1]*x;
                	}

		}
        }
	am=bm=1;
	if(blosum65mt) delete [] blosum65mt;
        if(am*bm==0) return 0;

        float c=a/sqrt(am*bm);
	if(TRES.logg) cerr<<n<<" "<<c<<" "<<seqngap[n]<<" "<<endl;
        return c;
}

float StrFmt::getsimplealgnweight(int n)  {

	Algn algn;

        float *blosum65mt=algn.getscoretable(matrix);

        if(blosum65mt==0) return 0;

	int len=strlen(seqngap);
 
        float a=0;
        float am=0;
	int n1,n2; 
	float b,c;
	
	 
        //for(int i=-5;i<len+5;i++)  {               
	for(int i=n-5;i<=n+5;i++) {
		if(fabs(i-n)>=5) continue;

		if(i<0||i>=len) {
			float x=1./(1+abs(i-n));
			x=pow(x,0.75);
			am+=2*x;
			continue;
		}
		StrFmt *t,*t0;
		for(t=this;t;t=t->more) 
		for(t0=t->more;t0;t0=t0->more) {	
			if(strcmp(t->token,"structure")!=0&&strcmp(t->token,"sequence")!=0) continue;
			if(strcmp(t0->token,"structure")!=0&&strcmp(t0->token,"sequence")!=0) continue;
			if(t->seqngap[i]=='-'&&t0->seqngap[i]=='-') continue;

			float x=1./(1+abs(i-n));
			x=pow(x,0.75);
                	if(t->seqngap[i]=='-') {                        	
				am+=2*x;
                	}
                	else if(t0->seqngap[i]=='-') {                        	
				am+=2*x;
                	}
                	else {                        	
				if(t->seqngap[i]==t0->seqngap[i]) a+=2*x;
				else {
					n1= TRES.getresid(t->seqngap[i]);
                        		n2= TRES.getresid(t0->seqngap[i]);
					if(n1!=-1&&n2!=-1) {
						//if(n1<n2) b=blosum65mt[n1*(n1+1)/2+n2];
                                		//else      b=blosum65mt[n2*(n2+1)/2+n1];
						if(n1>n2) b=blosum65mt[n1*(n1+1)/2+n2];
                                		else      b=blosum65mt[n2*(n2+1)/2+n1];
						c=blosum65mt[n1*(n1+1)/2+n1];
						c=c*blosum65mt[n2*(n2+1)/2+n2];
						c=sqrt(c);
						b=b/c;
						if(b<0) b=0;
					}
					else {
						b=0;
					}						
					a+=(b+1)*x;
				}
				am+=2*x;
                	}
		}
        }


	if(blosum65mt) delete [] blosum65mt;
 	if(am==0) return 0;
        c=a/am;
	if(TRES.logg) cerr<<n<<" "<<c<<" "<<seqngap[n]<<" "<<endl;
        return c;
}


float StrFmt::getsimplealgnbound(int n)  {

	Algn algn;

        float *blosum65mt=algn.getscoretable(matrix);

        if(blosum65mt==0) return 0;

	int len=strlen(seqngap);
 
        float a=0;
        float am=0;
	int n1,n2; 
	float b,c;
	
	                   
	for(int i=n-5;i<=n+5;i++) {
		if(fabs(i-n)>=5) continue;

		if(i<0||i>=len) {
			float x=1./(1+abs(i-n));
			x=pow(x,0.75);
			am+=2*x;
			continue;
		}
		for(StrFmt *t=this;t;t=t->more) 
		for(StrFmt *t0=t->more;t0;t0=t0->more) {	

			if(t->seqngap[i]=='-'&&t0->seqngap[i]=='-') continue;

			float x=1./(1+abs(i-n));
			x=pow(x,0.75);
                	if(t->seqngap[i]=='-') {                        	
				am+=2*x;
                	}
                	else if(t0->seqngap[i]=='-') {                        	
				am+=2*x;
                	}
                	else {                        	
				if(t->seqngap[i]==t0->seqngap[i]) a+=2*x;
				else {
					n1= TRES.getresid(t->seqngap[i]);
                        		n2= TRES.getresid(t0->seqngap[i]);
					if(n1!=-1&&n2!=-1) {
						//if(n1<n2) b=blosum65mt[n1*(n1+1)/2+n2];
                                		//else      b=blosum65mt[n2*(n2+1)/2+n1];
						if(n1>n2) b=blosum65mt[n1*(n1+1)/2+n2];
                                		else      b=blosum65mt[n2*(n2+1)/2+n1];
						c=blosum65mt[n1*(n1+1)/2+n1];
						c=c*blosum65mt[n2*(n2+1)/2+n2];
						c=sqrt(c);
						b=b/c;
						if(b<0) b=0;
					}
					else {
						b=0;
					}						
					a+=(b+1)*x;
				}
				am+=2*x;
                	}
		}
        }


	if(blosum65mt) delete [] blosum65mt;
 	if(am==0) return 0;
        c=a/am;
	if(TRES.logg) cerr<<n<<" "<<c<<" "<<seqngap[n]<<" "<<endl;
        return c;
}



float StrFmt::getalgnscore() {

	int len=strlen(seqngap);

	StrFmt *t=getparentStrFmt();
	StrFmt *s0=t->findsequencefmt();

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
			//if(n1<n2) a+=blosum65mt[n1*(n1+1)/2+n2];
			//else 	  a+=blosum65mt[n2*(n2+1)/2+n1];
			if(n1>n2) a+=blosum65mt[n1*(n1+1)/2+n2];
			else 	  a+=blosum65mt[n2*(n2+1)/2+n1];
		}
	}

	if(blosum65mt) delete [] blosum65mt;

	if(am*bm==0) return 0;

	float c=a/sqrt(am*bm);

	return c;
}


void StrFmt::setconservedres() {

	//prepare score
	float *score=0;

	if(score) delete [] score; score=0;
	int n=strlen(seqngap);
	score=new float[n];

	int i=0;
	for(i=0;i<n;i++) { score[i]=0; }

	//calculate score of each alignment

	StrFmt *m,*s0;


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

			//if(ci<c0) score[i]=blosum65mt[ci*(ci+1)/2+c0];
			//else      score[i]=blosum65mt[c0*(c0+1)/2+ci];
			if(ci>c0) score[i]=blosum65mt[ci*(ci+1)/2+c0];
			else      score[i]=blosum65mt[c0*(c0+1)/2+ci];
		}
	}
	if(blosum65mt)delete [] blosum65mt;
}

void StrFmt::setallconservedres() {

	StrFmt *t;
	for(t=this;t;t=t->next) {
		t->setconservedres();
	}
}

int StrFmt::isconservedpair(int n,int n0) {

	/*
	n=compare[n];
	n0=compare[n0];

	StrFmt *t,*t0;
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

int StrFmt::getcreditofhbond(Atm *a,Atm *a0) {
	
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

float *StrFmt::gethbondbounds(Res *r,Res *r0) {

	//take all hydrogens out including s-s bond and ionic bond

        StrFmt *mtt[1000];
	Atm	   *aaa[1000];
	float       aln[500];
	int         credit[500];

	StrFmt *t,*t0;

	int ntt=0;

        for(t=this;t;t=t->next) {
                StrFmt *s=t->findsequencefmt();
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
		//int n1,n2; 
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

int StrFmt::ishbondexist(Res *r, Res *r0) {
	
        StrFmt *t,*t0;
        for(t=this;t;t=t->next) {
		StrFmt *s=t->findsequencefmt();
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

void StrFmt::setsecstruct() {

	StrFmt *t,*t0;
        for(t=this;t;t=t->next) {
		for(t0=t;t0;t0=t0->more) {
			if(strcmp(t0->token,"structure")!=0) continue;
			t0->pdb->chn->clearhbond();
			t0->pdb->chn->header();
			t0->pdb->chn->buildhbond();
			t0->pdb->chn->setdsspstr();
			t0->pdb->chn->buildssbond();
			t0->pdb->chn->setthreestatesec();
			cerr<<"detecting secondary structure regions using dssp definition..."<<endl<<endl; 
			t0->pdb->chn->writesecondary(stderr);
			if(TRES.logg>3) cerr<<t0->pdb->chn->get2d()<<endl;
		}
	}
}

int  StrFmt::getstructurenumber(){
	
	StrFmt *t,*t0;

	int n=0;
        for(t0=this;t0;t0=t0->next) 
        for(t=t0;t;t=t->more) {
		if(strcmp(t->token,"structure")) continue;
		n++;
	}	
	return n;
}

void StrFmt::setnest() {

	StrFmt *t,*t0,*se;

 
        for(t0=this;t0;t0=t0->next) 
        for(t=t0;t;t=t->more) {		
		if(strcmp(t->token,"structure")) continue;
		se=t0->getsequenceStrFmt();
		t0->mutate=new Mutate();
		if(se->code) t0->mutate->code=strdup(se->code);
		t0->mutate->mpdb=new Pdb(t->pdb);
		t0->mutate->sqnto=strdup(t->seqngap); 
		t0->mutate->owner=t;
		t0->mutate->setseqnres();
		t0->mutate->mpdb->configure();
		//t->mutate->mpdb->write("start.pdb");
		t0->mutate->mpdb->chn->header();
		int n=strlen(t->seqngap);
		t0->mutate->match=new int[n+10];
		t0->mutate->compare=new int[n+10];
		t0->mutate->resn=new Res*[n+10];
		t0->mutate->seqn=new char[n+10];
		t0->mutate->resid=new int[n+10];
		for(int i=0;i<n;i++) {
			t0->mutate->match[i]=t->match[i];
			t0->mutate->compare[i]=t->compare[i];
			t0->mutate->resn[i]=t->resn[i];
			t0->mutate->seqn[i]=t->seqn[i];
			t0->mutate->resid[i]=-1;
		}
		t0->mutate->setmatch();
		t0->mutate->setdssp();
		t0->mutate->setboundscore();
	}	
}

void StrFmt::setnest0() {

	StrFmt *t,*t0;

 
        for(t0=this;t0;t0=t0->next) 
        for(t=t0;t;t=t->more) {
		if(strcmp(t->token,"structure")) continue;
		t->mutate=new Mutate();
		t->mutate->mpdb=new Pdb(pdb);
		t->mutate->sqnto=strdup(t->seqngap); 
		t->mutate->owner=t;
		t->mutate->setseqnres();
		t->mutate->mpdb->configure();
		//t->mutate->mpdb->write("start.pdb");
		t->mutate->mpdb->chn->header();
		int n=strlen(t->seqngap);
		t->mutate->match=new int[n+10];
		t->mutate->compare=new int[n+10];
		t->mutate->resn=new Res*[n+10];
		for(int i=0;i<n;i++) {
			t->mutate->match[i]=t->match[i];
			t->mutate->compare[i]=t->compare[i];
			t->mutate->resn[i]=t->resn[i];
		}
		t->mutate->setmatch();
		t->mutate->setdssp();
		t->mutate->setboundscore();
	}	
}

void StrFmt::setmutatedone() {

	StrFmt *t,*t0;
 
        for(t0=this;t0;t0=t0->next) 
        for(t=t0;t;t=t->more) {
		if(strcmp(t->token,"structure")) continue;
		t->mutate->domutate();		 
	}	
}


void StrFmt::setalldssp() {

	StrFmt *t,*t0;
	for(t0=this;t0;t0=t0->next)
        for(t=t0;t;t=t->more) {
		if(strcmp(t->token,"sequence")==0) t->setdssp();
	}
}

void StrFmt::setalloldgap() {

        StrFmt *t,*t0;
        for(t0=this;t0;t0=t0->next)
        for(t=t0;t;t=t->more) {
		if(strcmp(t->token,"sequence")==0) {
			if(t->oldgap) delete [] t->oldgap;
			t->oldgap=strdup(t->seqngap);
		}
		else if(strcmp(t->token,"structure")==0) {
			if(t->oldgap) delete [] t->oldgap;
			t->oldgap=strdup(t->seqngap);
		}
		else if(strcmp(t->token,"composite")==0) {
			if(t->oldgap) delete [] t->oldgap;
			t->oldgap=strdup(t->seqngap);
		}
		else if(strcmp(t->token,"secondary")==0) {
			if(t->oldgap) delete [] t->oldgap;
			t->oldgap=strdup(t->seqngap);
		}
		else {
			if(t->oldgap) delete [] t->oldgap;
			t->oldgap=strdup(t->seqngap);
		}
        }
}


void StrFmt::tuneallalign() {

	StrFmt *t,*t0;
	if(tune!=0) {
		cerr<<endl<<"tuning sequence alignment..."<<endl<<endl;
	}	
	for(t=this;t;t=t->next)  {
		for(t0=t;t0;t0=t0->more) {
			if(strcmp(t0->token,"sequence")!=0) continue;
			if(tune>=1) t0->tunealign();			
			if(tune>=2) t0->smoothinggap();
			if(tune>=3) t0->tuneterminal();
		}
	}
}

int StrFmt::tuneterminal(){

	cerr<<endl<<"remove gaps at terminal..."<<endl<<endl;

	StrFmt *parent=getparentStrFmt();
        StrFmt *se = parent->findstructurefmt();
        if(se==0) return 0;
	se->setmatch();
	setmatch();
	
	int nlen=strlen(seqngap);
	int nte=strlen(seqn);
	char *dssp=setnewdssp();
	dssp=setstate(dssp);
	int i,n,j;
	//float c;
	int numc=0;
	while(1) {
		n=0;
		int mlen=0;
		int sheet=0;
		int alpha=0;
		for(i=0;i<nlen;i++) {
			if(se->seqngap[i]=='-'&&seqngap[i]=='-') continue;
			if(se->seqngap[i]!='-'&&seqngap[i]!='-') {
				mlen++;
				if(dssp[i]=='e') sheet++;
				if(dssp[i]=='h') alpha++;				
				continue;
			}
			if(mlen)  break;		
		}
		numc++;
		if(numc>50) break;
		if(nte==mlen) break;
		if(i==nlen) continue;		 
		if(mlen<5||(mlen<10&&alpha<3&&sheet<3)) {
			if(se->seqngap[i]!='-'&&seqngap[i]=='-') {//delete
						
				for(j=i;j>0;j--) {
					if(seqngap[j-1]=='-') continue; 
					char c; 
					c=seqngap[j];
					seqngap[j]=seqngap[j-1];
					seqngap[j-1]=c;
					n=1; 
				}									 
			}
			else if(se->seqngap[i]=='-'&&seqngap[i]!='-') {
				
				for(j=i;j>0;j--) {
					if(se->seqngap[j-1]=='-') continue;   
					char c; 
					c=se->seqngap[j];
					se->seqngap[j]=se->seqngap[j-1];
					se->seqngap[j-1]=c;
					n=2;
				}									 
			}	
		}

		if(n==0) break; 
		if(dssp) delete [] dssp;dssp=0;
		se->setmatch();
		setmatch();		
 		dssp=setnewdssp(); 					 
		dssp=setstate(dssp);
	}
	
	//end terminal
	//int nte=strlen(seqn);
	while(1) {
		n=0;
		int mlen=0;
		int sheet=0;
		int alpha=0;
		for(i=nlen-1;i>=0;i--) {
			if(se->seqngap[i]=='-'&&seqngap[i]=='-') continue;
			if(se->seqngap[i]!='-'&&seqngap[i]!='-') {
				mlen++;
				if(dssp[i]=='e') sheet++;
				if(dssp[i]=='h') alpha++;				
				continue;
			}
			if(mlen)  break;		
		}
		if(nte==mlen) break;
		if(i==-1) continue;		 
		if(mlen<5||(mlen<10&&alpha<3&&sheet<3)) {
			if(se->seqngap[i]!='-'&&seqngap[i]=='-') {//delete
						
				for(j=i;j<nlen-1;j++) {
					if(seqngap[j+1]=='-') continue; 
					char c;		
					c=seqngap[j];
					seqngap[j]=seqngap[j+1];
					seqngap[j+1]=c;
					n=1;		
				}									 
			}
			else if(se->seqngap[i]=='-'&&seqngap[i]!='-') {
				
				for(j=i;j<nlen-1;j++) {
					if(se->seqngap[j+1]=='-') continue;
					char c;					 
					c=se->seqngap[j];
					se->seqngap[j]=se->seqngap[j+1];
					se->seqngap[j+1]=c;
					n=2;
				}									 
			}	
		}

		if(n==0) break; 
		if(dssp) delete [] dssp;dssp=0;
		se->setmatch();
		setmatch();		
 		dssp=setnewdssp(); 					 
		dssp=setstate(dssp);
	}
	if(dssp) delete [] dssp;dssp=0;	
	dssp=setnewdssp(); 
	dssp=setstate(dssp);
	if(TRES.logg) cerr<<dssp<<endl;
	if(TRES.logg) cerr<<se->seqngap<<endl;	
	if(TRES.logg) cerr<<seqngap<<endl;
	if(TRES.logg) cerr<<"end...."<<endl<<endl;		
	if(dssp) delete [] dssp;dssp=0;
	return 0;
} 

int StrFmt::tunealign() {
 
	cerr<<endl<<"remove gaps in secondary structure..."<<endl<<endl;

	StrFmt *parent=getparentStrFmt();
        StrFmt *se = parent->findstructurefmt();
        if(se==0) return 0;
	se->setmatch();
	setmatch();
	int nlen=strlen(seqngap);
	char *dssp=setnewdssp();
	dssp=setstate(dssp);
	int n,i,j;
	char c;
	n=2;
	int nout=0;
	int nzero=0;
	int numc=0;
	while(1) { 
		numc++;
		if(numc>50) break;
		n=0;				
		for(i=0;i<nlen;i++) {			
			
			if(dssp[i]!='h'&&dssp[i]!='e') continue; //only in the secondary region
			if(se->seqngap[i]=='-'&&seqngap[i]=='-') continue;
			if(se->seqngap[i]!='-'&&seqngap[i]!='-') continue;	
			if(i==0||i==nlen-1) continue;
			
			int n1,n2;		
			
			for(n1=i-1;n1>=0;n1--) {
				if(dssp[n1]!=dssp[i]) break;
			}

			for(n2=i+1;n2<nlen;n2++) {
				if(dssp[n2]!=dssp[i]) break;
			}

			//take care of terminal case
			
			if(nzero&&i==n1+1) {
				if(seqngap[i]!='-') continue; 
				if(dssp[n1]=='e'||dssp[n1]=='h')  continue;
				int n3,p1,p2,m; 
				for(n3=n1;n3>=0;n3--) {
					if(dssp[n3]=='e'||dssp[n3]=='h') break;
				}
				n3++;
				n3=max(n3,n1-10);
				if(n3<0) n3=0;
				p1=0;p2=0;
				for(m=n1;m>=n3;m--) {
					if(se->seqngap[m]!='-') p1++;
					if(seqngap[m]!='-') p2++;
					if(p2>=p1&&p2>0) break;
					//if(p2>0) break;
				}
				if(p2<p1) continue;
				if(p2==0) continue;
				//if(p2==p1&&p1==0) continue;
				for(m=n1;m>=n3;m--) {					
					if(seqngap[m]!='-') break;
				}
				c=seqngap[i];
				seqngap[i]=seqngap[m];
				seqngap[m]=c;
				n=1;
				break;
			}
			else if(nzero&&i==n2-1) {
				if(seqngap[i]!='-') continue; 
				if(dssp[n2]=='e'||dssp[n2]=='h')  continue;
				int n3,p1,p2,m;
				for(n3=n2;n3<nlen;n3++) {
					if(dssp[n3]=='e'||dssp[n3]=='h') break;
				}
				n3--;
				n3=min(n3,n2+10);
				if(n3>=nlen) n3=nlen-1;
				p1=0;p2=0;
				for(m=n2;m<=n3;m++) {
					if(se->seqngap[m]!='-') p1++;
					if(seqngap[m]!='-') p2++;
					if(p2>=p1&&p2>0) break;
					//if(p2>0) break;
				}
				if(p2<p1) continue;
				if(p2==0) continue;
				//if(p2==p1&&p1==0) continue;
				for(m=n2;m<=n3;m++) {					
					if(seqngap[m]!='-') break;
				}
				c=seqngap[i];
				seqngap[i]=seqngap[m];
				seqngap[m]=c;
				n=1;
				break;
			}

			if(n1==-1) n1=0;
			if(n2==nlen) n2=nlen-1;
 	
			float a=parent->getsimplealgnweight(n1);
			float b=parent->getsimplealgnweight(n2);

			int m=0;
			
			
			if(i==0||i==nlen-1) continue;
			if(dssp[i]!='-'&&dssp[i-1]=='-') {
				m=-1;
			}
			else if(dssp[i]!='-'&&dssp[i+1]=='-') {
				m=1;
			}
			else if(i-n1>2*(n2-i)) {
				m=1;
			}
			else if(2*(i-n1)<n2-i) {
				m=-1;
			}
			else if(a>b) {
				m=1;
			}
			else if(a<b) {
				m=-1;
			}
			else  if(i-n1>n2-i){
				m=1;
			}
			else if(i-n1<n2-i) {
				m=-1;
			}
			else m=1;
									

			if(se->seqngap[i]!='-'&&seqngap[i]=='-') {//delete
				if(m==1) {						 
					for(j=i;j<n2-1;j++) {
						if(seqngap[j+1]=='-') continue; 		
						c=seqngap[j];
						seqngap[j]=seqngap[j+1];
						seqngap[j+1]=c;
						n=1;						 
					}										 
				}
				else {			
					for(j=i;j>n1+1;j--) {
						if(seqngap[j-1]=='-') continue;  
						c=seqngap[j];
						seqngap[j]=seqngap[j-1];
						seqngap[j-1]=c;
						n=1;
					}					
				}
			}
			else if(se->seqngap[i]=='-'&&seqngap[i]!='-') {
				if(m==1) {
					
					for(j=i;j<n2-1;j++) {	
						if(se->seqngap[j+1]=='-') continue; 				 
						c=se->seqngap[j];
						se->seqngap[j]=se->seqngap[j+1];
						se->seqngap[j+1]=c;
						n=2;
					}					
				}
				else {
					
					for(j=i;j>n1+1;j--) {
						if(se->seqngap[j-1]=='-') continue;   
						c=se->seqngap[j];
						se->seqngap[j]=se->seqngap[j-1];
						se->seqngap[j-1]=c;
						n=2;
					}					
				}
			}
			if(n) nout++;
			if(n) break;
		}
 
		if(dssp) delete [] dssp;dssp=0;
		se->setmatch();
		setmatch();		
 		dssp=setnewdssp(); 
		dssp=setstate(dssp);		
		if(n==0&&nzero) break;
		else if(n) nzero=0;
		else nzero++;
	}	
	if(dssp) delete [] dssp;dssp=0;	
	dssp=setnewdssp(); 
	dssp=setstate(dssp);
	if(TRES.logg) cerr<<dssp<<endl;
	if(TRES.logg) cerr<<se->seqngap<<endl;	
	if(TRES.logg) cerr<<seqngap<<endl;
	if(TRES.logg) cerr<<"end...."<<endl<<endl;	
	delete [] dssp;dssp=0;
	return nout;
}

char *StrFmt::setstate(char *dssp) {

	if(dssp==0) return dssp;

	int n=strlen(seqngap);

	int i;

	for(i=0;i<n;i++) {
		if(dssp[i]!='e'&&dssp[i]!='h') dssp[i]='-';
	}

	return dssp;
}

void StrFmt::setotherfmt0(int n) {

	StrFmt *parent=getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
	
        if(se==0) return;

	if(other) delete other; other=0;
	
	if(mutate==0) return;

	if(mutate->mpdb==0) return;

	other=new StrFmt;

	Pdb *spdb =mutate->mpdb;
	 
	Res *r0;

	StrFmt *str=mutate->owner;

 
	for(r0=spdb->chn->res;r0;r0=r0->next) {

		if(r0->sec=='-') continue;

 		Res *r;
		for(r=r0;r;r=r->next) {
			if(r->sec=='-') break;
		}	
		if(r) r=r->last;
		
		Res *t0, *t;

		for(t0=r0->last;t0;t0=t0->last) {
			if(t0->sec!='-') break;
		}			
		if(t0==0) t0=spdb->chn->res;
		else t0=t0->next;

		for(t=r->next;t;t=t->next) {
			if(t->sec!='-') break;
		}
			
		if(t==0) t=spdb->chn->lastres();
		else t=t->last;

		int n1=0;
		int n2=0;

		Res *f;
		
		int low1,low2;
		float d=1000;
		low1=mutate->compare[t0->id0];
		d=parent->sitescore[low1];
		for(f=t0;f;f=f->next) {		
			if(f->id0>=r0->id0) break;
			int i=mutate->compare[f->id0];
			if(se->seqngap[i]=='-'&&str->seqngap[i]=='-') continue;
			if(parent->sitescore[i]<d) {
				d=parent->sitescore[i];
				low1=i;	
			}
			if(se->seqngap[i]!='-')  n1++;
			if(str->seqngap[i]!='-') n2++;
		}	
		
		

		n1=0;
		n2=0;
		
		low2=mutate->compare[r->id0];
		d=parent->sitescore[low2];
		for(f=r;f;f=f->next) {		
			if(f->id0>=t->id0) break;
			int i=mutate->compare[f->id0];
			if(se->seqngap[i]=='-'&&str->seqngap[i]=='-') continue;
			if(parent->sitescore[i]<d) {
				d=parent->sitescore[i];
				low2=i;	
			}
			if(se->seqngap[i]!='-')  n1++;
			if(str->seqngap[i]!='-') n2++;
		}	

	
 
	}	
}

void StrFmt::setotherfmt(int ner,float dcut,int shift) {

	StrFmt *parent=getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
	
        if(se==0) return;
	
	//if(other) delete other; other=0;
	
	if(mutate==0) return;

	if(mutate->mpdb==0) return;
	
	StrFmt *sy,*sy0;

	sy=0;
	sy0=0;

	for(sy=other;sy;sy=sy->next) {
		sy0=sy;
	}

	if(sy0==0)  {
		other=new StrFmt;
		sy=other;
	}
	else {
		sy0->next=new StrFmt;
		sy=sy0->next; 
	} 
	    
	Pdb *spdb=spdb=mutate->mpdb;
	mutate->sethbonddssp(4);  
	mutate->setdssp(1); 
	Res *r0;

	StrFmt *str=mutate->owner;

	StrFmt *s;

	int nlen=strlen(seqngap);

	s=sy;
	Algn algn;
	float *blosum65mt=algn.getscoretable(matrix);
	algn.opncost=3*blosum65mt[9];
	algn.gapcost=0.2;
	Strhandler cc;
	spdb->setlastresorder();
	mutate->owner->pdb->setlastresorder();
	for(r0=spdb->chn->res;r0;r0=r0->next) {

		if(r0->sec=='-') continue;

 		Res *r;
		for(r=r0;r;r=r->next) {
			if(r->sec=='-') break;			 
		}	
		if(r) r=r->last;
		
		Res *t0, *t;

		for(t0=r0->last;t0;t0=t0->last) {
			if(t0->sec!='-') break;
			if(r0->id0-t0->id0>shift+5) break;
		}			
		if(t0==0) t0=spdb->chn->res;
		else t0=t0->next;

		for(t=r->next;t;t=t->next) {
			if(t->sec!='-') break;
			if(t->id0-r->id0>shift+5) break;
		}
			
		if(t==0) t=spdb->chn->lastres();
		else t=t->last;
 
		//test
		Res *of1=0;
		Res *of2=0;
		//find the starting fitting residues
		for(of1=t0;of1;of1=of1->last) { 		
			int h=mutate->compare[of1->id0];
			//if(ner==1)  h=mutate->owner->match[h];
			//else	    h=mutate->match[h];
			h=mutate->owner->match[h];
			if(h==-1) continue;
			Res *fh1=0;
			//if(ner==1) {//owner
			fh1=mutate->owner->resn[h];	
			//}
			//else {
			//	fh1=mutate->resn[h];
			//}
			if(fh1==0) continue;
			float dis=of1->directrmsdanyway(fh1,0,3);
			if(dis<0.2) break; 
		}			
		
		//find the ending fitting residues
		for(of2=t;of2;of2=of2->next) { 		
			int h=mutate->compare[of2->id0];
			//if(ner==1)  h=mutate->owner->match[h];	
			//else	    h=mutate->match[h];
			h=mutate->owner->match[h];
			if(h==-1) continue;
			Res *fh1=0;
			//if(ner==1) {//owner
			fh1=mutate->owner->resn[h];	
			//}
			//else { //self
			//	fh1=mutate->resn[h];
			//}
			if(fh1==0) continue;
			float dis=of2->directrmsdanyway(fh1,0,3);
			if(dis<0.2) break; 
		}			
		if(of2==0||of1==0) {
			if(ner==1) {r0=r;continue;}
		}
		else if(of1&&of2) {
			if(ner==2) {r0=r;continue;}
		}
		//end
		

		char *ses,*strs;

		int n=t->id0-t0->id0+100;

		ses=new char[n];
		strs=new char[n];
 		
		//find out end residues id
		int n1=mutate->compare[t0->id0];
		int n2=mutate->compare[t->id0];
	
		int g1=mutate->compare[r0->id0];
		int g2=mutate->compare[r->id0];
		int stem[2];
		stem[0]=n1;
		stem[1]=n2;
		if(TRES.logg)mutate->printsegment(stem);
		//take out sequence
		int i; 
		int m1=0;
		int m2=0;
		for(i=n1;i<=n2;i++) {
			if(se->seqngap[i]!='-') ses[m2++]=se->seqngap[i];
			if(str->seqngap[i]!='-') strs[m1++]=str->seqngap[i];			
		}
		ses[m2]='\0';
		strs[m1]='\0';

		//set up matrix
		
		float **matrix=0;
		matrix=new float*[m1+2];
		for(i=0;i<m1;i++)matrix[i]=new float[m2]; 
 		matrix[m1]=0;
		matrix[m1+1]=0;
		int j;
		for(i=0;i<m1;i++)
		for(j=0;j<m2;j++) matrix[i][j]=0;
		//check the number residues in end loops
		int x1,x2;
		int y1,y2;

		x1=x2=0;
		for(i=n1;i<g1;i++) {			
			if(str->seqngap[i]!='-') x1++;	
			if(se->seqngap[i]!='-')  x2++; 	
		}
		
		y1=y2=0;
		for(i=g2+1;i<=n2;i++) {			
			if(str->seqngap[i]!='-') y1++;	
			if(se->seqngap[i]!='-')  y2++; 	
		}

		int mm=0;
		for(i=g1;i<=g2;i++) mm++;	

		//check if shift required
		int h1;

		int ii;

		

		float emax=0;
		float hox=0;
		for(i=g1;i<=g2;i++) {
			
			if(se->seqngap[i]=='-'&&str->seqngap[i]=='-') continue;
			if(se->seqngap[i]=='-') {
				int s1=TRES.getresid(se->seqngap[i]);
				if(s1==-1) continue;
				emax+=blosum65mt[s1*(s1+1)/2+s1];
			}
			else if(str->seqngap[i]=='-'){
				int s1=TRES.getresid(str->seqngap[i]);
				if(s1==-1) continue;
				emax+=blosum65mt[s1*(s1+1)/2+s1];
			}
			else {
				int s1=TRES.getresid(se->seqngap[i]);
				if(s1==-1) continue;
				int s2=TRES.getresid(str->seqngap[i]);
				if(s2==-1) continue;
				if(s1>s2) hox+=blosum65mt[s1*(s1+1)/2+s2];
				else	  hox+=blosum65mt[s2*(s2+1)/2+s1];
				emax+=max(blosum65mt[s1*(s1+1)/2+s1],blosum65mt[s2*(s2+1)/2+s2]);	
			}
		}
		if(emax==0) emax=1;
		hox=hox/emax;

		float d;
		if(t0->last)  d=TRES.distance(t0->atm->next,r0->atm->next);
		else          d=-1;

		for(ii=1;ii<=shift;ii++) {//left shift

			h1=x2-ii;
					 	 
			if(h1*3.3<d) continue;

			int j;
			for(i=0;i<m1;i++) 
			for(j=0;j<m2;j++) {
				matrix[i][j]=0;
				int s1=TRES.getresid(strs[i]);
				int s2=TRES.getresid(ses[j]);
				if(s1==-1||s2==-1) continue;
				//if(s1<s2) matrix[i][j]=blosum65mt[s1*(s1+1)/2+s2];
				//else	  matrix[i][j]=blosum65mt[s2*(s2+1)/2+s1];
				if(s1>s2) matrix[i][j]=blosum65mt[s1*(s1+1)/2+s2];
				else	  matrix[i][j]=blosum65mt[s2*(s2+1)/2+s1];				
			}	

 
			j=x2-ii-1;
			float pox=0;emax=0;
			for(i=x1;i<m1-y1;i++) {
				j++;
				if(j>=m2) continue;
				matrix[i][j]=10000;
				int s1=TRES.getresid(strs[i]);
				int s2=TRES.getresid(ses[j]);
				if(s1==-1&&s2==-1) continue;
				else if(s1==-1) {
					emax+=blosum65mt[s2*(s2+1)/2+s2];
				}
				else if(s2==-1) {
					emax+=blosum65mt[s1*(s1+1)/2+s1];
				}
				else {
					if(s1>s2) pox+=blosum65mt[s1*(s1+1)/2+s2];
					else	  pox+=blosum65mt[s2*(s2+1)/2+s1];
					emax+=max(blosum65mt[s1*(s1+1)/2+s1],blosum65mt[s2*(s2+1)/2+s2]);
				}
			} 
			if(emax==0) emax=1;
			pox=pox/emax;

			if(hox-pox>dcut) continue;			

			Algn aa;
			aa.opncost=3*blosum65mt[9];
			aa.gapcost=0.2; 
			aa.matrix=matrix;
			aa.setsequence(strs,ses); 
			aa.alignment();
			int slen=aa.getroutelength();
			char **result=aa.output(stderr);
			 		
			s->seqngap=new char[nlen+slen];
			s->token=strdup("structure");
			s->pdb=new Pdb(spdb);
			s->pdb->configure();
			s->pdb->setlastresorder();
			s->more=new StrFmt;
			s->more->seqngap=new char[nlen+slen];
			s->more->token=strdup("sequence");
			s->more->more=new StrFmt;
			s->more->more->seqngap=new char[nlen+slen];
			s->more->more->token=strdup("shift");
			Res *f1=0,*f2=0;

			//find the starting fitting residues
			for(f1=t0;f1;f1=f1->last) { 		
				int h=mutate->compare[f1->id0];
				if(ner==1)  h=mutate->owner->match[h];
				else	    h=mutate->match[h];
				if(h==-1) continue;
				Res *fh1=0;
				if(ner==1) {//owner
					fh1=mutate->owner->resn[h];	
				}
				else {
					fh1=mutate->resn[h];
				}
				if(fh1==0) continue;
				float dis=f1->directrmsdanyway(fh1,0,3);
				if(dis<0.2) break; 
			}			
			if(f1==0) f1=t0->chn->res;

			//find the ending fitting residues
			for(f2=t;f2;f2=f2->next) { 		
				int h=mutate->compare[f2->id0];
				if(ner==1)  h=mutate->owner->match[h];	
				else	    h=mutate->match[h];
				if(h==-1) continue;
				Res *fh1=0;
				if(ner==1) {//owner
					fh1=mutate->owner->resn[h];	
				}
				else { //self
					fh1=mutate->resn[h];
				}
				if(fh1==0) continue;
				float dis=f2->directrmsdanyway(fh1,0,3);
				if(dis<0.2) break; 
			}			
			if(f2==0) f2=t0->chn->lastres();

			Res *sf1=s->pdb->chn->isres0(f1->id0);
			Res *sf2=s->pdb->chn->isres0(f2->id0);
			if(sf1==0||sf2==0) {aa.matrix=0;continue;}
			if(ner==1) {
				Chn chn;				
				int h1=mutate->compare[f1->id0];
				int h2=mutate->compare[f2->id0];
				if(h1==-1||h2==-1) {aa.matrix=0;continue;}
				int s1=mutate->owner->match[h1];
				int s2=mutate->owner->match[h2];
				//if(s1==-1||s2==-1)  {aa.matrix=0;continue;}
				Res *ff1=0; if(s1!=-1)ff1=mutate->owner->resn[s1];
				Res *ff2=0; if(s2!=-1)ff2=mutate->owner->resn[s2];
				//if(ff1==0||ff2==0) {aa.matrix=0;continue;}
				chn.create(ff1,ff2->id0); 
				chn.setlastresorder();
				Res *df0=sf1->last;
				Res *df1=sf2->next;
				if(df0) { 		
					df0->next=chn.res;
					chn.lastres()->next=df1;
					chn.res=0;	
				}	
				else if(df1) {
					chn.lastres()->next=df1;
					s->pdb->chn->res=chn.res;					
					chn.res=0;
				}
				else {
					delete pdb->chn->res;
					pdb->chn->res=chn.res;
					chn.res=0;
				}
				sf2->next=0;
				if(sf1) delete sf1;sf1=0;
				Res *f;
				Res *f0=0;
				s->pdb->chn->pdb=s->pdb;
				int ne=0;int nd=0;
				for(f=s->pdb->chn->res;f;f=f->next) {
					f->chn=s->pdb->chn;
					f->last=f0;
					f->id0=ne;
					f->id=ne;
					ne++;
					Atm *a;
					for(a=f->atm;a;a=a->next) {
						a->id0=nd;
						nd++;
						a->res=f;
					}
				} 				
			}			
			s->pdb->configure();
			s->pdb->write("ss");

			int b1=mutate->compare[f1->id0];
			int b2=mutate->compare[f2->id0];

			int tt=0;
			//already done
			for(i=0;i<b1;i++) {
				s->seqngap[tt]=se->seqngap[i];
				s->more->seqngap[tt]=se->seqngap[i];
				s->more->more->seqngap[tt]='-';
				tt++;
			}

			//goback buffer segment
			for(i=b1;i<n1;i++) {
				if(ner) {
					s->seqngap[tt]=str->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				else {
					s->seqngap[tt]=se->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				tt++;
			}
			
			for(i=0;i<slen;i++) {			
				s->seqngap[tt]=result[0][i];
				s->more->seqngap[tt]=result[1][i];
				s->more->more->seqngap[tt]='r';
				tt++;
			}

			for(i=n2+1;i<b2;i++) {
				if(ner) {
					s->seqngap[tt]=str->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				else {
					s->seqngap[tt]=se->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				tt++;
			}

			for(i=b2+1;i<nlen;i++) {
				s->seqngap[tt]=se->seqngap[i];
				s->more->seqngap[tt]=se->seqngap[i];
				s->more->more->seqngap[tt]='-';
				tt++;
			}
			s->seqngap[tt]='\0';
			s->more->seqngap[tt]='\0';
			s->more->more->seqngap[tt]='\0';
			aa.matrix=0;
			s->next=new StrFmt;
			s=s->next;
			
		}

		//right shift
		if(r->next)  d=TRES.distance(t->atm->next,r->atm->next);
		else	     d=-1;

		for(ii=1;ii<=shift;ii++) {//right shift

			h1=y2-ii;
					 	 
			if(h1*3.3<d) continue;

			int j;
			for(i=0;i<m1;i++) 
			for(j=0;j<m2;j++) {
				matrix[i][j]=0;
				int s1=TRES.getresid(strs[i]);
				int s2=TRES.getresid(ses[j]);
				if(s1==-1||s2==-1) continue;
				//if(s1<s2) matrix[i][j]=blosum65mt[s1*(s1+1)/2+s2];
				//else	  matrix[i][j]=blosum65mt[s2*(s2+1)/2+s1];
				if(s1>s2) matrix[i][j]=blosum65mt[s1*(s1+1)/2+s2];
				else	  matrix[i][j]=blosum65mt[s2*(s2+1)/2+s1];				
			}	
 
			//j=m2-y2+ii-1;
			j=m2-y2+ii;	
			float pox=0;emax=0;
			for(i=m1-y1-1;i>=x1;i--) {
				j--;
				if(j<0) continue;
				matrix[i][j]=10000;
				int s1=TRES.getresid(strs[i]);
				int s2=TRES.getresid(ses[j]);
				if(s1==-1&&s2==-1) continue;
				else if(s1==-1) {
					emax+=blosum65mt[s2*(s2+1)/2+s2];
				}
				else if(s2==-1) {
					emax+=blosum65mt[s1*(s1+1)/2+s1];
				}
				else {
					if(s1>s2) pox+=blosum65mt[s1*(s1+1)/2+s2];
					else	  pox+=blosum65mt[s2*(s2+1)/2+s1];
					emax+=max(blosum65mt[s1*(s1+1)/2+s1],blosum65mt[s2*(s2+1)/2+s2]);
				}
			} 
					
			if(emax==0) emax=1;
			pox=pox/emax;

			if(hox-pox>dcut) continue;		

			Algn aa;
			aa.opncost=3*blosum65mt[9];
			aa.gapcost=0.2;  
			aa.matrix=matrix;
			aa.setsequence(strs,ses); 
			aa.alignment();
			int slen=aa.getroutelength();
			char **result=aa.output(stderr);
			 		
			s->seqngap=new char[nlen+slen];
			s->token=strdup("structure");
			s->pdb=new Pdb(spdb);
			s->pdb->configure();
			s->pdb->setlastresorder();
			s->more=new StrFmt;
			s->more->seqngap=new char[nlen+slen];
			s->more->token=strdup("sequence");
			s->more->more=new StrFmt;
			s->more->more->seqngap=new char[nlen+slen];
			s->more->more->token=strdup("shift");
			Res *f1=0,*f2=0;

			//find the starting fitting residues
			for(f1=t0;f1;f1=f1->last) { 		
				int h=mutate->compare[f1->id0];
				if(ner==1)  h=mutate->owner->match[h];
				else	    h=mutate->match[h];
				if(h==-1) continue;
				Res *fh1=0;
				if(ner==1) {//owner
					fh1=mutate->owner->resn[h];	
				}
				else {
					fh1=mutate->resn[h];
				}
				if(fh1==0) continue;
				float dis=f1->directrmsdanyway(fh1,0,3);
				if(dis<0.2) break; 
			}
			
			if(f1==0) f1=t0->chn->res;

			//find the ending fitting residues
			for(f2=t;f2;f2=f2->next) { 		
				int h=mutate->compare[f2->id0];
				if(ner==1)  h=mutate->owner->match[h];	
				else	    h=mutate->match[h];
				if(h==-1) continue;
				Res *fh1=0;
				if(ner==1) {//owner
					fh1=mutate->owner->resn[h];	
				}
				else { //self
					fh1=mutate->resn[h];
				}
				if(fh1==0) continue;
				float dis=f2->directrmsdanyway(fh1,0,3);
				if(dis<0.2) break; 
			}
			
			if(f2==0) f2=t0->chn->lastres();
		 
			Res *sf1=s->pdb->chn->isres0(f1->id0);
			Res *sf2=s->pdb->chn->isres0(f2->id0);
			if(sf1==0||sf2==0)  {aa.matrix=0;continue;}
			if(ner==0) {
				Chn chn;				
				int h1=mutate->compare[f1->id0];
				int h2=mutate->compare[f2->id0];
				if(h1==-1||h2==-1)  {aa.matrix=0;continue;}
				int s1=mutate->owner->match[h1];
				int s2=mutate->owner->match[h2];
				if(s1==-1||s2==-1)  {aa.matrix=0;continue;}
				Res *ff1=mutate->owner->resn[s1];
				Res *ff2=mutate->owner->resn[s2];
				if(ff1==0||ff2==0)  {aa.matrix=0;continue;}
				chn.create(ff1,ff2->id0); 
				chn.setlastresorder();
				Res *df0=sf1->last;
				Res *df1=sf2->next;
				if(df0) { 		
					df0->next=chn.res;
					chn.lastres()->next=df1;
					chn.res=0;	
				}	
				else if(df1) {
					chn.lastres()->next=df1;
					s->pdb->chn->res=chn.res;					
					chn.res=0;
				}
				else {
					delete pdb->chn->res;
					pdb->chn->res=chn.res;
					chn.res=0;
				}
				sf2->next=0;
				if(sf1) delete sf1;sf1=0;
				Res *f;
				Res *f0=0;
				s->pdb->chn->pdb=s->pdb;
				int ne=0;int nd=0;
				for(f=s->pdb->chn->res;f;f=f->next) {
					f->chn=s->pdb->chn;
					f->last=f0;
					f->id0=ne;
					f->id=ne;
					ne++;
					Atm *a;
					for(a=f->atm;a;a=a->next) {
						a->id0=nd;
						nd++;
						a->res=f;	
					}
				} 
				
			}			
			s->pdb->configure();
			s->pdb->write("ss");
			int b1=mutate->compare[f1->id0];
			int b2=mutate->compare[f2->id0];
			
			int tt=0;
			//already done
			for(i=0;i<b1;i++) {
				s->seqngap[tt]=se->seqngap[i];
				s->more->seqngap[tt]=se->seqngap[i];
				s->more->more->seqngap[tt]='-';
				tt++;
			}

			//goback buffer segment
			for(i=b1;i<n1;i++) {
				if(ner) {
					s->seqngap[tt]=str->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				else {
					s->seqngap[tt]=se->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				tt++;
			}
			
			for(i=0;i<slen;i++) {			
				s->seqngap[tt]=result[0][i];
				s->more->seqngap[tt]=result[1][i];
				s->more->more->seqngap[tt]='r';
				tt++;
			}

			for(i=n2+1;i<=b2;i++) {
				if(ner) {
					s->seqngap[tt]=str->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				else {
					s->seqngap[tt]=se->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				tt++;
			}

			for(i=b2+1;i<nlen;i++) {
				s->seqngap[tt]=se->seqngap[i];
				s->more->seqngap[tt]=se->seqngap[i];
				s->more->more->seqngap[tt]='-';
				tt++;
			}
			s->seqngap[tt]='\0';
			s->more->seqngap[tt]='\0';
			s->more->more->seqngap[tt]='\0';
			aa.matrix=0;
			s->next=new StrFmt;
			s=s->next;
		}
		cc.floatdel(matrix);
		cc.strdel(ses);
		cc.strdel(strs);
		r0=r;
	}	
	
	cc.floatdel(blosum65mt);
	if(other==0) return;

	StrFmt *c,*c0;
	
	
	if(other&&other->more==0) {delete other;other=0;}

	//delete nonexist
	c0=0;
	c=other;
	while(c) {
		if(c->pdb==0) {
			if(c0==0)  {
				other=c->next;
				c->next=0;
				delete c;
				c=other;
				continue;
			}
			else {
				c0->next=c->next;
				c->next=0;
				delete c;
				c=c0->next;
				continue;
			}	
		}
		else {
			char *s1=c->getseqn();
			char *s2=c->pdb->chn->getseqn();
			if(strcmp(s1,s2)) {
				if(c0==0)  {
					other=c->next;
					c->next=0;
					delete c;
					c=other;
					continue;
				}
				else {
					c0->next=c->next;
					c->next=0;
					delete c;
					c=c0->next;
					continue;
				}	
			}
		}		
		c0=c;
		c=c->next;
	}

	for(c=other;c;c=c->next) {
		if(c->next) { 
			if(c->next->more==0) {
				delete c->next;
				c->next=0;
			}
		}
	}
	//
	for(c=other;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
		if(strcmp(c0->token,"structure")==0&&c0->pdb&&c0->pdb->chn) {
			c0->cid=c0->pdb->chn->id;
			c0->start=c0->pdb->chn->res->id0+c0->pdb->chn->start;
			c0->end=c0->pdb->chn->lastres()->id0+c0->pdb->chn->start;
		}
	}
	//assign source code

	int nn=0;
	for(c=other;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
		
		if(strcmp(c0->token,"sequence")!=0) continue;
		c0->source='I';	
		c0->zipcode='-';
		char line[100];nn++;
		sprintf(line,"%s_shift_%i",se->code,nn);
		c0->code=strdup(line);	
		StrFmt *f=c->getparentStrFmt();
		f->code=strdup(line);			
	}

	
	if(other)other->tune=0;
	
	//
 	for(c=other;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
	    if(strcmp(c0->token,"sequence")!=0&&strcmp(c0->token,"structure")!=0)continue;
	    if(c0->token==0) {
		cerr<<"\n\n\nis it structure or sequence for the following pir:\n\n"<<endl;
		c0->printoutonly(stderr);
		exit(0);
	    }
	    if(strcmp(c0->token,"sequence")==0&&strchr("FI",c0->source)&&c0->code==0) {
		cerr<<"\n\n\nplease provide the name of the sequence:"<<endl;
		c0->printoutonly(stderr);
		exit(0);
	    }
	    if(strcmp(c0->token,"structure")==0&&(c0->pdb==0&&c0->code==0)) {
		cerr<<"\n\n\nplease provide the name of the structure:"<<endl;
		c0->printoutonly(stderr);
		exit(0);
	    }
	}

	
}

void StrFmt::setotherfmtloop(int ner,float dcut,int shift) {

	StrFmt *parent=getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
	
        if(se==0) return;
	
	//if(other) delete other; other=0;
	
	if(mutate==0) return;

	if(mutate->mpdb==0) return;
	
	StrFmt *sy,*sy0;

	sy=0;
	sy0=0;

	for(sy=other;sy;sy=sy->next) {
		sy0=sy;
	}

	if(sy0==0)  {
		other=new StrFmt;
		sy=other;
	}
	else {
		sy0->next=new StrFmt;
		sy=sy0->next; 
	} 
	    
	Pdb *spdb=spdb=mutate->mpdb;
	mutate->sethbonddssp(4);  
	mutate->setdssp(1); 
	Res *r0;

	StrFmt *str=mutate->owner;

	StrFmt *s;

	int nlen=strlen(seqngap);

	s=sy;
	Algn algn;
	float *blosum65mt=algn.getscoretable(matrix);
	algn.opncost=3*blosum65mt[9];
	algn.gapcost=0.2;
	Strhandler cc;
	spdb->setlastresorder();
	mutate->owner->pdb->setlastresorder();
	for(r0=spdb->chn->res;r0;r0=r0->next) {

		if(r0->sec=='-') continue;

 		Res *r;
		for(r=r0;r;r=r->next) {
			if(r->sec=='-') break;			 
		}	
		if(r) r=r->last;
		
		Res *t0, *t;

		for(t0=r0->last;t0;t0=t0->last) {
			if(t0->sec!='-') break;
			if(r0->id0-t0->id0>shift+5) break;
		}			
		if(t0==0) t0=spdb->chn->res;
		else t0=t0->next;

		for(t=r->next;t;t=t->next) {
			if(t->sec!='-') break;
			if(t->id0-r->id0>shift+5) break;
		}
			
		if(t==0) t=spdb->chn->lastres();
		else t=t->last;
 
		//test
		Res *of1=0;
		Res *of2=0;
		//find the starting fitting residues
		for(of1=t0;of1;of1=of1->last) { 		
			int h=mutate->compare[of1->id0];
			//if(ner==1)  h=mutate->owner->match[h];
			//else	    h=mutate->match[h];
			h=mutate->owner->match[h];
			if(h==-1) continue;
			Res *fh1=0;
			//if(ner==1) {//owner
			fh1=mutate->owner->resn[h];	
			//}
			//else {
			//	fh1=mutate->resn[h];
			//}
			if(fh1==0) continue;
			float dis=of1->directrmsdanyway(fh1,0,3);
			if(dis<0.2) break; 
		}			
		
		//find the ending fitting residues
		for(of2=t;of2;of2=of2->next) { 		
			int h=mutate->compare[of2->id0];
			//if(ner==1)  h=mutate->owner->match[h];	
			//else	    h=mutate->match[h];
			h=mutate->owner->match[h];
			if(h==-1) continue;
			Res *fh1=0;
			//if(ner==1) {//owner
			fh1=mutate->owner->resn[h];	
			//}
			//else { //self
			//	fh1=mutate->resn[h];
			//}
			if(fh1==0) continue;
			float dis=of2->directrmsdanyway(fh1,0,3);
			if(dis<0.2) break; 
		}			
		if(of2==0||of1==0) {
			if(ner==1) {r0=r;continue;}
		}
		else if(of1&&of2) {
			if(ner==2) {r0=r;continue;}
		}
		//end
		

		char *ses,*strs;

		int n=t->id0-t0->id0+100;

		ses=new char[n];
		strs=new char[n];
 		
		//find out end residues id
		int n1=mutate->compare[t0->id0];
		int n2=mutate->compare[t->id0];
	
		int g1=mutate->compare[r0->id0];
		int g2=mutate->compare[r->id0];
		int stem[2];
		stem[0]=n1;
		stem[1]=n2;
		if(TRES.logg)mutate->printsegment(stem);
		//take out sequence
		int i; 
		int m1=0;
		int m2=0;
		for(i=n1;i<=n2;i++) {
			if(se->seqngap[i]!='-') ses[m2++]=se->seqngap[i];
			if(str->seqngap[i]!='-') strs[m1++]=str->seqngap[i];			
		}
		ses[m2]='\0';
		strs[m1]='\0';

		//set up matrix
		
		float **matrix=0;
		matrix=new float*[m1+2];
		for(i=0;i<m1;i++)matrix[i]=new float[m2]; 
 		matrix[m1]=0;
		matrix[m1+1]=0;
		int j;
		for(i=0;i<m1;i++)
		for(j=0;j<m2;j++) matrix[i][j]=0;
		//check the number residues in end loops
		int x1,x2;
		int y1,y2;

		x1=x2=0;
		for(i=n1;i<g1;i++) {			
			if(str->seqngap[i]!='-') x1++;	
			if(se->seqngap[i]!='-')  x2++; 	
		}
		
		y1=y2=0;
		for(i=g2+1;i<=n2;i++) {			
			if(str->seqngap[i]!='-') y1++;	
			if(se->seqngap[i]!='-')  y2++; 	
		}

		int mm=0;
		for(i=g1;i<=g2;i++) mm++;	

		//check if shift required
		int h1;

		int ii;

		

		float emax=0;
		float hox=0;
		for(i=g1;i<=g2;i++) {
			
			if(se->seqngap[i]=='-'&&str->seqngap[i]=='-') continue;
			if(se->seqngap[i]=='-') {
				int s1=TRES.getresid(se->seqngap[i]);
				if(s1==-1) continue;
				emax+=blosum65mt[s1*(s1+1)/2+s1];
			}
			else if(str->seqngap[i]=='-'){
				int s1=TRES.getresid(str->seqngap[i]);
				if(s1==-1) continue;
				emax+=blosum65mt[s1*(s1+1)/2+s1];
			}
			else {
				int s1=TRES.getresid(se->seqngap[i]);
				if(s1==-1) continue;
				int s2=TRES.getresid(str->seqngap[i]);
				if(s2==-1) continue;
				if(s1>s2) hox+=blosum65mt[s1*(s1+1)/2+s2];
				else	  hox+=blosum65mt[s2*(s2+1)/2+s1];
				emax+=max(blosum65mt[s1*(s1+1)/2+s1],blosum65mt[s2*(s2+1)/2+s2]);	
			}
		}
		if(emax==0) emax=1;
		hox=hox/emax;

		float d;
		if(t0->last)  d=TRES.distance(t0->atm->next,r0->atm->next);
		else          d=-1;

		for(ii=1;ii<=shift;ii++) {//left shift

			h1=x2-ii;
					 	 
			if(h1*3.3<d) continue;

			int j;
			for(i=0;i<m1;i++) 
			for(j=0;j<m2;j++) {
				matrix[i][j]=0;
				int s1=TRES.getresid(strs[i]);
				int s2=TRES.getresid(ses[j]);
				if(s1==-1||s2==-1) continue;
				//if(s1<s2) matrix[i][j]=blosum65mt[s1*(s1+1)/2+s2];
				//else	  matrix[i][j]=blosum65mt[s2*(s2+1)/2+s1];
				if(s1>s2) matrix[i][j]=blosum65mt[s1*(s1+1)/2+s2];
				else	  matrix[i][j]=blosum65mt[s2*(s2+1)/2+s1];				
			}	

 
			j=x2-ii-1;
			float pox=0;emax=0;
			for(i=x1;i<m1-y1;i++) {
				j++;
				if(j>=m2) continue;
				matrix[i][j]=10000;
				int s1=TRES.getresid(strs[i]);
				int s2=TRES.getresid(ses[j]);
				if(s1==-1&&s2==-1) continue;
				else if(s1==-1) {
					emax+=blosum65mt[s2*(s2+1)/2+s2];
				}
				else if(s2==-1) {
					emax+=blosum65mt[s1*(s1+1)/2+s1];
				}
				else {
					if(s1>s2) pox+=blosum65mt[s1*(s1+1)/2+s2];
					else	  pox+=blosum65mt[s2*(s2+1)/2+s1];
					emax+=max(blosum65mt[s1*(s1+1)/2+s1],blosum65mt[s2*(s2+1)/2+s2]);
				}
			} 
			if(emax==0) emax=1;
			pox=pox/emax;

			if(hox-pox>dcut) continue;			

			Algn aa;
			aa.opncost=3*blosum65mt[9];
			aa.gapcost=0.2; 
			aa.matrix=matrix;
			aa.setsequence(strs,ses); 
			aa.alignment();
			int slen=aa.getroutelength();
			char **result=aa.output(stderr);
			 		
			s->seqngap=new char[nlen+slen];
			s->token=strdup("structure");
			s->pdb=new Pdb(spdb);
			s->pdb->configure();
			s->pdb->setlastresorder();
			s->more=new StrFmt;
			s->more->seqngap=new char[nlen+slen];
			s->more->token=strdup("sequence");
			s->more->more=new StrFmt;
			s->more->more->seqngap=new char[nlen+slen];
			s->more->more->token=strdup("shift");
			Res *f1=0,*f2=0;

			//find the starting fitting residues
			for(f1=t0;f1;f1=f1->last) { 		
				int h=mutate->compare[f1->id0];
				if(ner==1)  h=mutate->owner->match[h];
				else	    h=mutate->match[h];
				if(h==-1) continue;
				Res *fh1=0;
				if(ner==1) {//owner
					fh1=mutate->owner->resn[h];	
				}
				else {
					fh1=mutate->resn[h];
				}
				if(fh1==0) continue;
				float dis=f1->directrmsdanyway(fh1,0,3);
				if(dis<0.2) break; 
			}			
			if(f1==0) f1=t0->chn->res;

			//find the ending fitting residues
			for(f2=t;f2;f2=f2->next) { 		
				int h=mutate->compare[f2->id0];
				if(ner==1)  h=mutate->owner->match[h];	
				else	    h=mutate->match[h];
				if(h==-1) continue;
				Res *fh1=0;
				if(ner==1) {//owner
					fh1=mutate->owner->resn[h];	
				}
				else { //self
					fh1=mutate->resn[h];
				}
				if(fh1==0) continue;
				float dis=f2->directrmsdanyway(fh1,0,3);
				if(dis<0.2) break; 
			}			
			if(f2==0) f2=t0->chn->lastres();

			Res *sf1=s->pdb->chn->isres0(f1->id0);
			Res *sf2=s->pdb->chn->isres0(f2->id0);
			if(sf1==0||sf2==0) {aa.matrix=0;continue;}
			if(ner==1) {
				Chn chn;				
				int h1=mutate->compare[f1->id0];
				int h2=mutate->compare[f2->id0];
				if(h1==-1||h2==-1) {aa.matrix=0;continue;}
				int s1=mutate->owner->match[h1];
				int s2=mutate->owner->match[h2];
				//if(s1==-1||s2==-1)  {aa.matrix=0;continue;}
				Res *ff1=0; if(s1!=-1)ff1=mutate->owner->resn[s1];
				Res *ff2=0; if(s2!=-1)ff2=mutate->owner->resn[s2];
				//if(ff1==0||ff2==0) {aa.matrix=0;continue;}
				chn.create(ff1,ff2->id0); 
				chn.setlastresorder();
				Res *df0=sf1->last;
				Res *df1=sf2->next;
				if(df0) { 		
					df0->next=chn.res;
					chn.lastres()->next=df1;
					chn.res=0;	
				}	
				else if(df1) {
					chn.lastres()->next=df1;
					s->pdb->chn->res=chn.res;					
					chn.res=0;
				}
				else {
					delete pdb->chn->res;
					pdb->chn->res=chn.res;
					chn.res=0;
				}
				sf2->next=0;
				if(sf1) delete sf1;sf1=0;
				Res *f;
				Res *f0=0;
				s->pdb->chn->pdb=s->pdb;
				int ne=0;int nd=0;
				for(f=s->pdb->chn->res;f;f=f->next) {
					f->chn=s->pdb->chn;
					f->last=f0;
					f->id0=ne;
					f->id=ne;
					ne++;
					Atm *a;
					for(a=f->atm;a;a=a->next) {
						a->id0=nd;
						nd++;
						a->res=f;
					}
				} 				
			}			
			s->pdb->configure();
			s->pdb->write("ss");

			int b1=mutate->compare[f1->id0];
			int b2=mutate->compare[f2->id0];

			int tt=0;
			//already done
			for(i=0;i<b1;i++) {
				s->seqngap[tt]=se->seqngap[i];
				s->more->seqngap[tt]=se->seqngap[i];
				s->more->more->seqngap[tt]='-';
				tt++;
			}

			//goback buffer segment
			for(i=b1;i<n1;i++) {
				if(ner) {
					s->seqngap[tt]=str->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				else {
					s->seqngap[tt]=se->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				tt++;
			}
			
			for(i=0;i<slen;i++) {			
				s->seqngap[tt]=result[0][i];
				s->more->seqngap[tt]=result[1][i];
				s->more->more->seqngap[tt]='r';
				tt++;
			}

			for(i=n2+1;i<b2;i++) {
				if(ner) {
					s->seqngap[tt]=str->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				else {
					s->seqngap[tt]=se->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				tt++;
			}

			for(i=b2+1;i<nlen;i++) {
				s->seqngap[tt]=se->seqngap[i];
				s->more->seqngap[tt]=se->seqngap[i];
				s->more->more->seqngap[tt]='-';
				tt++;
			}
			s->seqngap[tt]='\0';
			s->more->seqngap[tt]='\0';
			s->more->more->seqngap[tt]='\0';
			aa.matrix=0;
			s->next=new StrFmt;
			s=s->next;
			
		}

		//right shift
		if(r->next)  d=TRES.distance(t->atm->next,r->atm->next);
		else	     d=-1;

		for(ii=1;ii<=shift;ii++) {//right shift

			h1=y2-ii;
					 	 
			if(h1*3.3<d) continue;

			int j;
			for(i=0;i<m1;i++) 
			for(j=0;j<m2;j++) {
				matrix[i][j]=0;
				int s1=TRES.getresid(strs[i]);
				int s2=TRES.getresid(ses[j]);
				if(s1==-1||s2==-1) continue;
				//if(s1<s2) matrix[i][j]=blosum65mt[s1*(s1+1)/2+s2];
				//else	  matrix[i][j]=blosum65mt[s2*(s2+1)/2+s1];
				if(s1>s2) matrix[i][j]=blosum65mt[s1*(s1+1)/2+s2];
				else	  matrix[i][j]=blosum65mt[s2*(s2+1)/2+s1];				
			}	
 
			//j=m2-y2+ii-1;
			j=m2-y2+ii;	
			float pox=0;emax=0;
			for(i=m1-y1-1;i>=x1;i--) {
				j--;
				if(j<0) continue;
				matrix[i][j]=10000;
				int s1=TRES.getresid(strs[i]);
				int s2=TRES.getresid(ses[j]);
				if(s1==-1&&s2==-1) continue;
				else if(s1==-1) {
					emax+=blosum65mt[s2*(s2+1)/2+s2];
				}
				else if(s2==-1) {
					emax+=blosum65mt[s1*(s1+1)/2+s1];
				}
				else {
					if(s1>s2) pox+=blosum65mt[s1*(s1+1)/2+s2];
					else	  pox+=blosum65mt[s2*(s2+1)/2+s1];
					emax+=max(blosum65mt[s1*(s1+1)/2+s1],blosum65mt[s2*(s2+1)/2+s2]);
				}
			} 
					
			if(emax==0) emax=1;
			pox=pox/emax;

			if(hox-pox>dcut) continue;		

			Algn aa;
			aa.opncost=3*blosum65mt[9];
			aa.gapcost=0.2;  
			aa.matrix=matrix;
			aa.setsequence(strs,ses); 
			aa.alignment();
			int slen=aa.getroutelength();
			char **result=aa.output(stderr);
			 		
			s->seqngap=new char[nlen+slen];
			s->token=strdup("structure");
			s->pdb=new Pdb(spdb);
			s->pdb->configure();
			s->pdb->setlastresorder();
			s->more=new StrFmt;
			s->more->seqngap=new char[nlen+slen];
			s->more->token=strdup("sequence");
			s->more->more=new StrFmt;
			s->more->more->seqngap=new char[nlen+slen];
			s->more->more->token=strdup("shift");
			Res *f1=0,*f2=0;

			//find the starting fitting residues
			for(f1=t0;f1;f1=f1->last) { 		
				int h=mutate->compare[f1->id0];
				if(ner==1)  h=mutate->owner->match[h];
				else	    h=mutate->match[h];
				if(h==-1) continue;
				Res *fh1=0;
				if(ner==1) {//owner
					fh1=mutate->owner->resn[h];	
				}
				else {
					fh1=mutate->resn[h];
				}
				if(fh1==0) continue;
				float dis=f1->directrmsdanyway(fh1,0,3);
				if(dis<0.2) break; 
			}
			
			if(f1==0) f1=t0->chn->res;

			//find the ending fitting residues
			for(f2=t;f2;f2=f2->next) { 		
				int h=mutate->compare[f2->id0];
				if(ner==1)  h=mutate->owner->match[h];	
				else	    h=mutate->match[h];
				if(h==-1) continue;
				Res *fh1=0;
				if(ner==1) {//owner
					fh1=mutate->owner->resn[h];	
				}
				else { //self
					fh1=mutate->resn[h];
				}
				if(fh1==0) continue;
				float dis=f2->directrmsdanyway(fh1,0,3);
				if(dis<0.2) break; 
			}
			
			if(f2==0) f2=t0->chn->lastres();
		 
			Res *sf1=s->pdb->chn->isres0(f1->id0);
			Res *sf2=s->pdb->chn->isres0(f2->id0);
			if(sf1==0||sf2==0)  {aa.matrix=0;continue;}
			if(ner==0) {
				Chn chn;				
				int h1=mutate->compare[f1->id0];
				int h2=mutate->compare[f2->id0];
				if(h1==-1||h2==-1)  {aa.matrix=0;continue;}
				int s1=mutate->owner->match[h1];
				int s2=mutate->owner->match[h2];
				if(s1==-1||s2==-1)  {aa.matrix=0;continue;}
				Res *ff1=mutate->owner->resn[s1];
				Res *ff2=mutate->owner->resn[s2];
				if(ff1==0||ff2==0)  {aa.matrix=0;continue;}
				chn.create(ff1,ff2->id0); 
				chn.setlastresorder();
				Res *df0=sf1->last;
				Res *df1=sf2->next;
				if(df0) { 		
					df0->next=chn.res;
					chn.lastres()->next=df1;
					chn.res=0;	
				}	
				else if(df1) {
					chn.lastres()->next=df1;
					s->pdb->chn->res=chn.res;					
					chn.res=0;
				}
				else {
					delete pdb->chn->res;
					pdb->chn->res=chn.res;
					chn.res=0;
				}
				sf2->next=0;
				if(sf1) delete sf1;sf1=0;
				Res *f;
				Res *f0=0;
				s->pdb->chn->pdb=s->pdb;
				int ne=0;int nd=0;
				for(f=s->pdb->chn->res;f;f=f->next) {
					f->chn=s->pdb->chn;
					f->last=f0;
					f->id0=ne;
					f->id=ne;
					ne++;
					Atm *a;
					for(a=f->atm;a;a=a->next) {
						a->id0=nd;
						nd++;
						a->res=f;	
					}
				} 
				
			}			
			s->pdb->configure();
			s->pdb->write("ss");
			int b1=mutate->compare[f1->id0];
			int b2=mutate->compare[f2->id0];
			
			int tt=0;
			//already done
			for(i=0;i<b1;i++) {
				s->seqngap[tt]=se->seqngap[i];
				s->more->seqngap[tt]=se->seqngap[i];
				s->more->more->seqngap[tt]='-';
				tt++;
			}

			//goback buffer segment
			for(i=b1;i<n1;i++) {
				if(ner) {
					s->seqngap[tt]=str->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				else {
					s->seqngap[tt]=se->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				tt++;
			}
			
			for(i=0;i<slen;i++) {			
				s->seqngap[tt]=result[0][i];
				s->more->seqngap[tt]=result[1][i];
				s->more->more->seqngap[tt]='r';
				tt++;
			}

			for(i=n2+1;i<=b2;i++) {
				if(ner) {
					s->seqngap[tt]=str->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				else {
					s->seqngap[tt]=se->seqngap[i];
					s->more->seqngap[tt]=se->seqngap[i];
					s->more->more->seqngap[tt]='-';
				}
				tt++;
			}

			for(i=b2+1;i<nlen;i++) {
				s->seqngap[tt]=se->seqngap[i];
				s->more->seqngap[tt]=se->seqngap[i];
				s->more->more->seqngap[tt]='-';
				tt++;
			}
			s->seqngap[tt]='\0';
			s->more->seqngap[tt]='\0';
			s->more->more->seqngap[tt]='\0';
			aa.matrix=0;
			s->next=new StrFmt;
			s=s->next;
		}
		cc.floatdel(matrix);
		cc.strdel(ses);
		cc.strdel(strs);
		r0=r;
	}	
	
	cc.floatdel(blosum65mt);
	if(other==0) return;

	
	if(other&&other->more==0) {delete other;other=0;}

	//delete nonexist
   	StrFmt *c0,*c;
	c0=0;
	c=other;
	while(c) {
		if(c->pdb==0) {
			if(c0==0)  {
				other=c->next;
				c->next=0;
				delete c;
				c=other;
				continue;
			}
			else {
				c0->next=c->next;
				c->next=0;
				delete c;
				c=c0->next;
				continue;
			}	
		}
		else {
			char *s1=c->getseqn();
			char *s2=c->pdb->chn->getseqn();
			if(strcmp(s1,s2)) {
				if(c0==0)  {
					other=c->next;
					c->next=0;
					delete c;
					c=other;
					continue;
				}
				else {
					c0->next=c->next;
					c->next=0;
					delete c;
					c=c0->next;
					continue;
				}	
			}
		}		
		c0=c;
		c=c->next;
	}

	for(c=other;c;c=c->next) {
		if(c->next) { 
			if(c->next->more==0) {
				delete c->next;
				c->next=0;
			}
		}
	}
	//
	for(c=other;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
		if(strcmp(c0->token,"structure")==0&&c0->pdb&&c0->pdb->chn) {
			c0->cid=c0->pdb->chn->id;
			c0->start=c0->pdb->chn->res->id0+c0->pdb->chn->start;
			c0->end=c0->pdb->chn->lastres()->id0+c0->pdb->chn->start;
		}
	}
	//assign source code

	int nn=0;
	for(c=other;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
		
		if(strcmp(c0->token,"sequence")!=0) continue;
		c0->source='I';	
		c0->zipcode='-';
		char line[100];nn++;
		sprintf(line,"%s_shift_%i",se->code,nn);
		c0->code=strdup(line);	
		StrFmt *f=c->getparentStrFmt();
		f->code=strdup(line);			
	}

	
	other->tune=0;
	
	//
 	for(c=other;c;c=c->next)
        for(c0=c;c0;c0=c0->more) {
	    if(strcmp(c0->token,"sequence")!=0&&strcmp(c0->token,"structure")!=0)continue;
	    if(c0->token==0) {
		cerr<<"\n\n\nis it structure or sequence for the following pir:\n\n"<<endl;
		c0->printoutonly(stderr);
		exit(0);
	    }
	    if(strcmp(c0->token,"sequence")==0&&strchr("FI",c0->source)&&c0->code==0) {
		cerr<<"\n\n\nplease provide the name of the sequence:"<<endl;
		c0->printoutonly(stderr);
		exit(0);
	    }
	    if(strcmp(c0->token,"structure")==0&&(c0->pdb==0&&c0->code==0)) {
		cerr<<"\n\n\nplease provide the name of the structure:"<<endl;
		c0->printoutonly(stderr);
		exit(0);
	    }
	}

	
}

char *StrFmt::getseqn() {

	int n=strlen(seqngap);
	int i;
	char *line=new char[n+1];
	int jj=0;
	for(i=0;i<n;i++) {
		if(seqngap[i]=='-') continue;
		line[jj++]=seqngap[i];
	}
	line[jj]='\0';
	return line;
}
void StrFmt::smoothinggap() {
 
	cerr<<endl<<"remove zig-zag gaps..."<<endl;

	StrFmt *parent=getparentStrFmt();
        StrFmt *se = parent->findstructurefmt();
        if(se==0) return;

	
	int nlen=strlen(seqngap);
	
	int i,m;
	


	char pt1[100],pt2[100];
	
 	int numc=0;
	re200:
		 
	numc++;
	if(numc>50) goto re250;			
	int n5;n5=0;

	for(i=0;i<nlen;i++) {			
			
		int n1=0;
		int n2=0;
		int n3=0;
		int n4=0;	
		int it=0;		
		for(m=i;m<nlen;m++) {
			if(m>=nlen) break;
			if(se->seqngap[m]=='-'&&seqngap[m]=='-') continue;
			if(it>5) continue;
			it++;
			if(se->seqngap[m]!='-') pt1[n1++]=se->seqngap[m];
			if(seqngap[m]!='-') pt2[n2++]=seqngap[m];
			if(se->seqngap[m]!='-'||seqngap[m]!='-') n3++;	
			if(n1==n2&&n1<n3) {
				n4=m;
				break;
			}
			else {
				n4=m;
			}
		}
		
		
		if(n1==n2&&n1<n3) {				
			n5=1;
			for(m=i;m<=n4;m++) {
				if(m>=nlen) break;
				se->seqngap[m]='-';
				seqngap[m]='-';
			}
			int nn=0;
			for(m=i;m<=n4;m++) {
				if(m>=nlen) break;
				if(nn<n1) {
					se->seqngap[m]=pt1[nn];
					seqngap[m]=pt2[nn];
				} 
				nn++;	
			}
		}
	}				
	se->setmatch();
	setmatch();

	if(n5) goto re200;

	re250:
	se->setmatch();
	setmatch();
}
 

void StrFmt::setmatch() {

	//reset the compare table from alignment to structure

        if(seqngap==0) return;

        int n=strlen(seqngap);

        int i;

        for(i=0;i<n;i++) {
		match[i]=-1;
		compare[i]=-1;
		//resn[i]=0;	
	}

        int m=0;

        for(i=0;i<n;i++) {
                if(seqngap[i]=='-') continue;
                match[i]=m;
                compare[m]=i;
		seqn[m]=seqngap[i];
                //resn[m]=pdb->chn->isres0(m);
                m++;
        }
	seqn[m]='\0';

	

        return;
}

void StrFmt::reallocatmatch(int nn) {

	//reset the compare table from alignment to structure

        if(seqngap==0) return;
	
	if(match) {
		match=(int *)realloc(match,sizeof(int)*nn);
		
	}
	if(compare) {
		compare=(int *)realloc(compare,sizeof(int)*nn);
	}
	if(seqngap) {
		int n=strlen(seqngap);
		seqngap=(char *)realloc(seqngap,sizeof(char)*nn);
		seqngap[n]='\0';
	}
	if(seqn) {
		int n=strlen(seqngap);
		seqn=(char *)realloc(seqn,sizeof(char)*nn);
		seqn[n]='\0';		
	}
	if(resn) {
		resn=(Res **)realloc(resn,sizeof(Res *)*nn);		
	}

	 
        int n=strlen(seqngap);
	
        int i;

        for(i=0;i<n;i++) {
		match[i]=-1;
		compare[i]=-1;
		//resn[i]=0;	
	}

        int m=0;

        for(i=0;i<n;i++) {
                if(seqngap[i]=='-') continue;
                match[i]=m;
                compare[m]=i;
		seqn[m]=seqngap[i];
                //resn[m]=pdb->chn->isres0(m);
                m++;
        }
	seqn[m]='\0';
	 
        return;
}

void StrFmt::setdssp() {

	
	StrFmt *parent=getparentStrFmt();

	StrFmt *str=parent->findstructurefmt();

	int n=strlen(seqngap);


	char *dssp=new char[n+1];

	int i;

	for(i=0;i<n;i++) dssp[i]='?';

	n=strlen(seqngap);

	for(i=0;i<n;i++) {
		int j=str->match[i];
		if(j==-1) continue;
		Res *r=str->resn[j];
		if(r==0) continue;
		
		dssp[i]=r->sec;
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
	for(i=0;i<n;i++) {
		int j=match[i];
		if(j==-1) continue;
		resn[j]->sec=dssp[i];
	}
	if(dssp) delete [] dssp; dssp=0;
	if(TRES.logg) cerr<<seqngap<<endl;
	if(TRES.logg) cerr<<str->seqngap<<endl;
	if(TRES.logg) cerr<<dssp<<endl;
}
char *StrFmt::setnewdsspsimple() {
	
	StrFmt *parent=getparentStrFmt();
	StrFmt *str=parent->findstructurefmt();
	
	int n=strlen(seqngap);

	char *dssp=new char[n+1];

	int i;

	for(i=0;i<n;i++) dssp[i]='-';
	if(str==0) return dssp;
	n=strlen(seqngap);

	for(i=0;i<n;i++) {
		int j=str->match[i];
		if(j==-1) continue;
		Res *r=str->resn[j];
		if(r==0) continue;		
		dssp[i]=r->sec;
	}
	
	dssp[n]='\0';
	if(TRES.logg) cerr<<str->seqngap<<endl;
	if(TRES.logg) cerr<<seqngap<<endl;
	if(TRES.logg) cerr<<dssp<<endl;

	return dssp;
}
char *StrFmt::setnewdssp() {
	
	StrFmt *parent=getparentStrFmt();
	StrFmt *str=parent->findstructurefmt();
	
	int n=strlen(seqngap);

	char *dssp=new char[n+1];

	int i;

	for(i=0;i<n;i++) dssp[i]='?';
	if(str==0) return dssp;
	n=strlen(seqngap);

	for(i=0;i<n;i++) {
		int j=str->match[i];
		if(j==-1) continue;
		Res *r=str->resn[j];
		if(r==0) continue;		
		dssp[i]=r->sec;
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
	if(TRES.logg) cerr<<str->seqngap<<endl;
	if(TRES.logg) cerr<<seqngap<<endl;
	if(TRES.logg) cerr<<dssp<<endl;

	return dssp;
}
void StrFmt::printhelp() {
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"nest is a homology model building program with the following capabilities:\n");
fprintf(stderr,"A. model building based on the given alignment\n");
fprintf(stderr,"B. construction of composite structure\n");  
fprintf(stderr,"C. model building based on multiple templates\n");
fprintf(stderr,"D. structure refinement\n");
fprintf(stderr,"\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"Usage: nest -seed num -fast num  -tune num -log str -opt num -nopt num  -out num file.pir\n");
fprintf(stderr,"-fast     fast mode.from 0-5; default is 3.with 5 fastest\n");
fprintf(stderr,"-tune     alignment tuning, from 0-3; defautl is 0.\n");
fprintf(stderr,"          0,no alignment tuning\n");
fprintf(stderr,"          1,remove gaps in secondary strucutre\n");
fprintf(stderr,"          2,remove zig-zag gaps  plus -tune 1;\n"); 
fprintf(stderr,"          3,remove terminal gaps plus -tune 2;\n");
fprintf(stderr,"-opt      refine mode, from 0-4; default is 1.with 4 thorough refinement\n");
fprintf(stderr,"          0, no minimization at all;\n");
fprintf(stderr,"          1, remove atom clashes; \n");
fprintf(stderr,"          2, refine in insertion and deletion regions, usually loop regions; \n");
fprintf(stderr,"          3, refine in all loop regions;\n");
fprintf(stderr,"          4, refine in all loop and secondary regions\n");
fprintf(stderr,"-out      output mode, from 0-2, defautl is 0. with 2 all intermediate output\n");
fprintf(stderr,"          0, only final structures written down;\n");
fprintf(stderr,"          1, intermediate structures also outputted;\n");
fprintf(stderr,"          2, more intermediate output\n");
fprintf(stderr,"-seed     set seed number, defautl is 18120.\n");		
fprintf(stderr,"-nopt     number of rounds of refinement.default is 1.\n");
fprintf(stderr,"-log      log file.default is standard screen output\n");
fprintf(stderr,"file.pir  pir alignment file\n");
}

void StrFmt::printrefinehelp() {
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"conref is a program to refine protein structures with constraints from the original structure\n");
fprintf(stderr,"\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu.\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"Usage: conref -prm num -seed num -fast num  -opt num -nopt num -obj num-num -rmsd float -out num -cid char file.pdb\n");
fprintf(stderr,"-fast     fast mode.from 0-5; default is 3.with 5 fastest\n");
fprintf(stderr,"-opt      refine mode, from 1-3; default is 3.with 3 most thorough refinement\n");
fprintf(stderr,"          1, remove atom clashes; \n");
fprintf(stderr,"          2, refine in all loop regions;\n");
fprintf(stderr,"          3, refine in all loop and secondary regions\n");
fprintf(stderr,"-seed     set seed number, defautl is 18120.\n");
fprintf(stderr,"-prm      force field parameter mode. from 1-4, default is 3.\n");
fprintf(stderr,"          1, charmm all atoms; 2, amber all atoms; 3, charmm heavy atoms; 4, amber heavy atoms\n");		
fprintf(stderr,"-nopt     number of rounds of refinement.default is 1.\n");
fprintf(stderr,"-out      output mode, from 0-2, defautl is 0. with 2 all intermediate output\n");
fprintf(stderr,"          0, only final structures written down;\n");  
fprintf(stderr,"          1, intermediate structures also outputted;\n");  
fprintf(stderr,"          2, more intermediates output\n");  
fprintf(stderr,"-cid      chain id.\n");
fprintf(stderr,"-obj      specify the segment to be minimized; default is whole chain.\n");
fprintf(stderr,"-rmsd     rmsd restraint.default is 2.0A\n");
fprintf(stderr,"file.pdb  pdb file\n");
}

void StrFmt::printminhelp() {
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"minst is a program to minimize protein structures\n");
fprintf(stderr,"using smoothing energy technique\n");
fprintf(stderr,"\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu.\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"Usage: minst -prm num -opt num -fast num -obj num-num -nopt num -out num -cid char file.pdb\n");
fprintf(stderr,"-opt      refine mode, from 1-3; default is 3.with 3 most thorough refinement\n");
fprintf(stderr,"          1, remove atom clashes; \n");
fprintf(stderr,"          2, refine in all loop regions;\n");
fprintf(stderr,"          3, refine in all loop and secondary regions\n");
fprintf(stderr,"-fast     fast mode.from 0-2; default is 0.with 2 fastest\n");
fprintf(stderr,"-obj      specify the segment to be minimized; default is whole chain.\n");
fprintf(stderr,"-prm      force field parameter mode. from 1-4, default is 3.\n");
fprintf(stderr,"          1, charmm all atoms; 2, amber all atoms; 3, charmm heavy atoms; 4, amber heavy atoms\n");		
fprintf(stderr,"-nopt     number of rounds of refinement.default is 1.\n");
fprintf(stderr,"-out      output mode, from 0-1, defautl is 0. with 2 all intermediate output\n");
fprintf(stderr,"          0, only final structures written down;\n");  
fprintf(stderr,"          1, intermediate structures also outputted;\n");   
fprintf(stderr,"-cid      chain id.\n");
fprintf(stderr,"file.pdb  pdb file\n");
}

void StrFmt::printhelpold() {
printf("Usage: nest -fast {on,off} -tune {on,off} -side {on,off} -loop {on,off} -at {on,off}\n"); 
printf("	    -refine {on,off} -seed {number} -ag {on,off} -w {number} -mout {on,off}\n");  
printf("	    -shift {number} -addh {on,off} alig_file\n");
printf("-fast     fast mode.default is off.\n");
printf("-tune     alignment tuning, defautl is on.\n");
printf("-side     perform side chain prediction on modeled structure, defautl is off.\n");
printf("-loop     perform loop prediction on modeled structure, defautl is off.\n");
printf("-refine   perform structure deep refinement on modeled structure, defautl is off.\n");
printf("-seed     set seed number, defautl is 18120.\n");
printf("-mout     output structures at each stage, defautl is on.\n");
printf("-shift    perform alignment shift on secondary region, defautl is 0.\n");
printf("-ag       align segment conformations from the composite template  or multiple templates, defautl is on.\n");
printf("-at       automatic multiple-template handling.\n");
printf("-addh     modeling in all atom type.\n");
printf("-w        window size in residue number of segment to be aligned.\n");

}
