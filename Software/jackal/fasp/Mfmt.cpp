#include"source.h"
Mfmt::Mfmt()
{
	aln=0;
	sec=1; 
	name=0;
}

Mfmt::~Mfmt()
{
	Strhandler cc;
	 
	aln=cc.strdel(aln);
	name=cc.strdel(name);
}

 
void Mfmt::setmfmt(StrFmt *fmt) {

	cerr<<"integrate alignments to show alignment difference..."<<endl;
	cerr<<"all the alignment should be between the same query sequence and the same structure"<<endl;
	cerr<<"under different alignment methods.."<<endl;
	
	if(fmt==0) return;


	StrFmt *s;

	cerr<<"check the alignment blocks..."<<endl;
	for(s=fmt;s;s=s->next) {
		StrFmt *s1=s->findstructurefmt();	
		StrFmt *s2=s->getStrFmt("secondary");
		if(s1==0&&s2==0) {
			s->printout(stderr);
			cerr<<"none existence of secondary and structure pirs\n";
			exit(0);
		}
		else if(s1&&s2) {
			s->printout(stderr);
			cerr<<"existence of both secondary and structure pirs\n";
                        exit(0);	
		}
	}

	int nsize=0;
	for(s=fmt;s;s=s->next) nsize++;
	cerr<<endl;
	cerr<<"the total number of alignment blocks:"<<nsize<<endl;
	cerr<<endl;
	nsize=nsize*2+1000;

	cerr<<"check existence of structure,sequence and secondary pirs.."<<endl;
	StrFmt *s0=fmt->findstructurefmt();
	for(s=fmt;s;s=s->next) {
		if(s0) break;
		s0=s->findstructurefmt();
	}
	Strhandler cc;
	for(s=fmt;s;s=s->next) {
		StrFmt *t=s->findsequencefmt();
		 
		if(t&&(t->seqngap==0||t->seqn==0||t->code==0)) {
			t->printoutonly(stderr);
			cerr<<"none existence of sequence or code name for sequence pir"<<endl;
			exit(0);
		}
		if(t&&strcmp(t->code,"query")==0) {
			t->printoutonly(stderr);
			cerr<<"code name for sequence pir could not be \"query\". \"query\" is the reserved code"<<endl;
			exit(0);
		}
		t=s->findstructurefmt();
		 
		if(t&&(t->seqngap==0||t->seqn==0)) {
			t->printoutonly(stderr);
			cerr<<"none existence of sequence for structure pir"<<endl;
			exit(0);
		}
		t=s->getStrFmt("secondary");
		if(t&&(t->seqngap==0)) {
			t->printoutonly(stderr);
			cerr<<"none existence of sequence for structure pir"<<endl;
			exit(0);
		}
		else if(t) {
			t->seqngap=cc.lowercase(t->seqngap);
		}
	}
	for(s=fmt;s;s=s->next) {
		StrFmt *str=s->getStrFmt("secondary"); 
		if(str==0) continue;
		int n=strlen(str->seqngap);
		int i;
		for(i=0;i<n;i++) {
			if(str->seqngap[i]!='e'&&str->seqngap[i]!='h') str->seqngap[i]='-';
		}
	}
	 
	s0=fmt->findsequencefmt(); 
	for(s=fmt;s;s=s->next) {
		if(s0) break;
		s0=fmt->findsequencefmt();
	}
	if(s0==0) {
		cerr<<"there is no sequence pirs.."<<endl;
		exit(0);
	}
	int nlen=strlen(s0->seqngap)+5000;	 

	int i;
	
	aln=new char*[nsize];
	for(i=0;i<nsize;i++) aln[i]=0;

	int j;
	for(j=0;j<3;j++) {
		aln[j]=new char[nlen];	
		for(i=0;i<nlen;i++) aln[j][i]='-';
	}
	int ne=strlen(s0->seqngap);
	for(i=0;i<ne;i++) aln[0][i]=s0->seqngap[i];
	s0=fmt->findstructurefmt();
	if(s0==0) {
		s0=s->getStrFmt("secondary"); 
		if(s0==0) {
			cerr<<"the alignment blocks must containt at least either structure or secondary pirs.."<<endl;
			exit(0);
		}
	}
	for(i=0;i<ne;i++) aln[1][i]=s0->seqngap[i]; 
	 
	char *dssp=s0->setnewdsspsimple();
	for(i=0;i<ne;i++) aln[2][i]=dssp[i]; 
	dssp=cc.strdel(dssp);
	for(i=0;i<3;i++) {
		aln[i][ne]='\0';
		int j;
		for(j=0;j<ne;j++) if(aln[i][j]=='-') aln[i][j]='Z';
	} 
	cerr<<endl;
	char *ss=new char[nlen];
	ss[0]='\0';
	char *ss1=new char[nlen];
	ss1[0]='\0';
 
	int nsize0=3;
	
	for(s=fmt->next;s;s=s->next) {
		ne=strlen(aln[0]);
		for(i=0;i<ne;i++) {
			if(aln[0][i]=='X'||aln[0][i]=='-') ss[i]='Z';
			else		                   ss[i]=aln[0][i];
		}
		ss[ne]='\0'; 
		StrFmt *sq=s->findsequencefmt();
		StrFmt *str=s->findstructurefmt();
		if(str==0) str=s->getStrFmt("secondary"); 
		if(str==0) {
			cerr<<"error. none existence of structure or secondary pirs..."<<endl;
			exit(0);
		}
		int ne1=strlen(sq->seqngap);
		for(i=0;i<ne1;i++) {
			if(sq->seqngap[i]=='X'||sq->seqngap[i]=='-') ss1[i]='Z';
			else		   			     ss1[i]=sq->seqngap[i];
		}
		ss1[ne1]='\0';
		Algn algn;
		 
		algn.setsequence(ss,ss1);
		algn.defalgnidz();
		int slen=algn.getroutelength();
		char **result=algn.output(stderr);

		for(i=0;i<nsize0+1;i++) {
			aln[nsize/2+i]=new char[nlen];
			aln[nsize/2+i][0]='\0';
		}

 		int n1=-1;
		int n2=-1;		
		
		int tot=0;
		for(i=0;i<slen;i++) {
			if(result[0][i]!='-'&&result[0][i]!='X')n1++;
			if(result[1][i]!='-'&&result[1][i]!='X')n2++;
			if(result[0][i]!='-'&&result[0][i]!='X'&&result[1][i]!='-'&&result[1][i]!='X') {
				int m1;
				for(m1=0;m1<nsize0-1;m1++) {
					aln[nsize/2+m1][tot]=aln[m1][n1];
				}	
				if(result[0][i]=='Z'&&result[1][i]!='Z') {
					aln[nsize/2+0][tot]=result[1][i];
				}	
				aln[nsize/2+nsize0-1][tot]=str->seqngap[n2]; 
				aln[nsize/2+nsize0][tot]=aln[nsize0-1][n1];
				tot++;
			}		
			else if(result[0][i]!='-'&&result[0][i]!='X') {
				int m1;
				for(m1=0;m1<nsize0-1;m1++) {
					aln[nsize/2+m1][tot]=aln[m1][n1];
				}			
				aln[nsize/2+nsize0-1][tot]='Z'; 
				aln[nsize/2+nsize0][tot]=aln[nsize0-1][n1];
				tot++;
			}	
			else if(result[1][i]!='-'&&result[1][i]!='X') {
				
				int m1;
				for(m1=0;m1<nsize0-1;m1++) {
					aln[nsize/2+m1][tot]='Z';
				}
				aln[nsize/2+0][tot]=result[1][i];			
				aln[nsize/2+nsize0-1][tot]=str->seqngap[n2]; 
				aln[nsize/2+nsize0][tot]='Z';
				tot++;
			}
		}
		
		for(i=0;i<nsize0;i++) { 		
 			if(aln[i]) {
				delete [] aln[i];aln[i]=0;
			}			
		}
		nsize0=nsize0+1;
		for(i=0;i<nsize0;i++) {
			aln[i]=aln[nsize/2+i];
			aln[nsize/2+i]=0;
			aln[i][tot]='\0';
		}
	}	

	for(i=0;i<nsize0;i++) {
		int n=strlen(aln[i]);
		int j;
		for(j=0;j<n;j++) if(aln[i][j]=='X'||aln[i][j]=='Z') aln[i][j]='-';
	}
	int nf=strlen(aln[0]);
	int nj=0;
	for(i=0;i<nf;i++) {
		int j;
		int m=0;
		for(j=0;j<nsize0;j++) {
			if(aln[j][i]=='-') m++;
		} 
		if(m==nsize0) continue;
		for(j=0;j<nsize0;j++) {
			aln[j][nj]=aln[j][i];
		}
		nj++;
	}
	for(i=0;i<nsize0;i++) {
		aln[i][nj]='\0';
	}
	aln[nsize0]=0;
	
	//start
	char **aln0=new char*[nsize0];
	for(i=0;i<nsize0;i++) {
		aln0[i]=aln[i];
		aln[i]=0;
	}
	for(i=0;i<nsize;i++) aln[i]=0;
	name=new char*[nsize];
        for(i=0;i<nsize;i++) name[i]=0;
        name[0]=strdup("query");
	aln[0]=aln0[0];
	int nx=1,ny=1;
	char line[1000];
	int nsec[10000];
	for(i=0;i<10000;i++) nsec[i]=0;
	for(s=fmt;s;s=s->next) {
		StrFmt *sq=s->findsequencefmt();
                name[nx]=strdup(sq->code);
		aln[nx]=aln0[ny];
		sq=s->findstructurefmt();
		if(sq==0) {
			nsec[nx]=1;
		}
		nx++;ny++;
		if(sq==0||sec==0) continue;
		sprintf(line,"%s-2d",sq->code);
		name[nx]=strdup(line);
		aln[nx]=new char[nlen];
		for(i=0;i<nlen;i++) aln[nx][i]='-';
		int nn=strlen(aln[nx-1]);
		int jj=0;
		for(i=0;i<nn;i++) {
			if(aln[nx-1][i]=='-'||aln[nx-1][i]=='X'||aln[nx-1][i]=='Z') continue;
			Res *r=sq->resn[jj];
			if(r) aln[nx][i]=r->sec;
			jj++;
		}
		aln[nx][nn]='\0';
		nx++;
	}
	for(i=ny;i<nsize0;i++) {
		aln0[i]=cc.strdel(aln0[i]);
	}
 	delete [] aln0;	
	nsize0=nx;
	aln[nsize0]=0;name[nsize0]=0;
	aln0=new char*[nsize0+1];
	char **name0=new char*[nsize0+1];
	for(i=0;i<nsize0;i++) {
		aln0[i]=aln[i];
		name0[i]=name[i];		
	}
	nx=0;
	for(i=0;i<nsize0;i++) {
		if(nsec[i]==1) {
			aln[nx]=aln0[i];
			name[nx]=name0[i];
			nx++;
		}
	}

	for(i=0;i<nsize0;i++) {
                if(nsec[i]==0) {
                        aln[nx]=aln0[i];
			name[nx]=name0[i];
			nx++;
                }
        }
	delete [] aln0;delete [] name0;
	aln[nsize0]=0;name[nsize0]=0;
}

void Mfmt::printhelp() {
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"difaln is a program to show the difference among multiple individual alignments\n");
fprintf(stderr,"\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"Usage: difaln -sec num file.pir\n");
fprintf(stderr,"-sec: 1 print out secondary structure\n");
fprintf(stderr,"-sec: 0 not print out secondary structure\n");
fprintf(stderr,"file.pir  pir alignment file in nest format\n");
}
 
void Mfmt::writeout(int f) {

	int n=0;
	while(aln[n]) n++;
	int i;
	 
	for(i=0;i<n;i++) {
		if(f%10==1&&i==n-1) continue;
		if(f/10==1) {
			cout<<">"<<name[i]<<endl;
			cout<<aln[i]<<endl;
		}		 
	}
	
	if(f/10==2) {
		int nj=0;
		for(i=0;i<n;i++) {
			nj=max(nj,strlen(name[i]));
		}
		nj=nj+8;
		for(i=0;i<n;i++) {
			char *s=name[i];
			s=new char[nj];
			strcpy(s,name[i]);
			int ni=strlen(s);
 			int j;			 
		}
	 	int ne=strlen(aln[0]);
		int j;
	}
}

void Mfmt::setfmt(StrFmt *fmt) {
	
	cerr<<"integrate alignments to show alignment difference..."<<endl;
	cerr<<"all the alignment should be between the same query sequence and the same structure"<<endl;
	cerr<<"under different alignment methods.."<<endl;
	char *sec=0;
	if(fmt==0) return;

	StrFmt *s;
	int nsize=0;
	for(s=fmt;s;s=s->next) nsize++;
	cerr<<endl;
	cerr<<"the total number of alignment blocks:"<<nsize<<endl;
	cerr<<endl;
	nsize=nsize+100;
	
	cerr<<"check existence of structure and sequence pirs.."<<endl;
	StrFmt *s0=fmt->findstructurefmt();
	for(s=fmt;s;s=s->next) {
		StrFmt *t=s->findsequencefmt();
		if(t->seqngap==0||t->seqn==0||t->code==0) {
			t->printoutonly(stderr);
			cerr<<"none existence of sequence or code name"<<endl;
			exit(0);
		}
		t=s->findstructurefmt();
		if(t->seqngap==0||t->seqn==0) {
			t->printoutonly(stderr);
			cerr<<"none existence of sequence"<<endl;
			exit(0);
		}
		if(strcmp(t->code,s0->code)) {
			t->printoutonly(stderr);
			s0->printoutonly(stderr);
			cerr<<"the structure should be the same..."<<endl;
			exit(0);
		} 
	}

	s0=fmt->findsequencefmt(); 
	int nlen=strlen(s0->seqngap)+5000;	 
	sec=new char[nlen];

	int i;
	for(i=0;i<nlen;i++) sec[i]='-';
	
	aln=new char*[nsize];
	for(i=0;i<nsize;i++) aln[i]=0;

	aln[0]=new char[nlen];
	for(i=0;i<nlen;i++) aln[0][i]='-';
	cerr<<endl;
	char *ss=strdup(s0->seqn);
	cerr<<"find the largest possible query sequence..."<<endl;	
	for(s=fmt;s;s=s->next) {
		Algn algn;
		StrFmt *sq=s->findsequencefmt();
		algn.setsequence(ss,sq->seqn);
		algn.defalgnid();
		int slen=algn.getroutelength();
		char **result=algn.output(stderr);
		if(ss) delete [] ss;ss=0;
		ss=new char[nlen];
		for(i=0;i<nlen;i++) ss[i]='-'; 
		for(i=0;i<slen;i++) {
			if(result[0][i]==result[1][i]) {
				if(result[0][i]=='-'||result[0][i]=='X')continue;
				ss[i]=result[0][i];
			}
			else if(result[0][i]!=result[1][i]) {
				
				if(result[0][i]!='-'&&result[0][i]!='X') {
					ss[i]=result[0][i];
				} 
				else if(result[1][i]!='-'&&result[1][i]!='X') {
					ss[i]=result[1][i];
				}
			} 
		}
		slen=0;	
		for(i=0;i<nlen;i++)  {
			if(ss[i]=='-') continue;
			ss[slen]=ss[i];slen++;
		}
		ss[slen]='\0';
	}
	cerr<<endl;
	cerr<<"the query sequence after combines from all alignments should be:"<<endl;
	cerr<<ss<<endl;
	cerr<<endl;
	//first 	 
	nsize=0;
	aln[nsize]=new char[nlen];
	for(i=0;i<nlen;i++) aln[0][i]='-';
	for(i=0;i<nlen;i++) aln[0][i]=ss[i];

	int *compare=new int[nlen];
	int *match=new int[nlen];
	for(i=0;i<nlen;i++) {
		compare[i]=-1;
		match[i]=-1;
	}
	int ne=0;
	for(i=0;i<nlen;i++) {
		if(aln[0][i]=='-') continue;
		match[i]=ne;
		compare[ne]=i;
		ne++;
	}

	for(s=fmt;s;s=s->next) {
		Algn algn;
		StrFmt *sq=s->findsequencefmt();
		StrFmt *str=s->findstructurefmt();
		algn.setsequence(ss,sq->seqn);
		algn.defalgnid();
		int slen=algn.getroutelength();
		int n1,n2;
		n1=0;
        	n2=0;
		char **result=algn.output(stderr);
		char *tt=new char[nlen];
		for(i=0;i<nlen;i++) tt[i]='-';
		for(i=0;i<slen;i++) {
			tt[i]=result[0][i];
		}
		for(i=0;i<slen;i++) {
			if(result[0][i]!='-'&&result[0][i]!='X')n1++;
			if(result[1][i]!='-'&&result[1][i]!='X')n2++;
		}
		slen=0;	
		for(i=0;i<nlen;i++)  {
			if(ss[i]=='-') continue;
			ss[slen]=ss[i];slen++;
		}
		ss[slen]='\0';
	}
}
