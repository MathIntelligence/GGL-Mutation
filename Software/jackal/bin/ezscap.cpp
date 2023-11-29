#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Pdb  pdb;
Scap scprd;
Rotate side;
char pdbname[100],rotname[200];
int nint,n,m;

//interpret the command line

nint=3;
rotname[0]='\0';
pdbname[0]='\0';
scprd.singletorsion=1;
scprd.colonyline=1;
scprd.colony=2;
scprd.ncolony=1;
scprd.nncolony=1;
scprd.nummore=100;
scprd.bmax=2;
scprd.tormax=2;
scprd.ring=1;
TRES.flag=100;
scprd.cutbb=100000000;
nint=1;
scprd.arbt=0;
m=0;
int resid=-1000;
char resoldname=' ';
char resnewname=' ';
char chnid='1';
for(n=1;n<argc;n+=2)
{
  if(argv[n][0]=='-'&&m) {
	cerr<<"error in arguments"<<endl;
	scprd.printeznewhelp();
        exit(0);
  }
  
  else if(!strcmp(argv[n],"-prm"))
  {
    int nt=atoi(argv[n+1]);
    if(nt==1) {
	nt=100;
    }
    else if(nt==2) {
	nt=200;
    }
    else if(nt==3) {
	nt=102;
    }
    else if(nt==4) {
	nt=202; 
    }
    else {
	cerr<<"-prm ranges from 1-4"<<endl;
	scprd.printeznewhelp();
	exit(0);
    }
    TRES.flag=nt;
  }
  
  else if(!strcmp(argv[n],"-out")) {
    int nt=atoi(argv[n+1]);
    if(nt==1) {
    	scprd.nout=1;
    }
    else if(nt==2) {
	scprd.nout=2;
    }
    else if(nt==0) {
	scprd.nout=0;
    }
    else {
	cerr<<"-out "<<nt<<" not defined!"<<endl;
        scprd.printeznewhelp();
        exit(0);
    }
  }
  
  else if(!strcmp(argv[n],"-rtm")) {
	int nt=atoi(argv[n+1]);
	if(nt==1) {
		strcpy(rotname,"side_large_rotamer");		
	}
	else if(nt==2) {
		strcpy(rotname,"side_mix_rotamer");	
	}
	else if(nt==3) {
		strcpy(rotname,"side_medium_rotamer");	
	}
	else if(nt==4) {
		strcpy(rotname,"side_small_rotamer");	
	}
	else {
		cerr<<"error: -rotm "<<argv[n+1]<<" not defined!"<<endl;
        	scprd.printeznewhelp();
        	exit(0);
	}
  }
  else if(!strcmp(argv[n],"-res")) {
	if(n+1==argc-1||argv[n+1]==0) {
		cerr<<"error in commands"<<endl;
		scprd.printeznewhelp();
    		exit(0);
	}
	char *s=strdup(argv[n+1]);
	Strhandler cc;
	s=cc.clearendchar(s,"\r\t\n ");
	if(s==0) {
		cerr<<"error in commands"<<endl;
		scprd.printeznewhelp();
    		exit(0);
	}
	resoldname=s[0];
	if(!(resoldname>='A'&&resoldname<='Z')) {
		cerr<<"error in residue name :"<<argv[n+1]<<endl;
		exit(0);
	}	
	if(strlen(s)<=1) {
		cerr<<"error in seting -res "<<s<<endl;
		scprd.printeznewhelp();
    		exit(0);
	}
	resid=atoi(s+1);
	
	int ii=strlen(s);
	if((s[ii-1]>='A'&&s[ii-1]<='Z')||(s[ii-1]>='a'&&s[ii-1]<='z')) {
		resnewname=s[ii-1];
	}
	if(s) delete [] s;s=0;
  }
  else if(!strcmp(argv[n],"-cid")) {
	if(n+1==argc-1||argv[n+1]==0) {
		cerr<<"error in commands"<<endl;
		scprd.printeznewhelp();
    		exit(0);
	}
	chnid=argv[n+1][0];
  }
  else if(argv[n][0]=='-') {
    cerr<<"do not know the meaning:"<<argv[n]<<endl;
    scprd.printeznewhelp();
    exit(0);
  }
  else
  {
     m++;
     if(m==1) strcpy(pdbname,argv[n]);
     else if(m>=2) {
		
		cerr<<"do not know the meaning:"<<argv[n]<<endl;
		exit(0);
     }
     n--;
  }
}

if(m==0) {
	cerr<<"no pdb file provided:"<<endl;
    	scprd.printeznewhelp();
    	exit(0);
}

//if(scprd.arbt<3) scprd.arbt=3;

if(rotname[0]=='\0') strcpy(rotname,"side_large_rotamer"); 

if(rotname[0]=='\0'||pdbname[0]=='\0'||argc==1) {
  scprd.printeznewhelp();
  exit(0);
}

TRES.read("tres");
TRES.smt=TRES.nsmt-1;
if(scprd.rott==0) TRES.readmore("back_small_rotamer",rotname);

pdb.read(pdbname,'0');
if(pdb.chn==0) {
cerr<<"error in pdb file:"<<pdbname<<endl;
exit(0);
}
if(chnid=='1') chnid=pdb.chn->id;
if(resnewname==' ') resnewname=resoldname;
if(TRES[resoldname]==0) {
cerr<<"residue name:"<<resoldname<<" is not standard"<<endl;
exit(0);
}
if(TRES[resnewname]==0) {
cerr<<"residue name:"<<resnewname<<" is not standard"<<endl;
exit(0);
}

pdb.setflgr(-99999);
 
scprd.pdb=&pdb;

scprd.setlist(chnid,resid,resoldname,resnewname);

if(scprd.includeself==0) pdb.transform(5,-99999);

scprd.hookside();

pdb.setflgr('A',-99999);
pdb.setflgr('P',-99999);
pdb.setflgr('G',-99999);

if(TRES.flag%100==0) {
	cerr<<"adding hydrogens"<<endl;
	pdb.header();
	pdb.addhatoms(1);
	pdb.configure();
}
scprd.scpred(pdbname,rotname,nint);
}
