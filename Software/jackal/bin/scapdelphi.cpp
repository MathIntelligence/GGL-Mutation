#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Pdb  pdb;
Scap scprd;
Rotate side;
char listname[200],pdbname[100],rotname[200];
int nint,n,m;

//interpret the command line

nint=3;
rotname[0]='\0';
pdbname[0]='\0';
 
//setup delphi
scprd.allcharges=2;
scprd.indi=4;
scprd.delphi=1;
scprd.delphiscale=0.5;
scprd.delphionself=0;
scprd.exdi=80;
scprd.delphiself=1;
//end

scprd.singletorsion=1;
scprd.colonyline=1;
scprd.colony=2;
scprd.ncolony=1;
scprd.nncolony=1;
scprd.nummore=40;
scprd.bmax=2;
scprd.tormax=2;
scprd.ring=1;
TRES.flag=100;
int seed=18120;

m=0;
for(n=1;n<argc;n+=2)
{
  if(argv[n][0]=='-'&&m) {
	cerr<<"error in arguments"<<endl;
	scprd.printnewhelp();
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
	scprd.printnewhelp();
	exit(0);
    }
    TRES.flag=nt;
  }
  else if(!strcmp(argv[n],"-min"))//minimization
  {
   int nt=atoi(argv[n+1]);
   if(nt==1) scprd.flag+=1000;
   else if(nt==2) scprd.flag+=1000*2;
   else if(nt==3) scprd.flag+=1000*3;
   else if(nt==4) scprd.flag+=1000*4;
   else if(nt==0) continue;
   else {
	cerr<<"-min "<<nt<<" not defined!"<<endl;
        scprd.printnewhelp();
        exit(0);
   }
  }
  else if(!strcmp(argv[n],"-seed")) {
        if(n+1==argc-1) {
               cerr<<"wrong in command. check."<<endl;
               scprd.printnewhelp();
               exit(0);
        }
        seed=atoi(argv[n+1]);
        if(seed<0) {
               cerr<<"seed number must be > 0"<<endl;
               exit(0);
        }
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
        scprd.printnewhelp();
        exit(0);
    }
  }
  else if(!strcmp(argv[n],"-ini")) { //number of calculations
    nint=atoi(argv[n+1]);
    if(nint<20) scprd.arbt=0;
    else        scprd.arbt=(nint-20)/2;
  }
  else if(!strcmp(argv[n],"-self")) {
    int nt=atoi(argv[n+1]);
    if(nt==1) {
	scprd.includeself=2;
    }
    else if(nt==2) {
	scprd.includeself=1;
    }
    else if(nt==0) {
	scprd.includeself=0;
    }
    else {
	cerr<<"error: -self "<<argv[n+1]<<" not defined!"<<endl;
        scprd.printnewhelp();
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
        	scprd.printnewhelp();
        	exit(0);
	}
  }
  else if(argv[n][0]=='-') {
    cerr<<"do not know the meaning:"<<argv[n]<<endl;
    scprd.printnewhelp();
    exit(0);
  }
  else
  {
     m++;
     if(m==1) strcpy(pdbname,argv[n]);
     else if(m>=2) {
		if(strstr(argv[n],"_scap.list")) {
			strcpy(listname,argv[n]);
		}
 		else {
			cerr<<"do not know the meaning:"<<argv[n]<<endl;
			exit(0);
		}
     }
     n--;
  }
}

if(m==0) {
	cerr<<"no pdb file provided:"<<endl;
    	scprd.printnewhelp();
    	exit(0);
}
srandom(seed);
if(scprd.arbt<3) scprd.arbt=3;

if(rotname[0]=='\0') strcpy(rotname,"side_large_rotamer"); 

if(rotname[0]=='\0'||pdbname[0]=='\0'||argc==1) {
  scprd.printnewhelp();
  exit(0);
}

TRES.read("tres");
TRES.smt=TRES.nsmt-1;
if(scprd.rott==0) TRES.readmore("back_small_rotamer",rotname);

pdb.read(pdbname,'0');

pdb.setflgr(10000);
pdb.setflgr('G',-99999);

scprd.pdb=&pdb;

scprd.readlist(listname);

if(scprd.includeself==0) pdb.transform(5,-99999);

scprd.hookside();

pdb.setflgr('A',-99999);
pdb.setflgr('P',-99999);
if(TRES.flag%100==0) {
	cerr<<"adding hydrogens"<<endl;
	pdb.header();
	pdb.addhatoms(1);
	pdb.configure();
}
char *s=0;
if(pdb.name) s=strdup(pdb.name);
else         s=strdup(pdbname);
scprd.scpred(s,rotname,nint);
}
