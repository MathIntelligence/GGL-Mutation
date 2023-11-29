#include"../fasp/head/source.h"

Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
PdbFix pdbfix;
//StrFmt mod;
Pdb *pdb=new Pdb();
int logg=0;
if(argc==1) {
	pdbfix.printhelpfix();
	exit(0);
}
TRES.flag=100;
//create log files
int lent=10000;
int n;
int fixwhat=1;
for(n=1;n<argc;n++) {
	 
	if(!strcmp(argv[n],"-prm"))
  	{
		if(n+1==argc-1) {
			cerr<<"wrong in command. check."<<endl;
			pdbfix.printhelpfix();
			exit(0);
		}
    		int nt=atoi(argv[n+1]);
		
    		if(nt==1) {
        		nt=100;
    		}
    		else if(nt==2) {
        		nt=102;
    		}
    		else {
        		cerr<<"-prm ranges from 1-4"<<endl;
        		pdbfix.printhelpfix();
        		exit(0);
    		}
    		TRES.flag=nt;
  	}
        else if(!strcmp(argv[n],"-logd")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        pdbfix.printhelpfix();
                        exit(0);
                }
                int nt=atoi(argv[n+1]);
                if(nt==999) logg=4;
                else if(nt==888) logg=3;
                else if(nt==777) logg=2;
                else if(nt==666) logg=1;
                else logg=0;
                TRES.logg=logg;
        }
	else if(!strcmp(argv[n],"-len")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
			pdbfix.printhelpfix();
                        exit(0);
                }
		lent=atoi(argv[n+1]);
        }			
	else if(!strcmp(argv[n],"-fast")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        pdbfix.printhelp();
                        exit(0);
                }
                int nt=atoi(argv[n+1]);
                if(nt<0||nt>5) {
                        cerr<<"option is from 0 to 5"<<endl;
                        cerr<<argv[n+1]<<endl;
                        exit(0);
                }
                pdbfix.fast=nt;
        }

	else if(!strcmp(argv[n],"-fix")){
		if(n+1==argc-1) {
			 cerr<<"wrong in command. check."<<endl;
			 pdbfix.printhelpfix();
			 exit(0);		
		}
		int nt=atoi(argv[n+1]);
		if(nt!=0&&nt!=1) {
			cerr<<"wrong in command. check."<<endl;
			pdbfix.printhelpfix();
			exit(0);
		}
		fixwhat=nt;	
	}				
	else if(argv[n][0]=='-') { 
		cerr<<"do not understand the meaning:"<<argv[n]<<endl;
		pdbfix.printhelpfix();
		exit(0);
	}
}

srandom(18120);
TRES.read("tres");


TRES.smt=1;
TRES.setnewengcoeff(0,0.5);
 
TRES.readmore("back_large_rotamer","side_large_rotamer");
TRES.prepareengtable();
Pdb *rotamer=new Pdb;
rotamer->readmore("side_small_rotamer");
TRES.rotamer->next->next=rotamer;

TRES.readdistpopular("bound_nearnext1168.dist","internal");
TRES.readdistpopular("bound_hbondind1168.dist","hbond");


pdb->read(argv[argc-1],'0');
char *name=0;
if(pdb->name) name=strdup(pdb->name);
else         {
	cerr<<"program detects error in protein name"<<endl;
	exit(0);
}
if(name==0){
	cerr<<"program detects error in protein name"<<endl;
        exit(0);
}
pdbfix.onlybackbone=0;
pdbfix.pdb=pdb; pdb=0;
pdbfix.pdb->configure();
pdbfix.checkpdb(pdbfix.pdb);
pdbfix.ready();
pdbfix.myfixnow();


pdbfix.pdb->info=pdbfix.pdborg->info;
pdbfix.pdborg->info=0;
pdbfix.pdb->endinfo=pdbfix.pdborg->endinfo;
pdbfix.pdborg->endinfo=0;
if(fixwhat) pdbfix.fixlostresidues(pdbfix.pdb,lent);
char nameo[1000];
sprintf(nameo,"%s_fix.pdb",name);
cerr<<"write down the final structure..."<<nameo<<endl;
pdbfix.pdb->write(nameo);
if(name) delete [] name;name=0;
}

