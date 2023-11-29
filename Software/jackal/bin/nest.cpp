#include"../fasp/head/source.h"

Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
StrFmt mod;
Pdb  pdb;
int out=0;
int tune=0;
int logg=1;
int fapr=3;
int seed=18120;
int sharp=1;
float rmsd=2.0;
int asloop=0;
int test=0;
TRES.flag=102; 
int seglen=2;
int nopt=1;
int restraint=1;
if(argc==1) {
	mod.printhelp();
	exit(0);
}

//create log files
int n;
for(n=1;n<argc;n++) {
	if(!strcmp(argv[n],"-log")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
			mod.printhelp();
                        exit(0);
                }
		if(argv[n+1]==0||strlen(argv[n+1])==0) {
			cerr<<"wrong in command. check."<<endl;
			mod.printhelp();
                        exit(0);
		}
		int nt=creat(argv[n+1],S_IRWXU|S_IROTH|S_IRGRP);
		dup2(nt,2);		 
        }
}

for(n=1;n<argc;n++) {
	 
	if(!strcmp(argv[n],"-tune")) {		 
		if(n+1==argc-1) {
			cerr<<"wrong in command. check."<<endl;
			mod.printhelp();
			exit(0);
		}
		int nt=atoi(argv[n+1]);
		if(nt<0||nt>3) {
			cerr<<"-tune "<<argv[n+1]<<" not defined!"<<endl;
			cerr<<"tune option ranges from 0-3"<<endl;
			exit(0);
		}
		tune=nt;
	}
	else if(!strcmp(argv[n],"-prm"))
  	{
		if(n+1==argc-1) {
			cerr<<"wrong in command. check."<<endl;
			mod.printhelp();
			exit(0);
		}
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
        		mod.printhelp();
        		exit(0);
    		}
    		TRES.flag=nt;
  	}
	else if(!strcmp(argv[n],"-nopt")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printhelp();
                        exit(0);
                }
                nopt=atoi(argv[n+1]);
        }
	else if(!strcmp(argv[n],"-restraint")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printhelp();
                        exit(0);
                }
                restraint=atoi(argv[n+1]);
        }
	else if(!strcmp(argv[n],"-log")) {
		continue;
	}
	else if(!strcmp(argv[n],"-test")) {		 
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printhelp();
                        exit(0);
                }
		test=atoi(argv[n+1]);               
        }
	else if(!strcmp(argv[n],"-seglen")) {		 
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printhelp();
                        exit(0);
                }
		seglen=atoi(argv[n+1]);               
        }
	else if(!strcmp(argv[n],"-out")) {		 
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printhelp();
                        exit(0);
                }
		int nt=atoi(argv[n+1]);
		if(nt<0||nt>2) {
			cerr<<"-out "<<argv[n+1]<<" not defined"<<endl;
			cerr<<"-out options ranges from 0-2"<<endl;
			exit(0);
		}
		out=nt;
        }
	else if(!strcmp(argv[n],"-fast")) {		 
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
			mod.printhelp();
                        exit(0);
                }
		int nt=atoi(argv[n+1]);
		if(nt<0||nt>5) {
			cerr<<"-fast "<<argv[n+1]<<" not defined!"<<endl;
			cerr<<"the fast option ranges from 0-5"<<endl;
			cerr<<"with 0 slowest and 5 fastest"<<endl;
			exit(0);
		}
		fapr=nt;
        }
	else if(!strcmp(argv[n],"-logd")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
			mod.printhelp();
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
	
	else if(!strcmp(argv[n],"-opt")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printhelp();
                        exit(0);
                }
		int nt=atoi(argv[n+1]);
		if(nt<0||nt>4) {
			cerr<<"-opt "<<argv[n+1]<<" not defined"<<endl;
			exit(0);
		}
                sharp=nt;
        }
  	else if(!strcmp(argv[n],"-rms")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printhelp();
                        exit(0);
                }
		float nt=atof(argv[n+1]);
		if(nt<0) {
			cerr<<"-rms"<<argv[n+1]<<" should >0"<<endl;
			exit(0);
		}
                rmsd=nt;
        }
	else if(!strcmp(argv[n],"-seed")) { 
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printhelp();
                        exit(0);
                }
                seed=atoi(argv[n+1]);
		if(seed<0) {
			cerr<<"seed number must be > 0"<<endl;
			exit(0);
		}
        }
	else if(!strcmp(argv[n],"-asloop")) { 
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printhelp();
                        exit(0);
                }
                asloop=atoi(argv[n+1]);		
        }
	else if(argv[n][0]=='-') { 
		cerr<<"do not understand the meaning:"<<argv[n]<<endl;
		mod.printhelp();
		exit(0);
	}
}

//cerr<<random()<<endl;
srandom(seed);
TRES.read("tres");
TRES.switchcharge(0);
TRES.surface(10,2.,-1);

TRES.smt=1;
TRES.setnewengcoeff(0,0.5);
 
TRES.readmore("back_small_rotamer","side_small_rotamer");
Pdb *rotamer=new Pdb;

rotamer->readmore("side_mix_rotamer");
rotamer->token=strdup("sidechain");
TRES.rotamer->next->next=rotamer;
TRES.prepareengtable();

TRES.readdistpopular("bound_nearnext1168.dist","internal");
TRES.readdistpopular("bound_hbondind1168.dist","hbond");


mod.tune=tune;
if(TRES.flag%10==0) mod.addh=1;
	
mod.readModelAlgnFmt(argv[argc-1]);
mod.prepare();
mod.setmutate("fapr",fapr);
mod.setmutate("seed",seed);
mod.setmutate("test",test);
mod.setmutate("sharp",sharp);
mod.setmutate("out",out);
mod.setmutate("rmsd",rmsd*100);
mod.setmutate("asloop",asloop);
mod.setmutate("seglen",seglen);
mod.setmutate("restraint",restraint);
//build
mod.exec();
int i; for(i=0;i<nopt;i++) mod.optimize();
mod.createComStrFmt();
mod.combine();
for(i=0;i<nopt;i++) mod.optimizecom();
mod.beforewritefinal();
mod.writefinal();
}

