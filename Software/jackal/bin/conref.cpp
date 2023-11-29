#include"../fasp/head/source.h"

Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
StrFmt mod;
 
int out=0;
int tune=0;
int logg=1;
int fapr=3;
int seed=18120;
int sharp=4;
float rmsd=2.0;
int asloop=0;
int test=0;
TRES.flag=102; 
int seglen=2;
int nopt=1;
int restraint=1;
int cid='1';
if(argc==1) {
	mod.printrefinehelp();
	exit(0);
}

int n;

for(n=1;n<argc;n++) {
	 
	if(!strcmp(argv[n],"-prm"))
  	{
		if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printrefinehelp();
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
        		mod.printrefinehelp();
        		exit(0);
    		}
    		TRES.flag=nt;
  	}
	else if(!strcmp(argv[n],"-nopt")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printrefinehelp();
                        exit(0);
                }
                nopt=atoi(argv[n+1]);
        }
	else if(!strcmp(argv[n],"-obj")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printrefinehelp();
                        exit(0);
                }
		mod.refstart=atoi(argv[n+1]);
    		char *poi=strstr(argv[n+1],"-")+1;
    		if(poi==0) {
       			cerr<<"error in command specifying start and end residues"<<endl;
       			mod.printhelp();
       			exit(0);
    		}
    		mod.refend=atoi(poi);                
        }
	else if(!strcmp(argv[n],"-cid")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printrefinehelp();
                        exit(0);
                }
		if(argv[n+1]==0&&strlen(argv[n+1])==0) {
			cerr<<"wrong in command. check."<<endl;
                        mod.printrefinehelp();
                        exit(0);
		}
                cid=argv[n+1][0];
		if(cid=='0') {
			cerr<<"chain id: 0 is not allowed."<<endl;
			exit(0);
		}
        }
	else if(!strcmp(argv[n],"-logd")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
			mod.printrefinehelp();
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
                        mod.printrefinehelp();
                        exit(0);
                }
		int nt=atoi(argv[n+1]);
		if(nt<=0||nt>3) {
			cerr<<"-opt "<<argv[n+1]<<" not defined"<<endl;
			exit(0);
		}
                sharp=nt;
		if(sharp>1) sharp=sharp+1;
        }
  	else if(!strcmp(argv[n],"-rmsd")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printrefinehelp();
                        exit(0);
                }
		float nt=atof(argv[n+1]);
		if(nt<0) {
			cerr<<"-rmsd"<<argv[n+1]<<" should >0"<<endl;
			exit(0);
		}
                rmsd=nt;
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
	else if(!strcmp(argv[n],"-seed")) { 
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        mod.printrefinehelp();
                        exit(0);
                }
                seed=atoi(argv[n+1]);
		if(seed<0) {
			cerr<<"seed number must be > 0"<<endl;
			exit(0);
		}
        }
	
	else if(argv[n][0]=='-') { 
		cerr<<"do not understand the meaning:"<<argv[n]<<endl;
		mod.printrefinehelp();
		exit(0);
	}
}

srandom(seed);
TRES.read("tres");


TRES.setnewengcoeff(0,0.5);
 
TRES.readmore("back_small_rotamer","side_small_rotamer");
Pdb *rotamer=new Pdb;

rotamer->readmore("side_mix_rotamer");
rotamer->token=strdup("sidechain");
TRES.rotamer->next->next=rotamer;
TRES.smt=1;
TRES.prepareengtable();

TRES.readdistpopular("bound_nearnext1168.dist","internal");
TRES.readdistpopular("bound_hbondind1168.dist","hbond");

mod.tune=tune;
if(TRES.flag%10==0) mod.addh=1;

if(cid=='1') {
	cerr<<"you did not specify the chain id"<<endl;
	cerr<<"the default is the first chain..."<<endl;
	cerr<<"refinement will be on the first chain.."<<endl;
}	

Pdb *pdb=new Pdb;
pdb->read(argv[argc-1],cid);
if(pdb->chn==0||pdb->chn->res==0) {
	cerr<<"there is no chain in pdb file:"<<argv[argc-1]<<endl;
	exit(0);
}
cid=pdb->chn->id;
mod.seqngap=pdb->chn->getseqn(); 
mod.code=strdup(pdb->name);
mod.token=strdup("structure");
mod.cid=pdb->chn->id;
mod.start=pdb->chn->res->oid;
mod.end=pdb->chn->lastres()->oid;
mod.pdb=pdb;
mod.onlyrefine=1;
mod.more=new StrFmt();
mod.more->start=pdb->chn->res->oid;
mod.more->end=pdb->chn->lastres()->oid;
mod.more->seqngap=pdb->chn->getseqn(); 
mod.more->code=new char[1000];
sprintf(mod.more->code,"%s",pdb->name);
mod.more->token=strdup("sequence");
mod.more->cid=pdb->chn->id;

//mod.readModelAlgnFmt(argv[argc-1]);
mod.prepare();
mod.setmutate("fapr",fapr);
mod.setmutate("seed",seed);
mod.setmutate("test",test);
mod.setmutate("sharp",sharp);
mod.setmutate("out",out);
mod.setmutate("rmsd",(int) (rmsd*100));
mod.setmutate("asloop",asloop);
mod.setmutate("seglen",seglen);
mod.setmutate("restraint",restraint);
//build
mod.exec();
int i; for(i=0;i<nopt;i++) mod.optimize();
mod.createComStrFmt();
mod.combine();
for(i=0;i<nopt;i++) mod.optimizecom();
mod.writefinal();
}

