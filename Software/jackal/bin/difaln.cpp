#include"../fasp/head/source.h"

Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Mfmt fmt;
StrFmt mod;
Pdb  pdb;
TRES.flag=102; 
/*
int out=1;
int tune=0;
int logg=1;
int fapr=3;
int seed=18120;
int sharp=1;
float rmsd=2.0;
int asloop=0;
int test=0;
int seglen=2;
int nopt=1;
int restraint=1;
*/
if(argc==1) {
	fmt.printhelp();
	exit(0);
}

//create log files
int n;
int sec=1; 
for(n=1;n<argc;n++) {
	  

	if(!strcmp(argv[n],"-sec")) {
                if(n+1==argc-1) {
                        cerr<<"wrong in command. check."<<endl;
                        fmt.printhelp();
                        exit(0);
                }
                int nt=atoi(argv[n+1]);
                if(nt<0||nt>1) {
                        cerr<<"-sec "<<argv[n+1]<<" not defined!"<<endl;
                        cerr<<"tune option ranges from 0-1"<<endl;
                        exit(0);
                }
                sec=nt;
        }
	else if(argv[n][0]=='-') { 
		cerr<<"do not understand the meaning:"<<argv[n]<<endl;
		fmt.printhelp();
		exit(0);
	}
}

//cerr<<random()<<endl;
//srandom(seed);
TRES.read("tres");


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


//mod.tune=tune;
//cerr<<TRES.flag<<endl;
if(TRES.flag%10==0) mod.addh=1;
	
mod.readModelAlgnFmt(argv[argc-1]);
mod.prepare();
fmt.sec=sec;
fmt.setmfmt(&mod);
fmt.writeout(12); 
}

