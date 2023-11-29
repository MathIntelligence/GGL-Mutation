#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Stralg algn;
//Pdb  native, pdb;
StrFmt mod;
//interpret the command line

if(argc<2) {
    TRES.printseqxyzhelp();
    exit(0);
}

int f=0;
int n;
TRES.flag=100;
for(n=1;n<argc;n+=2) {

	if(!strcmp(argv[n],"-f")) {
		if(n+1==argc-1) {
			cerr<<"error in command\n"<<endl;
			TRES.printseqxyzhelp();
                        exit(0);
		}
		f=atoi(argv[n+1]);			 	 
		if(f<0||f>1) {
			cerr<<"-f from 0-1"<<endl;
			exit(0);
		}
		if(f==0) TRES.flag=100;
		if(f==1) TRES.flag=102;
	}
	else if(argv[n][0]=='-') {                
                cerr<<"do not understand the meaning:"<<argv[n]<<endl;
                TRES.printnalgnhelp();
                exit(0);
        }
	else {
		if(n!=argc-1) {
			cerr<<"error in command. check."<<endl;
                        TRES.printseqxyzhelp();
                        exit(0);
		}
		break;
	}
}
TRES.read("tres");
mod.readModelAlgnFmt(argv[argc-1]);
mod.writeseqxyz();
}
