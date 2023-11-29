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
    TRES.printnalgnhelp();
    exit(0);
}
TRES.flag=102;
char line[1000];
line[0]='\0';
int f=0;
int n;
for(n=1;n<argc;n+=2) {

	if(!strcmp(argv[n],"-f")) {
		f=atoi(argv[n+1]);			 	 
		if(f<0||f>3) {
			cerr<<"-f from 0-3"<<endl;
			exit(0);
		}
	}
	else if(!strcmp(argv[n],"-o")) {	
		strcpy(line,argv[n+1]);			 	 
	}
	else if(argv[n][0]=='-') {                
                cerr<<"do not understand the meaning:"<<argv[n]<<endl;
                TRES.printnalgnhelp();
                exit(0);
        }
	else {
		if(n!=argc-1) {
			cerr<<"error in command. check."<<endl;
                        TRES.printnalgnhelp();
                        exit(0);
		}
		break;
	}
}
if(strlen(line)==0) {
	cerr<<"output file must be specified\n";
	//TRES.printalgnhelp();
    	exit(0);
}

TRES.read("tres");
mod.readModelAlgnFmt(argv[argc-1]);
mod.preparesuperimpose();
mod.structuresuperimpose(line,f);
}
