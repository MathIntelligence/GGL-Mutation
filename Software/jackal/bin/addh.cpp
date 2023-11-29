#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
//DistBound distbound;
Pdb  pdb;

if(argc<2) {
    TRES.printaddhhelp();
    exit(0);
}

int f=1;

int n=0;
for(n=1;n<argc;n+=2) {

	if(!strcmp(argv[n],"-k")) {
		f=atoi(argv[n+1]);
		if(f!=0&&f!=1) {
			cerr<<"error in -k option."<<endl;
			exit(0);
		}
		//n--;
	}
	else if(argv[n][0]=='-') {                
                cerr<<"do not understand the meaning:"<<argv[n]<<endl;
                TRES.printaddhhelp();
                exit(0);
        }
	else {
		if(n!=argc-1) {
			cerr<<"error in command. check."<<endl;
                        TRES.printaddhhelp();
                        exit(0);
		}
		break;
	}
}


TRES.flag=100;
TRES.read("tres");
pdb.read(argv[argc-1],'0');
pdb.header();
pdb.addhatoms(f);
pdb.configure();
char line[1000];
sprintf(line,"%s_addh.pdb",pdb.name);
pdb.configureatmid0();
pdb.setatmoid();
cerr<<"write out pdb:"<<line<<endl;
pdb.writeold(line);
}
