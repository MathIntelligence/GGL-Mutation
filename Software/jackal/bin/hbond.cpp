#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Pdb  pdb;
if(argc<2) {
    TRES.printhbondhelp();
    exit(0);
}
int f=0;
int n=0;
char *out=0;
for(n=1;n<argc;n+=2) {

	if(!strcmp(argv[n],"-h")) {
		f=atoi(argv[n+1]);			 	 
		if(!(f>=0&&f<=2)) {
			TRES.printhbondhelp();
                        exit(0);
		}
	}
	else if(!strcmp(argv[n],"-o")) {
                if(argv[n+1]==0) {
                        cerr<<endl<<endl<<"not defined: -a "<<argv[n+1]<<endl<<endl<<endl;
                        TRES.printchihelp();
                        exit(0);
                }
                out=strdup(argv[n+1]);
        }

	else if(argv[n][0]=='-') {                
                cerr<<"do not understand the meaning:"<<argv[n]<<endl;
                TRES.printhbondhelp();
                exit(0);
        }
	else {
		if(n!=argc-1) {
			cerr<<"error in command. check."<<endl;
                        TRES.printhbondhelp();
                        exit(0);
		}
		break;
	}
}

TRES.flag=100;
TRES.read("tres");
pdb.read(argv[argc-1],'0');
pdb.header();
Chn *c;
for(c=pdb.chn;c;c=c->next) if(f==0||f==1) c->buildhbond();
for(c=pdb.chn;c;c=c->next) if(f==0||f==2) c->buildssbond();
//char line[1000];
//sprintf(line,"%s.hbond",pdb.name);
if(out) {
        cerr<<"write output to:"<<out<<endl;
	pdb.writehbondlist(out);
}
else    pdb.writehbondlist(stderr);
exit(0);
}
