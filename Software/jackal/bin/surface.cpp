#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Pdb  pdb;
int n;
float probe=1.4;
int size=200;
int flag=100;

if(argc<2) {
    TRES.printsurfacehelp();
    exit(0);
}
char *out=0;
for(n=1;n<argc;n+=2) {
	if(!strcmp(argv[n],"-prm")) {
                if(n+1==argc) {
                        cerr<<"error in command. check."<<endl;
                        TRES.printsurfacehelp();
                        exit(0);
                }
                int seed=atoi(argv[n+1]);
		if(seed==1) {
			flag=100;
		}
		else if(seed==2) {
			flag=200;
		}
		else if(n==3) {
        		flag=102;
    		}
    		else if(n==4) {
        		flag=202;  
    		}	
		else {
			cerr<<"-fprm:"<<argv[n+1]<<" not defined"<<endl;
			TRES.printsurfacehelp();
                        exit(0);
		}
        }
	else if(!strcmp(argv[n],"-probe")) {
                if(n+1==argc) {
                        cerr<<"error in command. check."<<endl;
                        TRES.printsurfacehelp();
                        exit(0);
                }
                probe=atof(argv[n+1]);
		if(probe<0) {
			cerr<<"probe radius must be larger than 0\n";
			TRES.printsurfacehelp();
                        exit(0);
		}
        }
	else if(!strcmp(argv[n],"-out")) {
                if(argv[n+1]==0) {
                        cerr<<endl<<endl<<"not defined: -a "<<argv[n+1]<<endl<<endl<<endl;
                        TRES.printsurfacehelp();
                        exit(0);
                }
                out=strdup(argv[n+1]);
        }
	else if(!strcmp(argv[n],"-size")) {
                if(n+1==argc) {
                        cerr<<"error in command. check."<<endl;
                        TRES.printsurfacehelp();
                        exit(0);
                }
                size=atoi(argv[n+1]);
		if(size<100) {
			cerr<<"size must be larger than 100\n";
                        TRES.printsurfacehelp();
                        exit(0);
		}
        }
        else if(argv[n][0]=='-') {                
                cerr<<"do not understand the meaning:"<<argv[n]<<endl;
                TRES.printsurfacehelp();
                exit(0);
        }
	else {
		if(n!=argc-1) {
			cerr<<"error in command. check."<<endl;
                        TRES.printsurfacehelp();
                        exit(0);
		}
		break;
	}
}

TRES.flag=flag;
TRES.read("tres");
TRES.smt=1;
pdb.read(argv[argc-1],'0');
int m=sqrt(size/2);
pdb.setflg(1);
pdb.surface(m,probe);
//char *s=pdb.name;
//char line[1000];
//sprintf(line,"%s.area",s);
if(out) {
cerr<<"write out result:"<<out<<endl;
pdb.writesurface(out);
}
else {
pdb.writesurface(stderr);
}
//cerr<<"the output file:"<<line<<endl;
cerr<<"the total area:"<<pdb.area<<endl;
}
