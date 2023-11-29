#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Chiangle chiangle;
Pdb  pdb;
if(argc<2) {
    TRES.printchixyzhelp();
    exit(0);
}
int n=0;
char line[1000];
line[0]='\0';
int omega=-1;
int rotdir=1;
int nee=0;
for(n=1;n<argc;n+=2) {

	if(!strcmp(argv[n],"-t")) {
		if(n+1==argc||n+1==argc-1) {
			cerr<<"error in commands"<<endl;
			TRES.printchixyzhelp();
                	exit(0);
		}
		strcpy(line,argv[n+1]);
	}
	else if(!strcmp(argv[n],"-n")) {
                if(n+1==argc||n+1==argc-1) {
                        cerr<<"error in commands"<<endl;
                        TRES.printchixyzhelp();
                        exit(0);
                }
                omega=atoi(argv[n+1]);
		if(omega<0||omega>3) {
			cerr<<"error: -n "<<argv[n+1]<<endl;
			cerr<<"-n ranges from 0-3"<<endl;
			exit(0);
		}
        }
	else if(!strcmp(argv[n],"-r")) {
                if(n+1==argc||n+1==argc-1) {
                        cerr<<"error in commands"<<endl;
                        TRES.printchixyzhelp();
                        exit(0);
                }
                rotdir=atoi(argv[n+1]);
                if(rotdir<0||rotdir>1) {
                        cerr<<"error: -r "<<argv[n+1]<<endl;
                        cerr<<"-n ranges from 0-1"<<endl;
                        exit(0);
                }
        }
	else if(!strcmp(argv[n],"-d")) {
                if(n+1==argc||n+1==argc-1) {
                        cerr<<"error in commands"<<endl;
                        TRES.printchixyzhelp();
                        exit(0);
                }
                int nt=atoi(argv[n+1]);
                if(nt<0||nt>1) {
                        cerr<<"error: -d "<<argv[n+1]<<endl;
                        cerr<<"-n ranges from 0-1"<<endl;
                        exit(0);
                }
		chiangle.wavedir=nt;
        }
	else if(!strcmp(argv[n],"-e")) {
                if(n+1==argc||n+1==argc-1) {
                        cerr<<"error in commands"<<endl;
                        TRES.printchixyzhelp();
                        exit(0);
                }
                int nt=atoi(argv[n+1]);
               
		chiangle.endres=nt;nee++;
        }
	else if(!strcmp(argv[n],"-logg")) {
		if(n+1==argc||n+1==argc-1) {
                        cerr<<"error in commands"<<endl;
                        TRES.printchixyzhelp();
                        exit(0);
                }
		int nt=atoi(argv[n+1]);
		int logg=0;
		if(nt==999) logg=4;
                else if(nt==888) logg=3;
                else if(nt==777) logg=2;
                else if(nt==666) logg=1;
                else logg=0;
                TRES.logg=logg;	
	}
	else if(argv[n][0]=='-') {                
                cerr<<"do not understand the meaning:"<<argv[n]<<endl;
                TRES.printchixyzhelp();
                exit(0);
        }
	else {
		if(n!=argc-1) {
			cerr<<"error in command. check."<<endl;
                        TRES.printchixyzhelp();
                        exit(0);
		}
		break;
	}
}
if(nee==0) {
if(chiangle.wavedir==1) chiangle.endres=1000000;
else chiangle.endres=-1000000;
}
if(strlen(line)==0) {
	cerr<<"chi angle name must be specified.."<<endl;
	exit(0);
}
if(omega==-1) {
	cerr<<"-n must be set"<<endl;
	exit(0);
}
TRES.flag=100;
TRES.read("tres");
TRES.readmore("back_small_rotamer","side_small_rotamer");

chiangle.omega=omega;
chiangle.rotdir=rotdir;
chiangle.readone(line);
pdb.read(argv[argc-1],'0');
chiangle.pdb=&pdb;
if(pdb.chn==0||pdb.chn->res==0) {
	cerr<<"pdb file:"<<argv[argc-1]<<" does not exist or has no residues!"<<endl;
	exit(0);
}
chiangle.buildallstructure();
if(chiangle.pdb) {
	char line0[1000];
	sprintf(line0,"%s_rot.pdb",pdb.name);
	cerr<<"write out:"<<line0<<endl;
	Strhandler cc;
	chiangle.pdb->info=cc.strdel(chiangle.pdb->info);
	chiangle.pdb->endinfo=cc.strdel(chiangle.pdb->endinfo);
	chiangle.pdb->writeold(line0);
}
chiangle.pdb=0;
exit(0);
}
