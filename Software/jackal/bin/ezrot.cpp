#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Pdb  pdb;
if(argc<2) {
    TRES.printezchixyzhelp();
    exit(0);
}
int n=0;
 
int rotdir=1;
int id;
float ang=999;
for(n=1;n<argc;n+=2) {

	if(!strcmp(argv[n],"-a")) {
                if(n+1==argc||n+1==argc-1) {
                        cerr<<"error in commands"<<endl;
                        TRES.printezchixyzhelp();
                        exit(0);
                }
                id=atoi(argv[n+1]);		
        }
	else if(!strcmp(argv[n],"-b")) {
                if(n+1==argc||n+1==argc-1) {
                        cerr<<"error in commands"<<endl;
                        TRES.printezchixyzhelp();
                        exit(0);
                }
                ang=atof(argv[n+1]);		
        }
	else if(!strcmp(argv[n],"-r")) {
                if(n+1==argc||n+1==argc-1) {
                        cerr<<"error in commands"<<endl;
                        TRES.printezchixyzhelp();
                        exit(0);
                }
                rotdir=atoi(argv[n+1]);
                if(rotdir<0||rotdir>1) {
                        cerr<<"error: -r "<<argv[n+1]<<endl;
                        cerr<<"-n ranges from 0-1"<<endl;
                        exit(0);
                }

        }
	else if(argv[n][0]=='-') {                
                cerr<<"do not understand the meaning:"<<argv[n]<<endl;
                TRES.printezchixyzhelp();
                exit(0);
        }
	else {
		if(n!=argc-1) {
			cerr<<"error in command. check."<<endl;
                        TRES.printezchixyzhelp();
                        exit(0);
		}
		break;
	}
}

if(ang>500) {
	cerr<<"the angle must be set between 0, 360"<<endl;
	exit(0);
}

TRES.flag=100;
TRES.read("tres");
TRES.readmore("back_small_rotamer","side_small_rotamer");
Chiangle chiangle;
chiangle.omega=0;
chiangle.rotdir=rotdir;
pdb.read(argv[argc-1],'0');
if(pdb.chn==0||pdb.chn->res==0) {
	cerr<<"pdb file:"<<argv[argc-1]<<" does not exist or has no residues!"<<endl;
	exit(0);
}
chiangle.pdb=&pdb;
chiangle.setuprecord(id,ang);
if(chiangle.number==0)  exit(0);
chiangle.buildallstructure();
if(chiangle.pdb) {
	char line0[1000];
	sprintf(line0,"%s_rot.pdb",pdb.name);
	chiangle.pdb->writeold(line0);
	cerr<<"write down the new strucutre..."<<line0<<endl;
}
chiangle.pdb=0;
exit(0);
}
