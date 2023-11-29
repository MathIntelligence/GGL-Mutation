#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Stralg algn;
Pdb  native, pdb;

//interpret the command line

if(argc<3) {
    TRES.printalgnhelp();
    exit(0);
}

int resd=1;

char line[1000];
line[0]='\0';
int f=0;
int n;
for(n=1;n<argc;n+=2) {

	if(!strcmp(argv[n],"-f")) {
		f=atoi(argv[n+1]);			 	 
		if(f<0||f>3) {
			cerr<<"-f should be 0-3"<<endl;
			exit(0);
		} 
	}
	else if(!strcmp(argv[n],"-r")) {
                resd=atoi(argv[n+1]);
                if(f<0||f>1) {
                        cerr<<"-r should be 0-1"<<endl;
                        exit(0);
                }
        }
	else if(!strcmp(argv[n],"-o")) {	
		strcpy(line,argv[n+1]);			 	 
	}
	else if(argv[n][0]=='-') {                
                cerr<<"do not understand the meaning:"<<argv[n]<<endl;
                TRES.printalgnhelp();
                exit(0);
        }
	else {
		if(n!=argc-2) {
			cerr<<"error in command. check."<<endl;
                        TRES.printalgnhelp();
                        exit(0);
		}
		break;
	}
}
if(0&&strlen(line)==0) {
	cerr<<"output file must be specified\n";
	TRES.printalgnhelp();
    	exit(0);
}

TRES.read("tres");
pdb.read(argv[argc-2],'1');
native.read(argv[argc-1],'1');
algn.flag=1;
algn.superimpose(pdb.chn,native.chn,f);
float dd=algn.superimpose(pdb.chn,native.chn,f);
pdb.chn->id='A';
native.chn->id='B';
pdb.chn->next=native.chn;
if(line[0]!='\0')pdb.writeold(line);
pdb.chn->next=0;
if(resd==0) exit(0);
int nn=0;while(algn.alga[nn])nn++; 
int i=0;
for(i=0;i<nn;i++) {
	Atm *a;
	float rmsd=0;
	int   nrmsd=0;
	int flg=f;
	for(a=algn.alga[i]->atm;a;a=a->next) {
		if(flg==0&&a->tatm->id!=1) continue;
		if(flg==1&&a->tatm->name[1]!='H'&&a->tatm->id>=4) continue;
		if(flg==1&&a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>=3) continue;
		if(flg==2&&a->tatm->name[1]=='H') continue;
		if(flg==3&&a->tatm->id>2) continue;
		Atm *a0=(*algn.algb[i])[a->tatm->name];
		if(a0==0||a==0) continue;
		float d=TRES.distsqr(a->xyz,a0->xyz);
		rmsd+=d;
		nrmsd++;
	}
	if(nrmsd==0) nrmsd=1;
	cerr<<algn.alga[i]->name<<algn.alga[i]->oid<<"---"<<algn.algb[i]->name<<algn.algb[i]->oid<<":"<<sqrt(rmsd/nrmsd)<<endl;
}
cerr<<"the rmsd: "<<dd<<endl;
}
