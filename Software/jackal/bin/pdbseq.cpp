#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Pdb  pdb;
if(argc<2) {
    TRES.printxyzseqhelp();
    exit(0);
}
int f=0;
int n=0;
char ch='0';
int m=0;
char *out=0;
for(n=1;n<argc;n+=2) {

	if(!strcmp(argv[n],"-a")) {
		f=atoi(argv[n+1]);			 	 
	}
	else if(!strcmp(argv[n],"-o")) {
		out=strdup(argv[n+1]);
	}
	else if(!strcmp(argv[n],"-c")) {
		if(n+1==argc-1) {
			cerr<<"error in command. check."<<endl;
			exit(0);
		}
		int nf=strlen(argv[n+1]);
		if(nf!=1) {
			cerr<<"chain id can only be one character:"<<argv[n+1]<<endl;
			exit(0);
		}
		ch=argv[n+1][0];
        }
	else if(argv[n][0]=='-') {                
                cerr<<"do not understand the meaning:"<<argv[n]<<endl;
                TRES.printxyzseqhelp();
                exit(0);
        }
	else {
		m++;
		if(n!=argc-1) {
			cerr<<"error in command. check."<<endl;
                        TRES.printxyzseqhelp();
                        exit(0);
		}
		break;
	}
}
if(m!=1) {
	cerr<<"error in command. check."<<endl;
        TRES.printxyzseqhelp();
        exit(0);
}
TRES.flag=102;
TRES.read("tres");
pdb.read(argv[argc-1],ch);
if(f==0) {
	pdb.readseqcard();
}
else if(f==1) {
	pdb.setseqcard();
}
else if(f==2) {
	pdb.setseqcardnogap();
}
//char line[1000];
//sprintf(line,"%s.seq",pdb.name);
if(out) {
	cerr<<"write out result:"<<out<<endl;
	pdb.writeseq(out);
}
else pdb.writeseq(stderr);
exit(0);
}
