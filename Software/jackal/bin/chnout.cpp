#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Pdb  pdb;
TRES.flag=102;
if(argc<2) {
    TRES.printchnparserhelp();
    exit(0);
}
int n=0;
char ch='0';
int m=0;
char *out=0;
int ini=0;
char cidnew='-';
for(n=1;n<argc;n+=2) {
	if(!strcmp(argv[n],"-o")) {
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
	else if(!strcmp(argv[n],"-i")) {
                if(n+1==argc-1) {
                        cerr<<"error in command. check."<<endl;
                        exit(0);
                }
                ini=atoi(argv[n+1]);
        }
	else if(!strcmp(argv[n],"-n")) {
                if(n+1==argc-1) {
                        cerr<<"error in command. check."<<endl;
                        exit(0);
                }
                cidnew=argv[n+1][0];
        }
	else if(argv[n][0]=='-') {                
                cerr<<"do not understand the meaning:"<<argv[n]<<endl;
                TRES.printchnparserhelp();
                exit(0);
        }
	else {
		m++;
		if(n!=argc-1) {
			cerr<<"error in command. check."<<endl;
                        TRES.printchnparserhelp();
                        exit(0);
		}
		break;
	}
}
if(m!=1) {
	cerr<<"error in command. check."<<endl;
        TRES.printchnparserhelp();
        exit(0);
}
TRES.read("tres");
pdb.delchn=0;
pdb.read(argv[argc-1],ch);
//char line[1000];
//sprintf(line,"%s.seq",pdb.name);
Strhandler cc;
pdb.info=cc.strdel(pdb.info);
pdb.endinfo=cc.strdel(pdb.endinfo);
Chn *cu;
for(cu=pdb.chn;cu;cu=cu->next) cu->addresoid(ini);
for(cu=pdb.chn;cu;cu=cu->next) 
for(Res *r=cu->res;r;r=r->next)
for(Atm *a=r->atm;a;a=a->next) a->occup=1;

for(cu=pdb.chn;cu;cu=cu->next) {
if(cidnew!='-') {
cu->id=cidnew;
}
}
if(out==0) {
	pdb.writeold(stdout);
}
else {
	char line[1000];
	Chn *c;
	for(c=pdb.chn;c;c=c->next) {
		if(c->more==0) {
			if(c->id!=' ') { 
				sprintf(line,"%s_%c.pdb",out,c->id);
				cerr<<"output file:"<<line<<endl;
				c->writeold(line);
			}
			else {
				sprintf(line,"%s.pdb",out);
				cerr<<"output file:"<<line<<endl;
                                c->writeold(line);
			}
		}
		else {
			Chn *s;
			int n=0;
			for(s=c;s;s=s->more) {
				if(c->id!=' ') { 
                                	sprintf(line,"%s_%c_%i.pdb",out,c->id,n);
                                	s->writeold(line);
					cerr<<"output file:"<<line<<endl;
                        	}
                        	else {
                                	sprintf(line,"%s_%i.pdb",out,n);
                                	s->writeold(line);
					cerr<<"output file:"<<line<<endl;
                        	}
				n++;
			}
		}
	}
}
exit(0);
}
