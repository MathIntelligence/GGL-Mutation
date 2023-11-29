#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;

main(int argc,char *argv[])
{
Strhandler cc;
Loopy segen; 
Chiangle chiangle;
Pdb pdb;
pdb.ishet=1;
segen.revs=1;
segen.fapr=1;
char name[100],loopname[1000];
char filnam[100];
char *seqn=new char[1000];
name[0]='\0';
loopname[0]='\0';
filnam[0]='\0';
seqn[0]='\0';
TRES.flag=102;
TRES.smt=1;
srandom(18120);

if(argc<4) {
	segen.printhelp();
	exit(0);
}


int oid=1;
int kk=0;
int n;
for(n=1;n<argc;n+=2)
{
  if(!strcmp(argv[n],"-obj"))  //loops to be predicted
  {
    segen.start=atoi(argv[n+1]);
    char *pot=strstr(argv[n+1],"-");
    if(pot==0) { 
	cerr<<"error in specifying start and end residues of loop"<<endl;
	segen.printhelp();
	exit(0);
    }
    char *poi=strstr(argv[n+1],"-")+1;
    if(poi==0) {
       cerr<<"error in command specifying start and end residues"<<endl;
       segen.printhelp();
       exit(0);
    }
    segen.end=atoi(poi);
  }
  else if(!strcmp(argv[n],"-oid")){//use oid to judge start and end
	oid=atoi(argv[n+1]);	
  }
  else if(!strcmp(argv[n],"-seed")){
	int ranseed=atoi(argv[n+1]);
	if(ranseed<0) {
		cerr<<"random seed should be larger than 0"<<endl;
		exit(0);
	}
	srandom(ranseed);
  }
  else if(!strcmp(argv[n],"-logd")){
	int nt =atoi(argv[n+1]);
	if(nt==666) TRES.logg=1;
	else if(nt==777) TRES.logg=2;
        else if(nt==888) TRES.logg=3;
	else if(nt==999) TRES.logg=4;
	else {
		cerr<<"wrong in log"<<endl;
		exit(0);
	}
  }
  else if(!strcmp(argv[n],"-ini"))  //number of initials 
  {
    segen.arbt=atoi(argv[n+1]);
    if(segen.arbt==0) {
	cerr<<"number of initials should not be zero:"<<argv[n]<<endl;
        exit(0);
    }
    segen.part=segen.arbt/10;
    segen.part0=segen.arbt/2;
    if(segen.part<10 ) segen.part = segen.part0;
  }
  else if(!strcmp(argv[n],"-prm")) //topology
  {
    int nt=atoi(argv[n+1]);
    if(nt==1) TRES.flag=100;
    else if(nt==2)  TRES.flag=200;
    else if(nt==3)  TRES.flag=102;    
    else if(nt==4)  TRES.flag=202; 
    else {
	cerr<<"the -prm should be from 1-4."<<endl;
	cerr<<"-prm "<<argv[n+1]<<endl;
	exit(0);
    }
  }
  else if(!strcmp(argv[n],"-out")) //number of outputs
  {
    segen.revs=atof(argv[n+1]);
    if(segen.revs<=0) {
	cerr<<"number of output should not be zero:"<<argv[n]<<endl;
        exit(0);
    }
  }
  else if(!strcmp(argv[n],"-sec")) { //secondary structure
    segen.secd=argv[n+1][0];
    if(segen.secd=='c') segen.secd='-';
    if(segen.secd!='-'&&segen.secd!='h'&&segen.secd!='e') {
	cerr<<"secondary structure could only be e or h or c"<<endl;
        exit(0);
    }
  }
  else if(!strcmp(argv[n],"-cid")) { //chain id
    segen.cid=argv[n+1][0];
  }
  else if(!strcmp(argv[n],"-res")) {
    	strcpy(seqn,argv[n+1]);    //sequence
	seqn=cc.clearendchar(seqn," \n");
	if(seqn==0) {
		cerr<<argv[n]<<" is not a sequence"<<endl;	
		exit(0);	
	} 
        seqn=cc.uppercase(seqn);
  }
  else if(strstr(argv[n],"-fast")) {
	segen.fapr=atoi(argv[n+1]);
	if(segen.fapr!=0&&segen.fapr!=1) {
		cerr<<"fast option can only be 0 or 1"<<endl;
		exit(0);
	}
  }
  else if(argv[n][0]=='-') {
	cerr<<"do not know the meaning:"<<argv[n]<<endl;
	exit(0);
  }
  else
  {
    kk++;
    if(kk==1) strcpy(name,argv[n]);         //protein name
    if(kk>=2) {
		if(strstr(argv[n],"_loopy.chi")) {
                        strcpy(loopname,argv[n]);
			segen.databaseonly=1;
                }
		
                else if(strstr(argv[n],".output")){
			strcpy(filnam,argv[n]);
		}
		else {
                        cerr<<"do not know the meaning:"<<argv[n]<<endl;
                        exit(0);
                }
    }
    n--;
  }
}

if(segen.start==-1&&segen.end==-1) {
	cerr<<"no start and end numbers specified or they are specified wrong\n";
	exit(0);
}

if(segen.revs<1) segen.revs=1;
if(segen.revs>segen.part) segen.part=segen.revs*2;
if(segen.part>segen.part0) segen.part0=segen.part*2;
if(segen.part0>segen.arbt) segen.arbt=segen.part0*2;



TRES.read("tres");
TRES.readmore("back_small_rotamer","side_small_rotamer");
TRES.setnewengcoeff(0,0.5);
if(strlen(loopname)!=0) {
	chiangle.read(loopname);
	segen.chiangle=&chiangle;
}

seqn=cc.clearendchar(seqn,"\"");

pdb.read(name,'0');

//check pdb
if(pdb.chn==0||pdb.chn->res==0) {
	cerr<<"the pdb file :"<<name<<" contains no standard residues.."<<endl;
	exit(0);
}

//add hatoms
if(TRES.flag%100==0) {
	cerr<<"adding hydrogens"<<endl;
        pdb.header();
        pdb.addhatoms(1);
        pdb.configure();
}
if(segen.cid=='-') {
  cerr<<"loopy is going to use the chain with empty space ' ' as the default to define the loop"<<endl;
  segen.cid=' ';
}
 
segen.pdb=&pdb;

//check chain exist
Chn *chn=pdb[segen.cid];
if(chn==0) {
   cerr<<"chain: "<<segen.cid<<" does not exist!"<<endl;
   exit(0); 
}

//add missing residues
//check sequence
if(seqn) {
    	
    if(strlen(seqn)>segen.end-segen.start+1) {
        cerr<<"number of residues in loop from "<<segen.start<<" to "<<segen.end<<" does not match the sequence given"<<endl;
	cerr<<"there are residue insertion..."<<endl;
    }
    else if(strlen(seqn)<segen.end-segen.start+1){
	cerr<<"number of residues in loop from "<<segen.start<<" to "<<segen.end<<" does not match the sequence given"<<endl;
        cerr<<"there are residue deletion..."<<endl;

    } 

    for(int ii=0;ii<strlen(seqn);ii++) {
        if(!TRES[seqn[ii]]) {
                cerr<<seqn[ii]<<" is not a standard protein residue"<<endl;
                exit(0);
        }
    }
    if(strlen(seqn)<2) {
	cerr<<"the number of residue specified by '-res' must be larger than 2:"<<seqn<<endl;
	exit(0); 
    }
}

//check loop stems
Res *rfrom=chn->isresoid(segen.start-1);
Res *rend=chn->isresoid(segen.end+1);
if(rfrom==0&&rend==0) {
	cerr<<"the loop does not have stems"<<endl;
	exit(0);	
}

if(rfrom&&rfrom->isbackbonecomplete()==0) {
	cerr<<"the backbone of the loop stem:"<<rfrom->name<<rfrom->oid<<" not complete"<<endl;
	exit(0);
}

if(rend&&rend->isbackbonecomplete()==0) {
	cerr<<"the backbone of the loop stem:"<<rend->name<<rend->oid<<" not complete"<<endl;
	exit(0);
}

//check residue exists
Res *r;
r=chn->isresoid(segen.start);
if(r==0) {
	cerr<<"residue does not exist with id:"<<segen.start<<" in chain:"<<segen.cid<<endl;
	exit(0);
}
r=chn->isresoid(segen.end);
if(r==0) {
	cerr<<"residue does not exist with id:"<<segen.end<<" in chain:"<<segen.cid<<endl;
	exit(0);
}

if(seqn&&seqn[0]!='\0') {
cerr<<"replace the original segment with a new segment: "<<seqn<<endl;
chn->addresiduesonly(segen.start,segen.end,seqn);
}
 
if(oid) r=chn->isresoid(segen.start);
else r=chn->isres0(segen.start);
if(r==0) {
  cerr<<"residue: "<<segen.start<<" in chain: "<<segen.cid<<" does not exist!"<<endl;
  exit(0);
}
segen.start=r->id0;
int start0=r->oid;

if(oid) r=chn->isresoid(segen.end);
else    r=chn->isres0(segen.end);
if(r==0) {
  cerr<<"residue: "<<segen.end<<" in chain: "<<segen.cid<<" does not exist!"<<endl;
  exit(0);
}
segen.end=r->id0;
int end0=r->oid;

cerr<<endl;
cerr<<"the number of candidates loopy will generate:"<<segen.arbt<<endl;
//cerr<<"the number of candidates loopy will be kept based on energy:"<<segen.part<<endl;
cerr<<endl;

segen.ready();
float **xyz_all=segen.predt0();
segen.chiangle=0;
segen.next->chiangle=0;

int test=1;

if(test) {
	TRES.surface(10,2.,-1);
	segen.pdb->chn->setflg(segen.pdb->chn->res,1000000,1); 
	segen.pdb->chn->surface(10,1.4);

	FILE *fp;
	float d,fact,fact1; 
	int i,ii;
	Res *r1;

	fact=0;
	fact1=0;
	for(r=segen.pdb->chn->res;r;r=r->next)
	{
  		if(r->id0<segen.start||r->id0>segen.end)continue;
  		fact+=r->area;fact1+=r->tres->area;
	}
	fact=fact/fact1*100;
	r=(*segen.pdb->chn)[segen.start];r1=(*segen.pdb->chn)[segen.end];
	fact1=TRES.distance(r->atm->next->xyz,r1->atm->next->xyz);

	n=0;
	while(xyz_all[n])n++;

	if(filnam[0]!='\0')
	{
  		ii=0;
  		fp=fopen(filnam,"a");
		 
  		for(i=0;i<n;i++)
  		{
    			if(ii>=segen.revs) break;
    			segen.segment->chn->transfer(xyz_all[i],1);
    			d=segen.rmsd(2);
    			if(d<0.001) continue;
			
    			fprintf(fp,"%2i %s %c %5i %5i %8.3f %8.3f %8.3f\n",i+1,name,'0',start0,end0,fact,fact1,d);
    			ii++;
  		}
  		fclose(fp);
	}
}
re200:
int kep;
kep=0;
while(xyz_all[kep])kep++;
int i;
for(i=0;i<kep;i++)
{
	delete [] xyz_all[i];xyz_all[i]=0;
}
delete [] xyz_all;
exit(0);
}
