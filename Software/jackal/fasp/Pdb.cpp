#include"source.h"

void Pdb::initial() {
  token=0;
  info=0;
  endinfo=0;
  name = 0;;
  number=0;
  next=0;
  chn=0;
  area=0;
  isinfo=1;
  delchn=1;
  hetatm=0;
  ishet=0; //also read in hetatm atoms. 0 no, 1 yes
}

Pdb::Pdb()
{
  initial();
}

Pdb::Pdb(Pdb *s) {
  initial();
  if(s->name) name=strdup(s->name);
  number=s->number;
  ishet=s->ishet;
  delchn=s->delchn;
  if(s->chn){
	chn=s->chn->clonechain();
	Chn *c,*s1;
        for(c=chn;c;c=c->next) for(s1=c;s1;s1=s1->more) s1->pdb=this;
  }
  if(s->hetatm) {
	hetatm=s->hetatm->clonechain();
	Chn *c,*s1;
        for(c=hetatm;c;c=c->next) for(s1=c;s1;s1=s1->more) s1->pdb=this;
  }
}

Pdb *Pdb::clonepdb() {
 Pdb *t=new Pdb(this);
 if(next) t->next=next->clonepdb();
 else t->next=0;
 return t;
}

Pdb::~Pdb() 
{
 Strhandler cc;
 if(chn) {delete chn;chn=0;}
 if(next){delete next;next=0;}
 if(name){delete [] name;name=0;}
 if(token){delete [] token;token=0;}
 if(info) {
	info=cc.strdel(info);
 }
 if(endinfo) {
	endinfo=cc.strdel(endinfo);
 }
 if(hetatm) {
   delete hetatm;hetatm=0;
 }
}
 
Atom **Pdb::getAtomAll() {

Atom **aa=0;
Pdb *a=this;
int n=a->manyatm()*2+100;

n+=10;

aa=new Atom*[n];

for(int i=0;i<n;i++) { aa[i]=0; }

int nn=0;
for(Chn *cn=chn;cn;cn=cn->next)
for(Res *res=cn->res;res;res=res->next)
for(Atm *atm=res->atm;atm;atm=atm->next)
{
aa[nn]=new Atom(atm->xyz[0],atm->xyz[1],atm->xyz[2],atm->tatm->eng->radius);
aa[nn]->id=atm->id0;
nn++;
}

return aa;

}

Atm *Pdb::getatmbyoid(int n) {

	Chn *chn_temp;
	Atm *a=0;
  	for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next){
  		a=chn_temp->getatmbyoid(n);
		if(a) return a;
	}
	return 0;
}

void Pdb::write(char *s)
{
  FILE *fp;
  fp=fopen(s,"w");
  if(fp==0) {
        cerr<<"could not open file:" <<s<<endl;
        cerr<<"ignore the writing command..."<<endl;
        return;
  }
  write(fp);
  fflush(fp);
  fclose(fp);
}
void Pdb::writepdbused(char *s)
{
  FILE *fp;
  fp=fopen(s,"w");
  if(fp==0) {
        cerr<<"could not open file:" <<s<<endl;
        cerr<<"ignore the writing command..."<<endl;
        return;
  }
  writepdbused(fp);
  fflush(fp);
  fclose(fp);
}
void Pdb::writeold(char *s)
{
  FILE *fp;
  fp=fopen(s,"w");
  if(fp==0) {
        cerr<<"could not open file:" <<s<<endl;
        cerr<<"ignore the writing command..."<<endl;
        return;
  }
  int n=0;while(info&&info[n]) n++;
  int i;
  Strhandler cc;
  for(i=0;i<n;i++) {
        char *s=strdup(info[i]);
        s=cc.clearendchar(s,"\n"); 
        if(s==0) continue;
        fprintf(fp,"%s\n",s);
        if(s) delete [] s;s=0;
  }

  Chn *chn_temp;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next) {
  	chn_temp->writeold(fp);
	fprintf(fp,"TER\n");
  }
  for(chn_temp=hetatm;chn_temp;chn_temp=chn_temp->next) {
        chn_temp->writeold(fp);
        fprintf(fp,"TER\n");
  }

  n=0;while(endinfo&&endinfo[n]) n++;
  
  for(i=0;i<n;i++) {
        char *s=strdup(endinfo[i]);
        s=cc.clearendchar(s,"\n");
        if(s==0) continue;
        fprintf(fp,"%s\n",s);
        if(s) delete [] s;s=0;
  }


  fflush(fp);
  fclose(fp);
}
void Pdb::writerescard(char *s)
{
  FILE *fp;
  fp=fopen(s,"w");
  if(fp==0) {
        cerr<<"could not open file:" <<s<<endl;
        cerr<<"ignore the writing command..."<<endl;
        return;
  }
  int n=0;while(info&&info[n]) n++;
  int i;
  Strhandler cc;
  for(i=0;i<n;i++) {
        char *s=strdup(info[i]);
        s=cc.clearendchar(s,"\n"); 
        if(s==0) continue;
        fprintf(fp,"%s\n",s);
        if(s) delete [] s;s=0;
  }

  Chn *chn_temp;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next) {
  	chn_temp->writerescard(fp);
	fprintf(fp,"TER\n");
  }
  for(chn_temp=hetatm;chn_temp;chn_temp=chn_temp->next) {
        chn_temp->writerescard(fp);
        fprintf(fp,"TER\n");
  }

  n=0;while(endinfo&&endinfo[n]) n++;
  
  for(i=0;i<n;i++) {
        char *s=strdup(endinfo[i]);
        s=cc.clearendchar(s,"\n");
        if(s==0) continue;
        fprintf(fp,"%s\n",s);
        if(s) delete [] s;s=0;
  }


  fflush(fp);
  fclose(fp);
}


void Pdb::writeold(FILE *fp)
{
  if(fp==0)  return;
  int n=0;while(info&&info[n]) n++;
  int i;
  Strhandler cc;
  for(i=0;i<n;i++) {
        char *s=strdup(info[i]);
        s=cc.clearendchar(s,"\n"); 
        if(s==0) continue;
        fprintf(fp,"%s\n",s);
        if(s) delete [] s;s=0;
  }

  Chn *chn_temp;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
  {
  	chn_temp->writeold(fp);
	fprintf(fp,"TER\n");
  }
  for(chn_temp=hetatm;chn_temp;chn_temp=chn_temp->next)
  {
  	chn_temp->writeold(fp);
        fprintf(fp,"TER\n");
  }

  n=0;while(endinfo&&endinfo[n]) n++;
  for(i=0;i<n;i++) {
        char *s=strdup(endinfo[i]);
        s=cc.clearendchar(s,"\n");
        if(s==0) continue;
        fprintf(fp,"%s\n",s);
        if(s) delete [] s;s=0;
  }

}


void Pdb::setseqcard(){
  Chn *chn_temp;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
  chn_temp->setseqcard();
}

int Pdb::checkpdb(int flg){
  Chn *chn_temp;
  int nn=0;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
  nn+=chn_temp->checkpdb(flg);
  return nn;
}


void Pdb::setseqcardnogap(){
  Chn *chn_temp;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
  chn_temp->setseqcardnogap();
}


void Pdb::writeoldmore(char *s)
{
  FILE *fp;
  fp=fopen(s,"w");
  if(fp==0) {
	cerr<<"could not open file:" <<s<<endl;
	cerr<<"ignore the writing command..."<<endl;
	return;
  }
  Chn *chn_temp;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next) {
  	chn_temp->writeoldmore(fp);
	fprintf(fp,"TER\n");
  }
  for(chn_temp=hetatm;chn_temp;chn_temp=chn_temp->next) {
  	chn_temp->writeoldmore(fp);
	fprintf(fp,"TER\n");
  }

  fflush(fp);
  fclose(fp);
}

void Pdb::write(FILE *fp)
{
  int n=0;while(info&&info[n]) n++;
  int i;
  Strhandler cc;
  
  for(i=0;i<n;i++) {
	char *s=strdup(info[i]);
  	s=cc.clearendchar(s,"\n"); 
	if(s==0) continue;
	fprintf(fp,"%s\n",s);
	if(s) delete [] s;s=0;
  }
  
  Chn *chn_temp;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next) {
  	chn_temp->write(fp);
	fprintf(fp,"TER\n");
  }

  for(chn_temp=hetatm;chn_temp;chn_temp=chn_temp->next) {
        chn_temp->writeoldmore(fp);
        fprintf(fp,"TER\n");
  }
   
  n=0;while(endinfo&&endinfo[n]) n++;
  for(i=0;i<n;i++) {
        char *s=strdup(endinfo[i]);
        s=cc.clearendchar(s,"\n");
        if(s==0) continue;
        fprintf(fp,"%s\n",s);
        if(s) delete [] s;s=0;
  }
  
}

void Pdb::writepdbused(FILE *fp)
{
  int n=0;while(info&&info[n]) n++;
  //int i;
  Strhandler cc;
  /*
  for(i=0;i<n;i++) {
	char *s=strdup(info[i]);
  	s=cc.clearendchar(s,"\n"); 
	if(s==0) continue;
	fprintf(fp,"%s\n",s);
	if(s) delete [] s;s=0;
  }
  */
  Chn *chn_temp;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next) {
  	chn_temp->write(fp);
	fprintf(fp,"TER\n");
  }

  for(chn_temp=hetatm;chn_temp;chn_temp=chn_temp->next) {
        chn_temp->writeoldmore(fp);
        fprintf(fp,"TER\n");
  }
  /*  
  n=0;while(endinfo&&endinfo[n]) n++;
  for(i=0;i<n;i++) {
        char *s=strdup(endinfo[i]);
        s=cc.clearendchar(s,"\n");
        if(s==0) continue;
        fprintf(fp,"%s\n",s);
        if(s) delete [] s;s=0;
  }
  */
}

void Pdb::setoid(int n)
{
  Chn *chn_temp;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)chn_temp->setoid(n);
}


Chn *Pdb::operator[](int ch)
{
  Chn *chn_temp;
  int i;
  i=0;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
  {
    if(i==ch) return chn_temp; 
    i++;
  }
  return 0;
}

Chn*Pdb::ischain(int ch)
{
  Chn *chn_temp;
  int i;
  i=0;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
  {
    if(i==ch) return chn_temp;
    i++;
  }
  return 0;
}

Chn*Pdb::lastchain()
{
  Chn *chn_temp=0;
  if(chn==0) return 0;
  for(chn_temp=chn;chn_temp->next;chn_temp=chn_temp->next);
  return chn_temp;
}


Chn *Pdb::operator[](char ch)
{
 Chn *chn_temp;
 for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
 if(chn_temp->id==ch) return chn_temp;
 return 0;
}

Chn *Pdb::ischain(char ch)
{
 Chn *chn_temp;
 for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
 if(chn_temp->id==ch) return chn_temp;
 return 0;
}


int Pdb::getidnumber(Chn *s) {
  Chn *chn_temp;
  int i;
  i=0;
  for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
  {
    if(chn_temp==s) return i;
    i++;
  }
  return -1;
}
int Pdb::readmore(char *filnam)
{
FILE *fp;
int i,k,j,n;
char *rotamer,line[256];
Res *res_temp,*res_temp1,*res_temp2;
Tatm *tatm_temp;
Tres *tres_temp;
Atm *atm_temp;

if(filnam==0) return 0;
// find the file to open
if(name){delete [] name;name=0;}
i=strlen(filnam);
if(i==0) {cerr<<"error in file name";return 0;}
name=new char[i+1];
strcpy(name,filnam);
name[i]='\0';

//find directory
rotamer=RCS["library"];
fp=0;i=0;
//sprintf(line,"%s/%s",rotamer,filnam);
//i=strlen(rotamer);
sprintf(line,"%s",filnam);
i=-1;

while(((fp=fopen(line,"r"))==NULL))
{
 rotamer=rotamer+i+1;
 i=strlen(rotamer);
 if(*rotamer=='\0') 
 {
   cerr<<"warning:could not open: "<<filnam<<endl;
   cerr<<endl;
   cerr<<"check the jackal.dir file first in your current directory then in the environment path:JACKALDIR"<<endl;
   cerr<<"something wrong in specifying the library"<<endl;
   cerr<<"make sure the library points to the abosulte path"<<endl;
   exit(0);
 }
 sprintf(line,"%s/%s",rotamer,filnam);
}

cerr<<endl;
cerr<<"reading from file: "<<line<<endl;
cerr<<endl; 

if(fp==0) {
  cerr<<"warning:could not open: "<<filnam<<endl;
  cerr<<endl;
  cerr<<"check the jackal.dir file first in your current directory then in the environment path:JACKALDIR"<<endl;
  cerr<<"something wrong in specifying the library"<<endl;
  cerr<<"make sure the library points to the abosulte path"<<endl;

  exit(0);
}

res_temp=0;
atm_temp=0;
tatm_temp=0;
tres_temp=0;
k=1;

// read pdb line by line
chn=new Chn;
chn->pdb=this;

while(fgets(line,256,fp)!=NULL) 
{
  if(strlen(line)<30) continue; 
  if(line[0]=='!')continue;

  
  if(strncmp(line,"HEAD",4)==0) {
    i=atoi(line+4)-1;
    if(i==0) {res_temp->temp=new float[9];res_temp->nemp=9;}
    for(j=0;j<3;j++) res_temp->temp[i*3+j]=atof(line+8+8*j);
    if(i==2)k=1;
    continue;
  }
  else if(strncmp(line,"STAT",4)==0) {
    res_temp->rotdegree=atoi(line+8);
    res_temp->rotenergy=atof(line+16);
    res_temp->rotchance=atof(line+24);
    k=1;
    continue;
  }
  else 
  {
    if(k==1) tres_temp=TRES[line[5]];
    if(tres_temp==0) continue;
    if(res_temp==0&&chn->res==0)
    {
      chn->res=new Res;
      res_temp=chn->res;
      res_temp->tres=tres_temp;
      res_temp->name=tres_temp->name;
      k=0;
    }
    else if(k==1)
    {
      res_temp->next=new Res;
      res_temp=res_temp->next;
      res_temp->tres=tres_temp;
      res_temp->name=tres_temp->name;
      k=0;
    }
    tatm_temp=tres_temp->isatm(line);
    if(tatm_temp==0) continue;
    if(res_temp->isatm(line)) continue;
    if(res_temp->atm==0) { res_temp->atm=new Atm; atm_temp=res_temp->atm;}
    else {atm_temp->next=new Atm; atm_temp=atm_temp->next;}
    for(j=0;j<3;j++) atm_temp->xyz[j]=atof(line+8+8*j);
    atm_temp->tatm=tatm_temp;
    atm_temp->res=res_temp;
    strcpy(atm_temp->name,tatm_temp->name);
  }
}

//reorgnize the residue sequence!
//if(name) {delete [] name;name=0;}
//name=new char[100];
//strcpy(name,"rotamer");
fclose(fp);

Res **rr[30];

n=0;
for(res_temp=chn->res;res_temp;res_temp=res_temp->next)n++;

for(k=0;k<30;k++) 
{
  rr[k]=new Res*[n];
  for(i=0;i<n;i++) rr[k][i]=0; 
}

k=0;
res_temp1=chn->res;
for(res_temp=chn->res;res_temp;res_temp=res_temp->next)
{
  if(res_temp1->name!=res_temp->name) 
  {
   cerr<<"from "<<filnam<<": the number copies for residue: "<<res_temp1->name<<" "<<k<<endl;
   k=0;
  }
  rr[res_temp->name-'A'][k]=res_temp;
  k++;
  res_temp1=res_temp;
}

cerr<<"from "<<filnam<<": the number copies for residue: "<<res_temp1->name<<" "<<k<<endl;

chn->res=0;
chn->create(&TRES);
chn->configure();
 
for(res_temp=chn->res;res_temp;res_temp=res_temp->next) 
for(Atm *aa1=res_temp->atm;aa1;aa1=aa1->next) {
   if(aa1->tatm->id==0) aa1->bond[0]=0;
   else if(aa1->tatm->id==2) aa1->bond[1]=0;
}

for(k=0;k<30;k++) for(i=0;i<n;i++)
if(rr[k][i]) {rr[k][i]->more=0;rr[k][i]->next=0;}

//put the residue on more list
for(res_temp=chn->res;res_temp;res_temp=res_temp->next)
{
   for(k=0;k<n;k++)
   {
     res_temp2=rr[res_temp->name-'A'][k];
     if(res_temp2==0) break;
     if(k==0) {res_temp->more=res_temp2;res_temp1=res_temp->more;}
     else     {res_temp1->more=res_temp2;res_temp1=res_temp1->more;}
   }
}

//bonded, not including connectivity  

Atm *aa1,*aa2;
Atm  *bb1,*bb2;
for(res_temp=chn->res;res_temp;res_temp=res_temp->next)
{
  for(res_temp1=res_temp->more;res_temp1;res_temp1=res_temp1->more)
  { 
    for(aa1=res_temp->atm;aa1;aa1=aa1->next)
    {
       aa2=(*res_temp1)[aa1->tatm->id];
       if(aa2==0) continue;
       for(k=0;k<aa1->tatm->nbond;k++)
       {
         bb1=aa1->bond[k];
         if(bb1==0) continue;
         if(bb1->res->id==res_temp->id)
         {
           bb2=(*res_temp1)[bb1->tatm->id];
           aa2->bond[k]=bb2;
         }
         else
         {
           aa2->bond[k]=bb1;
         }
       }
    }
  }
}

for(res_temp=chn->res;res_temp;res_temp=res_temp->next)
{
  res_temp->nummore=res_temp->many();
}

for(i=0;i<30;i++) delete [] rr[i];
return 1;
}

int Pdb::isresidchar(char *line) {

	if(line==0) return 0;
	int n=strlen(line);
	int i;
	for(i=22;i<=26&&i<n;i++) {
		if(line[i]==' '||line[i]=='-') continue;
		if(line[i]>='0'&&line[i]<='9') continue;		 
		return 1;
	}
	return 0;
}

int Pdb::read(char *filnam,char chi)
{
Strhandler cc;
FILE *fp,*fp1;
int i,j;
char *pdb,line[256],line1[256],resname;
Chn *chn_temp;
Res *res_temp;
Tatm *tatm_temp;
Tres *tres_temp;
Atm *atm_temp;
Strhandler ceg;
if(filnam==0) return 0;
// find the file to open
if(name){delete [] name;name=0;}
i=strlen(filnam);
if(i==0) {cerr<<"error in file name";return 0;}
name=new char[i+1];
strcpy(name,filnam);
name[i]='\0';

char *ss=cc.lastindexof(name,"/");
if(ss&&ss[0]!='\0') {
   char *tt=strdup(ss);
   delete [] name;
   name=strdup(tt);
   delete [] tt;
   if(ceg.strncmplast(filnam,".pdb",4)==0) {
	ss=ceg.strstrlast(name,".pdb");    
   	if(ss) *ss='\0';
   }
}

pdb=RCS["pdb"];
fp=0;i=0; fp1=0;
//sprintf(line,"%s/%s",pdb,filnam);
//sprintf(line1,"%s/%s.pdb",pdb,filnam);
//i=strlen(pdb);
sprintf(line,"%s",filnam);

if(ceg.strncmplast(filnam,".pdb",4)==1) {
	sprintf(line1,"%s.pdb",filnam);
}
else {
	sprintf(line1,"%s",filnam);
}

i=-1;
while(((fp=fopen(line,"r"))==NULL)&&((fp1=fopen(line1,"r"))==NULL))
{
 pdb=pdb+i+1;
 if(pdb==0) {
	cerr<<"warning:could not open: "<<filnam<<endl;
	return 0;
 }
 i=strlen(pdb);
 if(*pdb=='\0') 
 {cerr<<"warning:could not open: "<<filnam<<endl;return 0;}
 sprintf(line,"%s/%s",pdb,filnam);
 sprintf(line1,"%s/%s.pdb",pdb,filnam);
}
if(fp1) {
	fp=fp1;
	cerr<<endl;
	cerr<<"reading pdb file:" <<line1<<endl;
	cerr<<endl;
}
else {
	cerr<<endl;		
	cerr<<"reading pdb file:" <<line<<endl;
	cerr<<endl;
}
chn_temp=0;
res_temp=0;
atm_temp=0;
tatm_temp=0;
tres_temp=0;
info=new char*[2000];
for(i=0;i<2000;i++) info[i]=0;

endinfo=new char*[2000];
for(i=0;i<2000;i++) endinfo[i]=0;

// read pdb line by line
int ninfo=0;int endninfo=0;
int startAtm=0;
while(fgets(line,256,fp)!=NULL) {
  if(strncmp(line,"TER",3)==0) {//new chain starts
	if(chn==0) //first chain
  	{
    		chn=new Chn;
		chn_temp=chn; 
  	}
  	else //the next new chain
  	{
    		chn_temp->next=new Chn; 
    		chn_temp=chn_temp->next;
  	} 
	continue;
  }
  if(strncmp(line,"ATOM  ",6)!=0&&strncmp(line,"HETATM",6)!=0&&startAtm==0) {
      if(isinfo==0) continue;
      if(ninfo>1990) continue;
      info[ninfo]=strdup(line);
      ninfo++;
      continue;
  }
  else if(strncmp(line,"ATOM  ",6)!=0&&startAtm==1&&strncmp(line,"HETATM",6)!=0) {
      if(isinfo==0) continue;
      if(endninfo>1990) continue;
      endinfo[endninfo]=strdup(line);
      endninfo++;
      continue;
  }
  else if(strncmp(line,"ATOM  ",6)!=0&&startAtm==1&&ishet==0) {
      cerr<<"ignored: "<<line;
      if(isinfo==0) continue;
      if(endninfo>1990) continue;
      endinfo[endninfo]=strdup(line);
      endninfo++;
      continue;
  }

  if(strncmp(line,"HETATM",6)==0&&ishet==0) {
      cerr<<"ignored: "<<line;
      continue; //later change this
  }
  

  if(strncmp(line,"ATOM  ",6)!=0&&strncmp(line,"HETATM",6)!=0) {
        continue;  
  }
 
  //is the chain wanted?
  startAtm=1;

 

  if(chi=='1')  chi=line[21];

  if(chi!='0')  {

	if(chi!=line[21]&&strncmp(line,"ATOM  ",6)==0) {
		cerr<<"ignroed! not the right chain:"<<line<<endl;
		continue;
	}
        else if(chi!=line[21]&&strncmp(line,"HETATM",6)==0&&line[21]!=' ') {
		cerr<<"ignroed! not the right chain:"<<line<<endl;
		continue;
        }
  }
  
//end 
  
//process the hetatm.
  if(strncmp(line,"HETATM",6)==0) {
	parsehetatmline(line);
        continue;
  }

//is the residue standard?
  tres_temp=TRES[line+17];
  if(tres_temp==0&&ishet) 
  {
    cerr<< "not standard residue, treated as hetatm:  "<<line<<endl;
    parsehetatmline(line);
    continue;
  }
  else if(tres_temp==0) {
    cerr<< "not standard residue:  "<<line<<endl;    
    continue;
  }

  resname=tres_temp->name;
  if(res_temp)
  { 
    //if(atoi(line+22)==res_temp->id&&resname!=res_temp->name) {//confusing aa 
    if(strncmp(line+22,res_temp->rescard,5)==0&&resname!=res_temp->name) {//confusing aa 
	cerr<<"ignored: "<<line<<endl;
    	continue; 
    }  
  }

//end

//convert charmm to pdb term
  if(strncmp(line+12," OT",3)==0)
  {
    line[14]=' ';
    line[15]=' ';
  } 
  else if(strncmp(line+12," HT",3)==0) 
  { 
    line[12]=' ';//line[15];
    line[15]=' ';
    line[14]='N';
  }
  else if(line[12]=='H')
  {
    strncpy(line1,line+12,4);
    line[12]=line1[3];
    line[13]=line1[0];
    line[14]=line1[1];
    line[15]=line1[2];
  }
  else if(line[13]=='H'&&line[12]==' ')
  {
    i=tres_temp->getuse(line+12,2,2);
    if(i==0)
    {
      i=tres_temp->getuse(line+12,2,1);
      if(i==0) {
	   cerr<<"ignored: "<<line<<endl;
	   continue;
      }
      else {line[12]=line[15];line[15]=' ';}
    }
  }
//end

//treat the last oxt atoms

  if(strncmp(line+13,"OXT",3)==0) {
	cerr<<"ignored: "<<line<<endl;
	continue;
  }

//end

//find if the atom belong to standard atoms of this residue

  if(strncmp(line+13,"H  ",3)==0){line[12]=' ';line[14]='N';}
  //if(res_temp&&atoi(line+22)==res_temp->id)//the same residue
  if(res_temp&&strncmp(line+22,res_temp->rescard,5)==0)//the same residue
  {
    i=res_temp->getuse(line+12,1,3);
    j=tres_temp->getuse(line+12,1,3);
    if(i>0) {
	if(line[16]!=' ') res_temp->addambgt(line);
    }
    if(i==j&&j!=0){
	/*res_temp->repeat=1;*/
	//res_temp->addambgt(line);
	cerr<<"ignored: "<<line<<endl;
	continue;
    } //no space for this atom
    if(line[13]!='H'&&i>0) {
	/*res_temp->repeat=1;*/
	cerr<<"ignored: "<<line<<endl;
	continue;
    } //this atom already exist!
    if(res_temp->isatm(line+12)) {
	/*res_temp->repeat=1;*/
	cerr<<"ignored: "<<line<<endl;
	continue;
    }// exist!
    tatm_temp=tres_temp->isatm(line+12,1,3,i);
    if(tatm_temp==0) 
    {
      i=res_temp->getuse(line+12,1,2);
      j=tres_temp->getuse(line+12,1,2); 
      if(i==j&&j!=0) {
	/*res_temp->repeat=1;*/
	cerr<<"ignored: "<<line<<endl;
	continue;
      } //no space for this atom
      if(j!=1) {
	cerr<<"ignored: "<<line<<endl;
	continue; //two atoms or none match this atom
      }
      tatm_temp=tres_temp->isatm(line+12,1,2,i);         
      if(tatm_temp==0) {
	cerr<<"ignored: "<<line<<endl;
	continue;
      }
    }
    if((*res_temp)[tatm_temp->id]) {
	/*res_temp->repeat=1;*/
	cerr<<"ignored: "<<line<<endl;
	continue;
    }
  }
  else
  {
    tatm_temp=tres_temp->isatm(line+12,1,3,0); 
    if(tatm_temp==0) 
    {
      i=tres_temp->getuse(line+12,1,2);
      if(i!=1) {
	cerr<<"ignored: "<<line<<endl;
	continue; //two or none such kind of atoms
      }
      tatm_temp=tres_temp->isatm(line+12,1,2,0);
      if(tatm_temp==0) {
	cerr<<"not standard atom. ignored: "<<line<<endl;
	continue; //not standard atom found!
      }
    }
  }

  if(tatm_temp==0) {
	cerr<<"not standard atoms,ignored: "<<line<<endl;
	continue;  // not standard atoms
  }
 
//end
  
  if(chn==0) //first chain
  {
    chn=new Chn;chn_temp=chn; 
  }
  else if(line[21]!=chn_temp->id)//the next new chain
  {
    chn_temp->next=new Chn; 
    chn_temp=chn_temp->next;
  } 
  else if(res_temp&&res_temp->id>atoi(line+22)&&isresidchar(line)==0)// the next line residue id is smaller
  {
    chn_temp->next=new Chn;
    chn_temp=chn_temp->next;
  }
  chn_temp->pdb=this;
  if(chn_temp->res==0)// the chain id not assigned 
  { 
    chn_temp->id=line[21];
    chn_temp->res=new Res; // the first residue of the chain
    res_temp=chn_temp->res;  
    chn_temp->pdb=this;
  }
  //else if( res_temp->id!=atoi(line+22)) // the next new residue
  else if(res_temp&&strncmp(line+22,res_temp->rescard,5)!=0)
  {
    res_temp->next=new Res;
    res_temp=res_temp->next;
  }
  if(res_temp->atm==0) // the residue id not assigned
  { 
    res_temp->id=atoi(line+22);
    strncpy(res_temp->rescard,line+22,5); //for resid in card!
    res_temp->rescard[5]='\0';
    res_temp->oid=res_temp->id;
    res_temp->name=resname;
    res_temp->tres=tres_temp;
    res_temp->chn=chn_temp;
    res_temp->atm=new Atm; // the first atom of the residue
    atm_temp=res_temp->atm; 
  } 
  else if(res_temp->isatm(line+13)==0) // the next new atom
  {
    atm_temp->next=new Atm;
    atm_temp=atm_temp->next;
  }
  else {
	/*res_temp->repeat=1;*/
	cerr<<"redundant atoms.ignored: "<<line<<endl;
	continue;
  } // the same atom
  if(atm_temp->id0==0 ) // the new atom
  {  
    atm_temp->id0=atoi(line+6);
    atm_temp->oid=atm_temp->id0;
    atm_temp->res=res_temp;    
  }
  for(i=0;i<3;i++) atm_temp->xyz[i]=atof(line+30+8*i);
  atm_temp->occup=atof(line+54);
  atm_temp->bfact=atof(line+60);
  strncpy(atm_temp->name,tatm_temp->name,4);//strncpy(atm_temp->name,line+12,4);
  atm_temp->name[4]='\0';
  atm_temp->tatm=tatm_temp;
}
  fclose(fp);
  cleanemptychn();
  indexreswithcard();
  configure();
  checkpdbchainid();
  configure();
  //checkhetatmchainid(); 
  hetatmconfigure();
  cerr<<endl<<"set up force field parameters for hetero atoms..." <<endl;
  TRES.setenghetatm("hetatm.prm");
  
  return 1;
}

void Pdb::cleanemptychn() {

	

}
int Pdb::isconnected(Res *a,Res *b) {

	if(a==0||b==0) return 0;
	Atm *aa[4];
	Atm *bb[4];	

	//first
	aa[0]=a->isatm(" C  ");
  	aa[1]=a->isatm(" O  ");
	aa[2]=a->isatm(" CA ");
	aa[3]=a->isatm(" N  ");
	//next
	bb[0]=b->isatm(" C  ");
  	bb[1]=b->isatm(" O  ");
	bb[2]=b->isatm(" CA ");
	bb[3]=b->isatm(" N  ");
	
	int i,j;
	int ii,jj;
	ii=0;
	jj=0;
	for(i=0;i<4;i++) { 
	  	j=i;
		if(aa[j]==0||bb[i]==0) continue;
		float d=TRES.distance(bb[i],aa[j]);	
		ii++;	 		
		//N with the previous residue
		if(d<1.0) jj++;
	} 
	if(ii==jj) return -1; //identical residues

	//first
	aa[0]=a->isatm(" C  ");
  	aa[1]=a->isatm(" O  ");
	aa[2]=a->isatm(" CA ");
	
	//next
	bb[0]=b->isatm(" N  ");
	bb[1]=b->isatm(" CA ");
	ii=0;jj=0;
	for(i=0;i<2;i++)  
	for(j=0;j<3;j++) { 	
		if(aa[j]==0||bb[i]==0) continue;
		float d=TRES.distance(bb[i],aa[j]);
		ii++;		 
		//N with the previous residue
		if(i==0&&j==0&&d<2.3&&d>0.7)  jj++;
		else if(i==0&&j==1&&d<3.3&&d>1.2) jj++;
		else if(i==0&&j==2&&d<3.5&&d>1.2) jj++;
		//CA with the previous residue
		else if(i==1&&j==0&&d<3.5&&d>1.2) jj++;
		else if(i==1&&j==1&&d<3.7&&d>1.3) jj++;
		else if(i==1&&j==2&&d<4.8&&d>1.9) jj++;

		
	} 

	if(ii==jj) return 1;

	return 0;
}

void Pdb::indexreswithcard() {	
	Chn *c;
	for(c=chn;c;c=c->next) {	
		indexreswithcard(c); 
	}
}

void Pdb::adaptindex(Res *r,int n) {

	Res *s;
	if(r==0) return;

	if(r->next==0) {
		r->oid=n;
		r->id=n;
		return;
	}

	if(r->rescard[4]!=' '||r->next->rescard[4]!=' ') {
		r->oid=n;
		r->id=n;
		return; 
	}
	for(s=r->next;s;s=s->next) {
		s->oid=s->oid-r->oid+n;
		s->id=s->oid;	
	}
	r->oid=n;
	r->id=n;
}

void Pdb::indexreswithcard(Chn *c) {
  	if(c==0) return;
	Res *r;
	for(r=c->res;r;r=r->next) {
		if(r->next==0) continue;
		if(r->next->id>r->id) continue;
		int nn=isconnected(r,r->next);
		if(nn==1) {
			adaptindex(r->next,r->id+1);
			if(TRES.logg>3)cerr<<"reindex the residue id: "<<r->next->name<<r->next->oid<<r->next->rescard[4]; 
			//r->next->id=r->id+1;
			//r->next->oid=r->next->id;
			if(TRES.logg>3)cerr<<"-->"<<r->next->name<<r->next->oid<<endl;
		}
		else if(nn==0) {
			adaptindex(r->next,r->id+2);
			if(TRES.logg>3)cerr<<"reindex the residue id with chain break: "<<r->next->name<<r->next->oid<<r->next->rescard[4]; 
			//r->next->id=r->id+2;
			//r->next->oid=r->next->id;
			if(TRES.logg>3)cerr<<"-->"<<r->next->name<<r->next->oid<<endl;
		} 		
		else if(nn==-1) {//identical
			if(TRES.logg>3)cerr<<"the two residues at the same place and should be the same residues, ignore the second: ";
			if(TRES.logg>3)cerr<<r->name<<r->oid<<"--"<<r->next->name<<r->next->oid<<r->next->rescard[4]<<endl;
			Res *rr=r->next;
			r->next=rr->next;
			rr->next=0;
			delete rr;			
		}
	}
	//writerescard("ses");
}

void Pdb::setrescard(int f) {

	Chn *c;

	for(c=chn;c;c=c->next) {
		Res *r;
		for(r=c->res;r;r=r->next) {
			Res *s=r;			 
			if(f&&strlen(s->rescard)==5) cerr<<"indexing from old to new..."<<s->rescard<<" "<<c->id<<"--";
			else if(f) cerr<<"indexing from old to new..."<<"....."<<" "<<c->id<<"--";
			sprintf(s->rescard,"%4d ",s->oid);
			s->rescard[5]='\0';
			if(f)cerr<<s->rescard<<endl; 
		}
	}
}

int Pdb::checkrescard() {

	Chn *c;

        for(c=chn;c;c=c->next) {
                Res *r;
                for(r=c->res;r;r=r->next) {
			if(strlen(r->rescard)<5) return 1;
                }
        }
	return 0;
}

void Pdb::parsehetatmline(char *line) {
  Res *r;
  Chn *c=0;
  int rid=atoi(line+22); 
  //find the chain with the same id
  for(c=hetatm;c;c=c->next) {  
	if(c->id==line[21]) {
		break;
	}
  }
  
  //the first chain
  if(c==0&&hetatm==0) {
	hetatm=new Chn;
	c=hetatm;
  }
  else if(c==0) { //new chain
	for(c=hetatm;c->next;c=c->next);
	c->next=new Chn;
	c=c->next;
  }
  else if(c) {//the old chain,detect residue id jumps.
	for(r=c->res;r&&r->next;r=r->next);
	if(r&&rid<r->oid) {//new chain since the residue seqid jumps to smaller
		c->next=new Chn;
		c=c->next;
	}
  }
  c->pdb=this;
  c->id=line[21];
  c->ishet=1;
  //find residue
  
  
  for(r=c->res;r;r=r->next) {
	if(r->oid==rid&&strncmp(r->tres->name3,line+17,3)==0)  {
  		break;
   	}
  }
  
  //the first residue
  if(r==0&&c->res==0) {
	c->res=new Res;
	r=c->res;
  }
  else if(r==0) { //the new residue
	for(r=c->res;r->next;r=r->next);
	r->next=new Res;
	r=r->next;
  }
  r->chn=c;
  r->id=rid;
  r->oid=rid;
  r->name='-';

  Tres *t;
  //create the tres if the new res type
  if(TRES.hetm==0&&r->tres==0) {
	TRES.hetm=new Tres;
	t=TRES.hetm;
	strncpy(t->name3,line+17,3);
	t->name3[3]='\0';
  }
  else if(r->tres==0) {
	t=TRES.hetm->addhetm(line+17);
  }
  else if(r->tres) {
	t=r->tres;
  }
  r->tres=t;
  t->name='-';
  t->id=rid;
  //detect atoms if old or new
  Atm *a;

  for(a=r->atm;a;a=a->next) {
	if(strncmp(a->name,line+13,4)==0) {
		cerr<<"ignore, redundant! "<<line<<endl;
		return;
	}
  }

  Tatm *ta=t->addhetmatm(line+12);

  if(r->atm==0) {
	r->atm=new Atm;
	a=r->atm;
  }
  else {
  	for(a=r->atm;a->next;a=a->next);
  	a->next=new Atm;
  	a=a->next;
  }
  a->res=r;
  a->id0=atoi(line+6);
  a->oid=a->id0;
  a->tatm=ta;  

  int i;
  for(i=0;i<3;i++) a->xyz[i]=atof(line+30+8*i);
  a->occup=atof(line+54);
  a->bfact=atof(line+60);
  strncpy(a->name,ta->name,4); 
  a->name[4]='\0';
  a->tatm=ta;
}

void Pdb::checkpdbchainid() {

  //Res *r;
  Chn *c,*c1;

  while(1) {
	int n=0;
  	for(c=chn;c;c=c->next) {
		if(c->res==0) {
			cerr<<endl;
			cerr<<"remove empty chain.."<<c->id<<endl;
			cerr<<endl;
			removechain(c);
			n++;
		}
  	}
	if(n==0) break;
  }

  cerr<<endl;
  
  while(1)  {
  int ng=0;
  for(c=chn;c;c=c->next) {
  	for(c1=c->next;c1;c1=c1->next) {
		if(c->id==c1->id&&c->res&&c1->res) {
			cerr<<"there are two chains having the same id:"<<c->id<<" only keep the first one!"<< endl;
			if(delchn==1) removechain(c1);
			else relinkchain(c,c1);
			ng++;
			goto re200;
			/*
			cerr<<"the possible reason is that:"<<endl;
			Res *r0=c->lastres();
			cerr<<c1->res->name<<c1->res->oid<<" must have larger id than "<<r0->name<<r0->oid<<endl;
			exit(0);
			*/
		}
	}         	
  }
  re200:
  if(ng==0) break;
  }
}


void Pdb::checkhetatmchainid() {

  //Res *r;
  Chn *c,*c1;

  for(c=hetatm;c;c=c->next) {
	if(c->res==0) {
		cerr<<endl;
		cerr<<"remove empty hetatm chain.."<<c->id<<endl;
		cerr<<endl;
		removehetatmchain(c);
	}
  }

  cerr<<endl;
  
  while(1)  {
  int ng=0;
  for(c=hetatm;c;c=c->next) {
  	for(c1=c->next;c1;c1=c1->next) {
		if(c->id==c1->id&&c->res&&c1->res) {
			cerr<<"there are two hetatm chains having the same id:"<<c->id<<" only keep the first one!"<< endl;
			if(delchn==1) removehetatmchain(c1);
			else relinkhetatmchain(c,c1);
			ng++;
			goto re200;
			/*
			cerr<<"the possible reason is that:"<<endl;
			Res *r0=c->lastres();
			cerr<<c1->res->name<<c1->res->oid<<" must have larger id than "<<r0->name<<r0->oid<<endl;
			exit(0);
			*/
		}
	}         	
  }
  re200:
  if(ng==0) break;
  }
}

void Pdb::configure() {
 
 Chn *chn_temp;
 Res *res_temp;
 Atm *atm_temp;
 Tatm *tatm_temp;

 number=0;
 for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
 {
   chn_temp->pdb=this;
   if(chn_temp->res==0) continue;
   //chn_temp->start=chn_temp->res->id;
   //for(res_temp=chn_temp->res;res_temp;res_temp=res_temp->next) res_temp->id-=chn_temp->start;
   chn_temp->addresid(chn_temp->start);
   chn_temp->start=0;
   chn_temp->configure();
   number++;
 }

//detect the number of residue and atoms
int n=0,nn=0,rr=0;
for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
for(res_temp=chn_temp->res;res_temp;res_temp=res_temp->next)
{
  res_temp->id0=rr++;
  //atom id
  for(atm_temp=res_temp->atm;atm_temp;atm_temp=atm_temp->next)
  {
    atm_temp->id0=n++;
  }

  //atom id0 with tres complete supposed
  for(tatm_temp=res_temp->tres->tatm;tatm_temp;tatm_temp=tatm_temp->next)
  {
    atm_temp=(*res_temp)[tatm_temp->id];
    if(atm_temp)atm_temp->id=nn;
    nn++;
  }
}
 
}


void Pdb::hetatmconfigure() {
 
 Chn *chn_temp;
 Res *res_temp;
 Atm *atm_temp;
  

 int nn=0,rr=0;
 //int rrid=0;
 for(chn_temp=chn;chn_temp;chn_temp=chn_temp->next)
 for(res_temp=chn_temp->res;res_temp;res_temp=res_temp->next) {
	rr=res_temp->id0;//rrid=res_temp->id;
 	for(atm_temp=res_temp->atm;atm_temp;atm_temp=atm_temp->next) {
		nn=atm_temp->id0;
 	}
 }
 rr++;nn++;//rrid++;
 
 for(chn_temp=hetatm;chn_temp;chn_temp=chn_temp->next) {
	if(chn_temp->res==0) continue;
	//chn_temp->start=chn_temp->res->oid-rrid;	
	chn_temp->start=chn_temp->res->oid;
 	for(res_temp=chn_temp->res;res_temp;res_temp=res_temp->next) {
		res_temp->id0=rr++;res_temp->id=res_temp->oid-chn_temp->start;
		//rrid=res_temp->id;
 		for(atm_temp=res_temp->atm;atm_temp;atm_temp=atm_temp->next)
 		{
   			atm_temp->id0=nn++;atm_temp->id=atm_temp->id0;
 		}
 	}
	//rrid++;
 }
}



void Pdb::convrt(Chn *s,FILE *fp) 
//convert the pdb file format to charmm format
{
Res *r;
Atm *a,*a1;
int i,m,j,k,nm,p;
nm=0;
for(r=s->res;r;r=r->next)
{
for(a=r->atm;a;a=a->next)
{
  if(a->tatm->name[1]=='H') continue;
  if(r->name=='I'&&a->tatm->id==7)a->name[3]=' ';
  j=0;k=0;
  for(i=0;i<4;i++) if(a->name[i]!=' ')j++;
  for(i=0;i<4;i++) if(a->bond[i]&&a->bond[i]->tatm->name[1]=='H') k++;
  if(k==0) continue;
  m=0;
  for(i=0;i<4;i++)
  {    
    a1=a->bond[i];
    if(a1==0) continue;
    if(a1->tatm->name[1]!='H')continue;
    m++;
    if(j==1)  
    {
      if(a->tatm->name[1]=='N')
      {
        sprintf(a1->name," HN \0");break;
      }
      else
      {
       cerr<<"confusing atom with hydrogen:"<<a->name<<endl;
      }
    }
    else if(j==2)
    {
      if(k>1) sprintf(a1->name," H%c%i",a->tatm->name[2],m);
      else sprintf(a1->name," H%c ",a->tatm->name[2]);
    }
    else if(j==3)
    {
      if(k>1) sprintf(a1->name,"H%c%c%i",a->tatm->name[2],a->tatm->name[3],m);
      else sprintf(a1->name," H%c%c",a->tatm->name[2],a->tatm->name[3]);
    }

  }//i==4;
}//a
if(fp)
{
 for(a=r->atm;a;a=a->next)
 {
  if(a->tatm->name[1]=='H') continue;
  nm++;p=a->id0;a->id0=nm;a->write(fp);a->id0=p;
  for(i=0;i<4;i++)
  {
   a1=a->bond[i];   
   if(a1==0) continue;
   if(a1->tatm->name[1]=='H') 
   {nm++;p=a->id0;a1->id0=nm;a1->write(fp);a->id0=p;}
  }
 }
}
}//r
}

int  Pdb::dssp0(char *filnam0)
{
FILE *fp;
char  line[256],*dssp;
int k;
char ch;
Chn *chn_temp;
Res *r;
dssp=RCS["dssp"];
sprintf(line,"%s/%s.dssp",dssp,filnam0);
fp=fopen(line,"r");
if(fp==0) return 0;

while(fgets(line,256,fp)!=NULL)
{
  if(strstr(line,"#  RESIDUE AA STRUCTURE")!=NULL) break;
}

ch='?';
while(fgets(line,256,fp)!=NULL)
{
  if(strchr(line,'!')) continue;
  k=atoi(line+5);
  if(ch!=line[11])
  {
    ch=line[11];
    chn_temp=(*this)[ch];
    if(chn_temp==0) {ch='?';continue;}
    for(r=chn_temp->res;r;r=r->next)r->sec='?';
  }
  k=k-chn_temp->start;
  r=chn->isres(k);
  if(r) 
  {
    r->sec=line[16];
    if(r->sec!=' ')r->sec+=32;
    continue;
  }
}

fclose(fp);

return 1;

}
void Pdb::dssp(char *filnam0)
{
FILE *fp,*fp1;
char  line1[256],line[256],*dssp,*bin;
int k,i;
char ch;
char filnam[1000],filnam1[1000];
char  *pdb,*pdbfile;
Chn *chn_temp;
Res *r;

//locate pdb file directory
pdb=RCS["pdb"];
fp=0;i=0; fp1=0;
sprintf(line,"%s/%s",pdb,filnam0);
sprintf(line1,"%s/%s.pdb",pdb,filnam0);
i=strlen(pdb);
while(((fp=fopen(line,"r"))==NULL)&&((fp1=fopen(line1,"r"))==NULL))
{
 pdb=pdb+i+1;
 i=strlen(pdb);
 if(*pdb=='\0')
 {cerr<<"warning:could not open: "<<filnam<<endl;return ;}
 sprintf(line,"%s/%s",pdb,filnam0);
 sprintf(line1,"%s/%s.pdb",pdb,filnam0);
}
if(fp) {pdbfile=line;fclose(fp);}
else if(fp1) {pdbfile=line1;fclose(fp1);}

strcpy(filnam,pdbfile);

//locate dssp directory

dssp=RCS["dssp"];
bin=RCS["bin"];
sprintf(filnam1,"%s/%s.dssp",dssp,filnam0);
sprintf(line,"%s/dssp -na %s %s",bin,filnam,filnam1); 
system(line); // calculating secondary structure
fflush(stdout);
if((fp=fopen(filnam1,"r"))==NULL) 
{
  cerr<<"could not open:"<<filnam1<<endl;
  return;
}

while(fgets(line,256,fp)!=NULL)
{
  if(strstr(line,"#  RESIDUE AA STRUCTURE")!=NULL) break;
}
 
ch='?';
while(fgets(line,256,fp)!=NULL)
{
  if(strchr(line,'!')) continue;
  k=atoi(line+5);
  if(ch!=line[11]) 
  {
    ch=line[11];
    chn_temp=(*this)[ch]; 
    if(chn_temp==0) {ch='?';continue;}
    for(r=chn_temp->res;r;r=r->next)r->sec='?';
  }
  k=k-chn_temp->start;
  r=chn->isres(k);
  if(r)  
  {
    r->sec=line[16];
    if(r->sec!=' ')r->sec+=32;
    continue;
  } 
}

fclose(fp);
}

void Pdb::setoff(float *xyz,float cut)
{
  Res *r;
  Atm *a;
  float d;
  for(Chn *ch=chn;ch;ch=ch->next)
  for(r=ch->res;r;r=r->next)
  for(a=r->atm;a;a=a->next)
  {
    d=TRES.distance(a->xyz,xyz);
    if(d<cut) a->flag=1;
  }
}

void Pdb::setflgr(int f) {

  for(Chn *ch=chn;ch;ch=ch->next)
        ch->setflgr(ch->res,1000000,f);

}

void Pdb::setflgr(char c,int f) {

  for(Chn *ch=chn;ch;ch=ch->next)
        ch->setflgr(c,f);

}

void Pdb::setflg(int f) {
  
  for(Chn *ch=chn;ch;ch=ch->next) 
 	ch->setflg(ch->res,1000000,f);

}
void Pdb::setoff(Res *r0,int n,float cut)
{
 Res *r;
 Atm *a;
 Lattice lat;
 int i;
 float d;

 lat.putoff();
 lat.grdsiz=2.0;
 lat.ready(this);
 lat.grdsiz=2;
 lat.radall=15;
 lat.flag=0;
 lat.puton(this);
 lat.putoff(r0,n);
 setflg(0);
 //setflg(res,1000000,0);
 r0->chn->setflg(r0,n,1);

 for(r=r0;r;r=r->next)
 {
   if(r->id0>n) break;
   for(a=r->atm;a;a=a->next)
   {
     lat.getcell(a,cut);
     for(i=0;i<lat.nget;i++)
     {
       d=TRES.distance(a->xyz,lat.obtain[i]->atm->xyz);
       if(d>cut) continue;
       lat.putoff(lat.obtain[i]->atm);
       lat.obtain[i]->atm->flag=1;
     }
   }
 }
}

void Pdb::addsidechain() {

 for(Chn *c=chn;c;c=c->next) c->addsidechain();

}


void Pdb::allresid0() {

 int n=0;
 for(Chn *c=chn;c;c=c->next)
 for(Res *r=c->res;r;r=r->next) {
     r->id0=n++;
 }
} 

int Pdb::manyres() {
 
 int n=0;
 for(Chn *c=chn;c;c=c->next)
 for(Res *r=c->res;r;r=r->next) {
	n++;
 }
 return n;
}
int Pdb::getresnum() {

 int n=0;
 for(Chn *c=chn;c;c=c->next) n+=c->number;
 return n;
}

void Pdb::dihedral(FILE *fp){
 for(Chn *c=chn;c;c=c->next) c->dihedral(fp);
}

void Pdb::dihedraloption(FILE *fp,int m){
 //R   10 -     176.5   -63.3   -45.2   177.7    63.3    67.5  -178.9    -2.2
 if(m==0) {
        fprintf(fp,"!            OMEGA   PHI     PSI      X1       X2      X3     X4        X5\n");
 }
 else if(m==1) {
        fprintf(fp,"!            PHI     PSI      X1       X2      X3     X4        X5\n");
 }
 else if(m==2) {
        fprintf(fp,"!            PHI     PSI\n");
 }
 else if(m==3) {
        fprintf(fp,"!             X1       X2      X3     X4        X5\n");
 }
 else if(m==4) {
        fprintf(fp,"!            OMEGA\n");
 }

 for(Chn *c=chn;c;c=c->next) c->dihedraloption(fp,m);
}

void Pdb::dihedraloption(char *s,int m){
 if(s==0) return;
 FILE *fp=fopen(s,"w");
 if(fp==0) return; 
 //R   10 -     176.5   -63.3   -45.2   177.7    63.3    67.5  -178.9    -2.2
 if(m==0) {
	fprintf(fp,"!            OMEGA   PHI     PSI      X1       X2      X3     X4        X5\n");
 }
 else if(m==1) {
	fprintf(fp,"!            PHI     PSI      X1       X2      X3     X4        X5\n");
 }
 else if(m==2) {
	fprintf(fp,"!            PHI     PSI\n");
 }
 else if(m==3) {
	fprintf(fp,"!             X1       X2      X3     X4        X5\n");
 }
 else if(m==4) {
	fprintf(fp,"!            OMEGA\n");
 }
 for(Chn *c=chn;c;c=c->next) c->dihedraloption(fp,m);
 fclose(fp);
}

void Pdb::transfer(float *s,int f){

Res *r;
int i;
i=0;
for(Chn *c=chn;c;c=c->next)
for(r=c->res;r;r=r->next)
{
 r->transfer(s+i,f);
 i+=r->tres->number*3;
}


}

void Pdb::giveresmoreid(){

 for(Chn *s=chn;s;s=s->next)
 for(Res *r=s->res;r;r=r->next)r->giveid();

}

Res * Pdb::isres(int n) {
 Chn *s;
 Res *r;
 for(s=chn;s;s=s->next)
 for(r=s->res;r;r=r->next) {
    if(r->id0==n) return r;
 }
 return 0;
}

int Pdb::manyatm() {

 Chn *s;
 Res *r;
 Atm *a;
 int n=0;
 for(s=chn;s;s=s->next)
 for(r=s->res;r;r=r->next) 
 for(a=r->atm;a;a=a->next) n++;
 return n;
}

int Pdb::maxatmid0() {

 Chn *s;
 Res *r;
 Atm *a;
 int n=0;
 for(s=chn;s;s=s->next)
 for(r=s->res;r;r=r->next)
 for(a=r->atm;a;a=a->next) n=max(n,a->id0);
 return n;

}

int Pdb::maxresid0() {

 Chn *s;
 Res *r;
 int n=0;
 for(s=chn;s;s=s->next)
 for(r=s->res;r;r=r->next)
 n=max(n,r->id0);
 return n;

}

void Pdb::transform(int f,int g) {

Chn *s;
for(s=chn;s;s=s->next) s->transform(f,g);
}

void Pdb::transform(int f) {

Chn *s;
for(s=chn;s;s=s->next) s->transform(f);
}

void Pdb::center() {

Chn *s;
Res *r;
Atm *a;
float coo[3];
int i,n;
for(i=0;i<3;i++) coo[i]=0;
n=0;
for(s=chn;s;s=s->next)
for(r=s->res;r;r=r->next)
for(a=r->atm;a;a=a->next) {
for(i=0;i<3;i++) coo[i]+=a->xyz[i];
n++;
}
for(i=0;i<3;i++) coo[i]=coo[i]/n;

for(s=chn;s;s=s->next)
for(r=s->res;r;r=r->next)
for(a=r->atm;a;a=a->next) {
for(i=0;i<3;i++) a->xyz[i]-=coo[i];
}

}

int Pdb::manychain() {

	int n=0;

	Chn *c;

	for(c=chn;c;c=c->next) n++;

	return n;
}

int Pdb::ischainexist(Chn* c){

	Chn *s;
	for(s=chn;s;s=s->next) if(s==c) return 1;

	return 0;
}

int Pdb::ishetatmchainexist(Chn* c){

	Chn *s;
	for(s=hetatm;s;s=s->next) if(s==c) return 1;

	return 0;
}

void Pdb::removechain(Chn* c) {

	if(ischainexist(c)==0) return;

	Chn *s,*s0;
	
	s0=0;

	for(s=chn;s;s=s->next){
		if(s==c) {
			if(s0==0) {	//first chain	
				chn=c->next;
				c->next=0;
				delete c;
				return;
			}
			else {
				s0->next=c->next;
				c->next=0;
				delete c;
				return;
			}
		}
		s0=s;
	}

}
void Pdb::removehetatmchain(Chn* c) {

	if(ishetatmchainexist(c)==0) return;

	Chn *s,*s0;
	
	s0=0;

	for(s=hetatm;s;s=s->next){
		if(s==c) {
			if(s0==0) {	//first chain	
				hetatm=c->next;
				c->next=0;
				delete c;
				return;
			}
			else {
				s0->next=c->next;
				c->next=0;
				delete c;
				return;
			}
		}
		s0=s;
	}

}
void Pdb::relinkchain(Chn *c0,Chn *c) {

	if(ischainexist(c)==0) return;

	Chn *s,*s0;
	
	s0=0;

	for(s=chn;s;s=s->next){
		if(s==c) {
			if(s0==0) {	//first chain	
				chn=c->next;
				c->next=0;
				addmorechain(c0,c);
				//delete c;
				return;
			}
			else {
				s0->next=c->next;
				c->next=0;
				addmorechain(c0,c);
				//delete c;
				return;
			}
		}
		s0=s;
	}
}
void Pdb::relinkhetatmchain(Chn *c0,Chn *c) {

	if(ishetatmchainexist(c)==0) return;

	Chn *s,*s0;
	
	s0=0;

	for(s=hetatm;s;s=s->next){
		if(s==c) {
			if(s0==0) {	//first chain	
				hetatm=c->next;
				c->next=0;
				addmorechain(c0,c);
				//delete c;
				return;
			}
			else {
				s0->next=c->next;
				c->next=0;
				addmorechain(c0,c);
				//delete c;
				return;
			}
		}
		s0=s;
	}
}
void Pdb::addmorechain(Chn *c0,Chn *c) {
	Chn *s;
	for(s=c0;s->more;s=s->more);
	s->more=c;
}


void Pdb::readseqcard() {
Chn *c;
for(c=chn;c;c=c->next)readseqcard(c->id);
}


void Pdb::readseqcard(char cid) {

	Strhandler cc;

	int n=cc.gettotnum(info);

	int i,j,nn,m,k;

	int len=strlen("SEQRES");

	Chn *cn;
	char c;
	char *s,**t;

	for(i=0;i<n;i++) {

		if(strncmp(info[i],"SEQRES",len)!=0) continue;
		
		if(cid!='0'&&info[i][11]!=cid) continue;
		
		c=info[i][11];

		cn=ischain(c);

		if(cn==0) continue;

		if(cn->seqcard==0) {
			nn=atoi(info[i]+12)+1000;
			cn->seqcard=new char[nn];
			for(j=0;j<nn;j++) cn->seqcard[j]='\0';
		}
		s=cc.getsegment(info[i], 19,70);
		s=cc.clearendchar(s," \n\r");
						
		t=cc.pairbytoken(s," ");
		
		m=cc.gettotnum(t);
	
		for(j=0;j<m;j++) {
			Tres *tr=TRES[t[j]];
			if(tr==0) c='X';
			else      c=tr->name;
			k=strlen(cn->seqcard);
			cn->seqcard[k]=c;
		}
		s=cc.strdel(s);
		t=cc.strdel(t);
	}
}

float Pdb::getarea() {
float a=0;
for(Chn *e=chn;e;e=e->next) a+=e->getarea();
return a;
}

void Pdb::writedelphi(char *file,Res *r,int n,int start,int end) {

	char line[1000];

	sprintf(line,"%s.pdb",file);

	writepdbused(line);

	sprintf(line,"%s.crg",file);

	FILE *fp=fopen(line,"w");
	fprintf(fp,"atom__resnumbc_charge_\n");
	for(Res *rr=r;rr;rr=rr->next) {
		if(rr->id0>=n) break;
		for(Atm *aa=rr->atm;aa;aa=aa->next) {
			if(aa->tatm->eng->charge==0) continue;
			if(aa->tatm->name[1]=='H') {
				if(aa->tatm->bond[0]->id<start||aa->tatm->bond[0]->id>=end) continue;
			}
			else {
				if(aa->tatm->id<start||aa->tatm->id>=end) continue;
			}
			fprintf(fp,"%s  %s%4i %8.3f\n",aa->name,rr->tres->name3,rr->id+rr->chn->start,aa->tatm->eng->charge);
		}
	}
	fclose(fp);
}

void Pdb::writeoutrotamer(FILE *fpp,char c) {

	
	Res *r=chn->isres(c,0);
	
	if(r==0) return;

	if(fpp==0) return;

	int n=0;
	for(Res *rr=r->more;rr;rr=rr->more) {

		rr->write(fpp,4);
		fflush(fpp);
		Res *rtt=rr;
		fprintf(fpp,"HEAD 1 %8.3f %8.3f %8.3f\n",
    		rtt->temp[0],rtt->temp[1],rtt->temp[2]);
    		fprintf(fpp,"HEAD 2 %8.3f %8.3f %8.3f\n",
    		rtt->temp[3],rtt->temp[4],rtt->temp[5]);
    		fprintf(fpp,"HEAD 3 %8.3f %8.3f %8.3f\n",
    		rtt->temp[6],rtt->temp[7],rtt->temp[8]);
    		fflush(fpp);
	}

	cerr<<"total for residue :"<<r->name<<" "<<"is "<<n<<endl;
}

Atm *Pdb::getatmbyid0(int d)
{

	for(Chn *c=chn;c;c=c->next)
	for(Res *r=c->res;r;r=r->next)
	for(Atm *a=r->atm;a;a=a->next)
	if(a->id0==d) return a;
	return 0;
}

void Pdb::configureatmid0() {

	int n=0;
	for(Chn *c=chn;c;c=c->next)
        for(Res *r=c->res;r;r=r->next)
        for(Atm *a=r->atm;a;a=a->next) {
		a->id0=n;
		n++;
	}
}

void Pdb::configuremoreatmid0() {

        int n=0;
	Res *rr,*r;
	Atm *a;
	Chn *c;
        for(c=chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {
		int nr=0;
		for(rr=r;rr;rr=rr->more) {
			nr++;
			if(nr>53) continue;
        		for(a=r->atm;a;a=a->next) {
                		a->id0=n;
                		n++;
        		}
		}
	}
}

void Pdb::setatmoid() {

        for(Chn *c=chn;c;c=c->next)
        for(Res *r=c->res;r;r=r->next)
        for(Atm *a=r->atm;a;a=a->next) {
                a->oid=a->id0;
        }
}

void Pdb::setresid0back() {

        for(Chn *c=chn;c;c=c->next)
        for(Res *r=c->res;r;r=r->next)
	{
                r->id0=r->oid;
        }
}


void Pdb::setatmid0back() {

        for(Chn *c=chn;c;c=c->next)
        for(Res *r=c->res;r;r=r->next)
        for(Atm *a=r->atm;a;a=a->next) {
                a->id0=a->oid;
        }
}


void Pdb::setlastresorder() {

	for(Chn *c=chn;c;c=c->next) c->setlastresorder();

}

void Pdb::setthreestatesec() {
	for(Chn *c=chn;c;c=c->next) c->setthreestatesec();
}

void Pdb::deletecontact() {
        for(Chn *c=chn;c;c=c->next) c->deletecontact();
}

void Pdb::setallnear() {
	Chn *c;
	Res *r;
	Atm *a;
	for(c=chn;c;c=c->next)
	for(r=c->res;r;r=r->next)
	for(a=r->atm;a;a=a->next) {
        	a->allnear(3,0);
	}
}	

void Pdb::getmaxxyz(float *coo) {

	int i;
	for(i=0;i<3;i++) coo[i]=chn->res->atm->xyz[i];

	Chn *c;
        Res *r;
        Atm *a;
        for(c=chn;c;c=c->next)
        for(r=c->res;r;r=r->next)
        for(a=r->atm;a;a=a->next) {
		for(i=0;i<3;i++) coo[i]=max(coo[i],a->xyz[i]);
        }
}

void Pdb::addxyz(float *coo) {

        int i;

        Chn *c;
        Res *r;
        Atm *a;
        for(c=chn;c;c=c->next)
        for(r=c->res;r;r=r->next)
        for(a=r->atm;a;a=a->next) {
                for(i=0;i<3;i++) a->xyz[i]=coo[i];
        }
}

void Pdb::surface(int m,float probe){

Chn *c;
for(c=chn;c;c=c->next) c->surface(m,probe);
area=getarea();
}

void Pdb::writesurface(char *ss) {
Chn *c;
Res *res;
Atm *atm;
FILE *fp=fopen(ss,"w");
fprintf(fp,"!ATM       RES             RAD      AREA\n");
for(c=chn;c;c=c->next) {
for(res=c->res;res;res=res->next)  {
for(atm=res->atm;atm;atm=atm->next){
fprintf(fp,"%s %5i %s %5i %c %8.3f %8.3f\n",atm->name,atm->oid,res->tres->name3,res->oid,c->id,atm->tatm->eng->radius,atm->area);
}
}
}
 
fprintf(fp,"\n\n");
fprintf(fp,"!RES          AREA\n");
for(c=chn;c;c=c->next) {
for(res=c->res;res;res=res->next)  {
fprintf(fp,"%s %5i %c %8.3f\n",res->tres->name3,res->oid,c->id,res->area);
}
}

fprintf(fp,"\n\n");
if(manychain()>1) {
fprintf(fp,"!Chain\n");
for(c=chn;c;c=c->next) {
fprintf(fp,"%c %8.3f\n",c->id,c->area);
}
}
fprintf(fp,"total area of pdb:%8.3f\n",area);
fclose(fp);
}

void Pdb::writesurface(FILE *fp) {
Chn *c;
Res *res;
Atm *atm;
//FILE *fp=fopen(ss,"w");
fprintf(fp,"!ATM       RES             RAD      AREA\n");
for(c=chn;c;c=c->next) {
for(res=c->res;res;res=res->next)  {
for(atm=res->atm;atm;atm=atm->next){
fprintf(fp,"%s %5i %s %5i %c %8.3f %8.3f\n",atm->name,atm->oid,res->tres->name3,res->oid,c->id,atm->tatm->eng->radius,atm->area);
}
}
}

fprintf(fp,"\n\n");
fprintf(fp,"!RES          AREA\n");
for(c=chn;c;c=c->next) {
for(res=c->res;res;res=res->next)  {
fprintf(fp,"%s %5i %c %8.3f\n",res->tres->name3,res->oid,c->id,res->area);
}
}

fprintf(fp,"\n\n");
if(manychain()>1) {
fprintf(fp,"!Chain\n");
for(c=chn;c;c=c->next) {
fprintf(fp,"%c %8.3f\n",c->id,c->area);
}
}
fprintf(fp,"total area of pdb:%8.3f\n",area);
fclose(fp);
}

char *Pdb::getname() {

if(name==0) return 0;

char *ss=strdup(name);

Strhandler cc;

int n=cc.strncmplast(ss,".pdb",4);

if(n==1) return ss;

char *s=cc.strstrlast(ss,".pdb");

if(s) *s='\0';

return ss;
}

char *Pdb::getname(char *name0) {

if(name0==0) return 0;

char *ss=strdup(name0);

Strhandler cc;

int n=cc.strncmplast(ss,".pdb",4);

if(n==1) return ss;

char *s=cc.strstrlast(ss,".pdb");

if(s) *s='\0';

return ss;
}
 
void Pdb::addhatoms(int f){
Chn *c;
for(c=chn;c;c=c->next) c->addhatoms(f);
}
void Pdb::header(){
Chn *c;
for(c=chn;c;c=c->next) c->header();
}

void Pdb::headerpdbfix(){
Chn *c;
for(c=chn;c;c=c->next) c->headerpdbfix(c->res,c->lastres()->id0+10000);
}

void Pdb::writehbondlist(char *s) {
if(s==0) return;
FILE *fp;
fp=fopen(s,"w");
if(fp==0) return;
 
Chn *c;
fprintf(fp,"!RES ID ATM RES ID  ATM    DIS\n");
for(c=chn;c;c=c->next) c->writehbondlist(fp);
for(c=chn;c;c=c->next) c->writessbondlist(fp);
fclose(fp);
}

void Pdb::writehbondlist(FILE *fp) {
if(fp==0) return;
 
Chn *c;
fprintf(fp,"!RES ID ATM RES ID  ATM    DIS\n");
for(c=chn;c;c=c->next) c->writehbondlist(fp);
for(c=chn;c;c=c->next) c->writessbondlist(fp);
 
}

void Pdb::writeseq(char *s) {
if(s==0) return;
FILE *fp=fopen(s,"w");
if(fp==0) return;
Chn *c;
for(c=chn;c;c=c->next)
{
        fprintf(fp,">P1;%s\n",name);
        int n1,n2;
        if(c->seqcard==0) {n1=0;n2=0;}
        else {n2=strlen(c->seqcard);n1=1;}
        fprintf(fp,"sequence:%s:%i:%c:%i:%c\n",name,n1,c->id,n2,c->id);
        if(c->seqcard==0) continue;
        fprintf(fp,"%s\n",c->seqcard);
}
fclose(fp);
}
void Pdb::writeseq(FILE *fp) {
if(fp==0) return;
Chn *c;
for(c=chn;c;c=c->next)
{
	fprintf(fp,">P1;%s\n",name);
	int n1,n2;
	if(c->seqcard==0) {n1=0;n2=0;}
	else {n2=strlen(c->seqcard);n1=1;}
	fprintf(fp,"sequence:%s:%i:%c:%i:%c\n",name,n1,c->id,n2,c->id);
	if(c->seqcard==0) continue;
	fprintf(fp,"%s\n",c->seqcard);
}
}
void Pdb::setallatmid0() {

	Chn *c;
	Res *r;
	Res *r1;
	Atm *a;
	
	int n=0;
	for(c=chn;c;c=c->next) 
	for(r=c->res;r;r=r->next)
	for(r1=r;r1;r1=r1->more) 
	for(a=r1->atm;a;a=a->next) {
		n++;
		a->id0=n;
	}
}
