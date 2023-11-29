#include"source.h"
Rotamer::Rotamer()
{
}
Rotamer::~Rotamer()
{
}

void Rotamer::calcrotamer(int argc, char **argv) {
FILE *fp;
char line[200],line0[200];
Res *rrr[1000000];
int nnn;
Pdb *pdb,*pdb0; 
Chn *chn;
TRES.flag=100;
TRES.read("tres");
Res *r,*r1,*r2;
int i,j,k,m,n,kk;
float y,z; 
Atm *a,*a1;
Tatm *tatm;
char c1;
int step;

step=atoi(argv[2]);
fp=fopen(argv[1],"r");
pdb0=new Pdb;pdb=pdb0;
//read whole pdb
i=0;j=0;
while(fgets(line,200,fp)!=NULL)
{
  if(strlen(line)<10) continue;
  if(i&&j){pdb->next=new Pdb;pdb=pdb->next;}
  line[4]='\0';
  c1=line[5];
  if(c1<'A'||c1>'Z') c1='1';
  sprintf(line0,"%s_mod",line);
  cerr<<line0<<" "<<i<<endl;
  j=pdb->read(line0,c1);
  i++;
}

//header it! 
i=0;j=0;
for(pdb=pdb0;pdb;pdb=pdb->next)
{ 
  for(chn=pdb->chn;chn;chn=chn->next)
  {
  if(chn->res==0) continue;
  i++;
  j+=chn->number;
  chn->header();
  chn->dihedral((FILE*)0);
  }
}

cerr<<"the number of chains: "<<i<<" "<<endl;
cerr<<"the number of total residues: "<<j<<endl;

//kick out the residues not standard or could not be linked
nnn=0;
for(pdb=pdb0;pdb;pdb=pdb->next)
{
  Res *r1=0;
  for(chn=pdb->chn;chn;chn=chn->next)
  {
  for(r=chn->res;r;r=r->next)
  { 
    if(r1==0) {r1=r;continue;}
    if(r->number!=r->tres->number)  {r1=r;continue;}
    
    if(r1->atm->next->xyz[0]!=r->temp[0]) {
	        r1=r;
		continue;
    }
    int ne=0;
    for(a=r->atm;a;a=a->next) ne++;
    if(ne!=r->tres->number)  {r1=r;continue;}
    if(r->nemp==0) { r1=r;continue;}
    k=0;
    for(a=r->atm;a;a=a->next)
    {
      if(a->tatm->rotate==1)
      if(a->chi>400) {k=1;break;}
    }
    if(k) {r1=r;continue;}
    if(r)rrr[nnn++]=r; 
    r1=r;
  }  
  }
} 

if(selfenergy==1) {
float *temp=new float[nnn+100];
int   *temp1=new int[nnn+100];
for(i=0;i<nnn;i++) {
	r=rrr[i];
	float e=r->clashnoring(5,0);
	temp[i]=e;
	r->rotenergy=e;
	
}
Qsort cc;

cc.sort(temp,nnn,temp1);
Res *rrr1[1000000]; 
for(i=0;i<nnn;i++) {
	int j=temp1[i];
	rrr1[i]=rrr[j];
}
for(i=0;i<nnn;i++) {
	rrr[i]=rrr1[i];
}
}
else if(selfenergy==2) {
float *temp=new float[nnn+100];
int   *temp1=new int[nnn+100];
for(i=0;i<nnn;i++) {
	r=rrr[i];
	float e=r->clashnoring(5,0);
	temp[i]=e;
	r->rotenergy=e;
	
}
}

for(i=0;i<nnn;i++)
{
  rrr[i]->next=0;
  r=rrr[i];
  m=0;
  r->rotdegree=step;
  for(a=r->atm;a;a=a->next)
  {
    if(a->tatm->rotate==1)
    {
      a->good=a->chi;
      if(a->good<0) a->good=360+a->good;
      if(a->tatm->balance==1)
      {
        if(a->good>180) a->good=a->good-180;
      }
      if(a->good<0){cerr<<"strange!"<<endl;}
      j=a->good;
      n=step;
      /*
      if(r->tres->name=='R'||r->tres->name=='K') {
	if(n<15&&n>5) n=15;
      }
      */
      if(m<2) n=500;
      
      else if(m>=4&&r->tres->name=='R') {
	if(n<15&&n>5) n=20;
      }
      else if(m>=4&&r->tres->name=='K') {
        if(n<15&&n>5) n=20;
      }
      
      else if(m>=4&&(r->tres->name=='Q'||r->tres->name=='E'||r->tres->name=='M')) {
        if(n<15&&n>5) n=20;
      }
      
      else if(a->tatm->name[1]=='O'&&(r->tres->name=='Y'||r->tres->name=='T'||r->tres->name=='S')) {
	if(n<15&&n>5) n=20;
      }
      /*
      else if((r->tres->name=='R'||r->tres->name=='K')) {
        if(fabs(step-10)<0.1) n=step+10;
        //else if(fabs(step-20)<0.1) n=step+10;
      }
      */
      else  
      {
        if(strchr("FWY",r->name))
        if(m<4)n=step;
      }
      a->area=n;
      m++;     
    }
  }
}

//strcpy(c.pdb->name,"resmer");
cerr<<"the number of total residues after kick is: "<<nnn<<endl;

float tem[10];
int temp[10],x,temp1[10];
Res **inres;
int mmm;
int **angle,emp[7];
int *have;
int n1,n2,n3;
int m1;
inres=new Res*[nnn];

angle=new int*[nnn];
for(j=0;j<nnn;j++) angle[j]=new int[7];
 
 
have=new int[nnn];


FILE *fpp[11]; 

for(i=0;i<11;i++)
{
sprintf(line,"%s_side_rotamer%d_%d",argv[1],100-i*2,step);
fpp[i]=fopen(line,"w");
}



int rt[11],rrt[11][30],rrrt[30];
for(j=0;j<11;j++)rt[j]=0;
for(i=0;i<30;i++)
{
  cerr<<char(i+'A')<<endl;
  for(j=0;j<11;j++)rrt[j][i]=0;
  rrrt[i]=0;
  mmm=0;
  for(j=0;j<nnn;j++)   
  {
    if((rrr[j]->name-'A')==i) inres[mmm++]=rrr[j];
  }
  rrrt[i]=mmm;
  Qsort tt;
  Res *rtt[100000],*rtt1[100000];
  int ntt,jn,jj,temp1[100000];
  float temp[100000];
  float fe;
  //
  for(j=0;j<mmm;j++) {
	inres[j]->flag=1;
	Atm *a;
	for(a=inres[j]->atm;a;a=a->next) {
		a->occup=a->area;
		a->area=0;
	}
  }
  
  if(inres[0]->name!='A') {
  	calcrotamer0(inres,mmm,10);
  	for(j=0;j<mmm;j++) {
		temp[j]=-inres[j]->rotchance;
		inres[j]->rotchance=0;
  	} 
  	tt.sort(temp,mmm,temp1);
  	for(j=0;j<mmm;j++) {
		int jj=temp1[j];
		rtt[j]=inres[jj];
  	}
  	for(j=0;j<mmm;j++) {
		inres[j]=rtt[j];
  	}
  }
  for(j=0;j<mmm;j++) {
	inres[j]->rotchance=0;
	inres[j]->flag=0;
	Atm *a;
	for(a=inres[j]->atm;a;a=a->next) {
		a->area=a->occup;
		a->occup=0;
	}
  }
  //
  for(j=0;j<mmm;j++) for(k=0;k<7;k++) angle[j][k]=0;
  for(k=0;k<7;k++) emp[k]=0;

  for(j=0;j<mmm;j++)
  {
    k=0;m1=0;
    for(a=inres[j]->atm;a;a=a->next)
    {
      if(a->tatm->rotate==1) 
      {  
        angle[j][k]=a->good;emp[k]=a->area+0.2;
        k++;
        //if(k>=6) break;
      }
    }
    if(m1&&m1!=k) cerr<<"m1 not equal k"<<endl;
    m1=k;
  }

  for(j=0;j<mmm;j++)
  {
    if(inres[j]==0) continue;
    inres[j]->flag+=1;
    n2=0;
    for(n1=0;n1<mmm;n1++)
    {
      if(n1==j) continue;
      if(inres[n1]==0) continue;
      if(abs(angle[j][0]-angle[n1][0])<emp[0]) {have[n2++]=n1;}//inres[j]->flag+=1;}
      else if(360-abs(angle[j][0]-angle[n1][0])<emp[0]) {have[n2++]=n1;}
    } 
    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][1]-angle[have[n1]][1])<emp[1]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][1]-angle[have[n1]][1])<emp[1]) {have[n2++]=have[n1];}
    if(k==2) {for(n1=0;n1<n2;n1++){inres[have[n1]]=0;inres[j]->flag+=1;} continue;} 
    

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][2]-angle[have[n1]][2])<emp[2]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][2]-angle[have[n1]][2])<emp[2]) {have[n2++]=have[n1];}
    if(k==3) {for(n1=0;n1<n2;n1++){inres[have[n1]]=0;inres[j]->flag+=1;}continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][3]-angle[have[n1]][3])<emp[3]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][3]-angle[have[n1]][3])<emp[3]) {have[n2++]=have[n1];}
    if(k==4) {for(n1=0;n1<n2;n1++){inres[have[n1]]=0;inres[j]->flag+=1;}continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][4]-angle[have[n1]][4])<emp[4]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][4]-angle[have[n1]][4])<emp[4]) {have[n2++]=have[n1];}
    if(k==5) {for(n1=0;n1<n2;n1++){inres[have[n1]]=0;inres[j]->flag+=1;}continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][5]-angle[have[n1]][5])<emp[5]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][5]-angle[have[n1]][5])<emp[5]) {have[n2++]=have[n1];}
    if(k==6) {for(n1=0;n1<n2;n1++){inres[have[n1]]=0;inres[j]->flag+=1;}continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][6]-angle[have[n1]][6])<emp[6]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][6]-angle[have[n1]][6])<emp[6]) {have[n2++]=have[n1];}
    if(k==7) {for(n1=0;n1<n2;n1++){inres[have[n1]]=0;inres[j]->flag+=1;}continue;}

  }
  
 
  ntt=0;
  for(j=0;j<mmm;j++)
  {
    if(inres[j]!=0)
    {
      rtt1[ntt++]=inres[j];
    }
  }
  calcrotamer0(rtt1,ntt,20);  
  

  fe=0; int fe0=0; 
  for(j=0;j<ntt;j++)
  {
    temp[j]=-rtt1[j]->rotchance; 
    fe+=rtt1[j]->rotchance;
    fe0+=rtt1[j]->flag;
  }
  //cerr<<"....."<<fe<<" "<<fe0<<endl; 
  
  for(j=0;j<ntt;j++) {
	if(fe) rtt1[j]->rotchance=rtt1[j]->rotchance*1.0/fe0;
  }
   
  for(j=0;j<ntt;j++) {
        temp[j]=-8.31*0.3*log(rtt1[j]->rotchance);
  }
  tt.sort(temp,ntt,temp1);

  jj=0;
  Atm *aaa;

  for(j=0;j<ntt;j++)
  {
    jn=temp1[j];
    rtt[jj++]=rtt1[jn];
    //for(aaa=rtt1[jn]->atm;aaa;aaa=aaa->next)aaa->good=1;
  }
  rtt[jj]=0;
  if(jj!=ntt) {
	cerr<<" error:"<<jj<<" "<<ntt<<endl;
  }
  jj=ntt;
  for(j=0;j<ntt;j++) rtt1[j]=rtt[j];
  
  calcrotamer1(rtt,jj,20);
  for(j=0;j<ntt;j++) rtt[j]=rtt1[j];
  //....
  calcrotamer0(rtt1,ntt,40);  
  

  fe=0; fe0=0; 
  for(j=0;j<ntt;j++)
  {
    temp[j]=-rtt1[j]->rotchance; 
    fe+=rtt1[j]->rotchance;
    fe0+=rtt1[j]->flag;
  }
  //cerr<<"....."<<fe<<" "<<fe0<<endl; 
  
  for(j=0;j<ntt;j++) {
	if(fe) rtt1[j]->rotchance=rtt1[j]->rotchance*1.0/fe0;
  }
   
  for(j=0;j<ntt;j++) {
        temp[j]=-8.31*0.3*log(rtt1[j]->rotchance);
  }
  tt.sort(temp,ntt,temp1);

  jj=0;
 // Atm *aaa;

  for(j=0;j<ntt;j++)
  {
    jn=temp1[j];
    rtt[jj++]=rtt1[jn];
    //for(aaa=rtt1[jn]->atm;aaa;aaa=aaa->next)aaa->good=1;
  }
  rtt[jj]=0;
  if(jj!=ntt) {
	cerr<<" error:"<<jj<<" "<<ntt<<endl;
  }
  jj=ntt;
  for(j=0;j<ntt;j++) rtt1[j]=rtt[j];
 
  //....
  calcrotamer1(rtt,jj,40);
   for(j=0;j<ntt;j++) rtt[j]=rtt1[j];
 
 //....
  calcrotamer0(rtt1,ntt,60);  
  

  fe=0; fe0=0; 
  for(j=0;j<ntt;j++)
  {
    temp[j]=-rtt1[j]->rotchance; 
    fe+=rtt1[j]->rotchance;
    fe0+=rtt1[j]->flag;
  }
  //cerr<<"....."<<fe<<" "<<fe0<<endl; 
  
  for(j=0;j<ntt;j++) {
	if(fe) rtt1[j]->rotchance=rtt1[j]->rotchance*1.0/fe0;
  }
   
  for(j=0;j<ntt;j++) {
        temp[j]=-8.31*0.3*log(rtt1[j]->rotchance);
  }
  tt.sort(temp,ntt,temp1);

  jj=0;
 // Atm *aaa;

  for(j=0;j<ntt;j++)
  {
    jn=temp1[j];
    rtt[jj++]=rtt1[jn];
    //for(aaa=rtt1[jn]->atm;aaa;aaa=aaa->next)aaa->good=1;
  }
  rtt[jj]=0;
  if(jj!=ntt) {
	cerr<<" error:"<<jj<<" "<<ntt<<endl;
  }
  jj=ntt;
  for(j=0;j<ntt;j++) rtt1[j]=rtt[j];
 
  //....

 calcrotamer1(rtt,jj,60);
  for(j=0;j<ntt;j++) rtt[j]=rtt1[j];

 //....
  calcrotamer0(rtt1,ntt,90);  
  

  fe=0; fe0=0; 
  for(j=0;j<ntt;j++)
  {
    temp[j]=-rtt1[j]->rotchance; 
    fe+=rtt1[j]->rotchance;
    fe0+=rtt1[j]->flag;
  }
  //cerr<<"....."<<fe<<" "<<fe0<<endl; 
  
  for(j=0;j<ntt;j++) {
	if(fe) rtt1[j]->rotchance=rtt1[j]->rotchance*1.0/fe0;
  }
   
  for(j=0;j<ntt;j++) {
        temp[j]=-8.31*0.3*log(rtt1[j]->rotchance);
  }
  tt.sort(temp,ntt,temp1);

  jj=0;
 // Atm *aaa;

  for(j=0;j<ntt;j++)
  {
    jn=temp1[j];
    rtt[jj++]=rtt1[jn];
    //for(aaa=rtt1[jn]->atm;aaa;aaa=aaa->next)aaa->good=1;
  }
  rtt[jj]=0;
  if(jj!=ntt) {
	cerr<<" error:"<<jj<<" "<<ntt<<endl;
  }
  jj=ntt;
  for(j=0;j<ntt;j++) rtt1[j]=rtt[j];
 
  //....

 calcrotamer1(rtt,jj,90);   

   for(j=0;j<ntt;j++) rtt[j]=rtt1[j];

 //....
  calcrotamer0(rtt1,ntt,120);  
  

  fe=0; fe0=0; 
  for(j=0;j<ntt;j++)
  {
    temp[j]=-rtt1[j]->rotchance; 
    fe+=rtt1[j]->rotchance;
    fe0+=rtt1[j]->flag;
  }
  //cerr<<"....."<<fe<<" "<<fe0<<endl; 
  
  for(j=0;j<ntt;j++) {
	if(fe) rtt1[j]->rotchance=rtt1[j]->rotchance*1.0/fe0;
  }
   
  for(j=0;j<ntt;j++) {
        temp[j]=-8.31*0.3*log(rtt1[j]->rotchance);
  }
  tt.sort(temp,ntt,temp1);

  jj=0;
 // Atm *aaa;

  for(j=0;j<ntt;j++)
  {
    jn=temp1[j];
    rtt[jj++]=rtt1[jn];
    //for(aaa=rtt1[jn]->atm;aaa;aaa=aaa->next)aaa->good=1;
  }
  rtt[jj]=0;
  if(jj!=ntt) {
	cerr<<" error:"<<jj<<" "<<ntt<<endl;
  }
  jj=ntt;
  for(j=0;j<ntt;j++) rtt1[j]=rtt[j];
 
  //....

 calcrotamer1(rtt,jj,120);

   for(j=0;j<ntt;j++) rtt[j]=rtt1[j];

 //....
  calcrotamer0(rtt1,ntt,150);  
  

  fe=0; fe0=0; 
  for(j=0;j<ntt;j++)
  {
    temp[j]=-rtt1[j]->rotchance; 
    fe+=rtt1[j]->rotchance;
    fe0+=rtt1[j]->flag;
  }
  //cerr<<"....."<<fe<<" "<<fe0<<endl; 
  
  for(j=0;j<ntt;j++) {
	if(fe) rtt1[j]->rotchance=rtt1[j]->rotchance*1.0/fe0;
  }
   
  for(j=0;j<ntt;j++) {
        temp[j]=-8.31*0.3*log(rtt1[j]->rotchance);
  }
  tt.sort(temp,ntt,temp1);

  jj=0;
 // Atm *aaa;

  for(j=0;j<ntt;j++)
  {
    jn=temp1[j];
    rtt[jj++]=rtt1[jn];
    //for(aaa=rtt1[jn]->atm;aaa;aaa=aaa->next)aaa->good=1;
  }
  rtt[jj]=0;
  if(jj!=ntt) {
	cerr<<" error:"<<jj<<" "<<ntt<<endl;
  }
  jj=ntt;
  for(j=0;j<ntt;j++) rtt1[j]=rtt[j];
 
  //....

 calcrotamer1(rtt,jj,150);
     for(j=0;j<ntt;j++) rtt[j]=rtt1[j];

 //....
  calcrotamer0(rtt1,ntt,170);  
  

  fe=0; fe0=0; 
  for(j=0;j<ntt;j++)
  {
    temp[j]=-rtt1[j]->rotchance; 
    fe+=rtt1[j]->rotchance;
    fe0+=rtt1[j]->flag;
  }
  //cerr<<"....."<<fe<<" "<<fe0<<endl; 
  
  for(j=0;j<ntt;j++) {
	if(fe) rtt1[j]->rotchance=rtt1[j]->rotchance*1.0/fe0;
  }
   
  for(j=0;j<ntt;j++) {
        temp[j]=-8.31*0.3*log(rtt1[j]->rotchance);
  }
  tt.sort(temp,ntt,temp1);

  jj=0;
 // Atm *aaa;

  for(j=0;j<ntt;j++)
  {
    jn=temp1[j];
    rtt[jj++]=rtt1[jn];
    //for(aaa=rtt1[jn]->atm;aaa;aaa=aaa->next)aaa->good=1;
  }
  rtt[jj]=0;
  if(jj!=ntt) {
	cerr<<" error:"<<jj<<" "<<ntt<<endl;
  }
  jj=ntt;
  for(j=0;j<ntt;j++) rtt1[j]=rtt[j];
 
  //....

 calcrotamer1(rtt,jj,170);

 //....
  calcrotamer0(rtt1,ntt,40);  
  

  fe=0; fe0=0; 
  for(j=0;j<ntt;j++)
  {
    temp[j]=-rtt1[j]->rotchance; 
    fe+=rtt1[j]->rotchance;
    fe0+=rtt1[j]->flag;
  }
  //cerr<<"....."<<fe<<" "<<fe0<<endl; 
  
  for(j=0;j<ntt;j++) {
	if(fe) rtt1[j]->rotchance=rtt1[j]->rotchance*1.0/fe0;
  }
   
  for(j=0;j<ntt;j++) {
        temp[j]=-8.31*0.3*log(rtt1[j]->rotchance);
  }
  tt.sort(temp,ntt,temp1);

  jj=0;
 // Atm *aaa;

  for(j=0;j<ntt;j++)
  {
    jn=temp1[j];
    rtt[jj++]=rtt1[jn];
    //for(aaa=rtt1[jn]->atm;aaa;aaa=aaa->next)aaa->good=1;
  }
  rtt[jj]=0;
  if(jj!=ntt) {
	cerr<<" error:"<<jj<<" "<<ntt<<endl;
  }
  jj=ntt;
  for(j=0;j<ntt;j++) rtt1[j]=rtt[j];
 
  //....

  for(j=0;j<ntt;j++)
  {
    temp[j]=-(rtt1[j]->rotdegree+rtt1[j]->rotchance);      
  }
 
   
  tt.sort(temp,ntt,temp1);

  jj=0;
 

  for(j=0;j<ntt;j++)
  {
    jn=temp1[j];
    rtt[jj++]=rtt1[jn];
    //for(aaa=rtt1[jn]->atm;aaa;aaa=aaa->next)aaa->good=1;
  }
  rtt[jj]=0;


  int iu,ix,iy,iz,ig;
  float uu;
  float e0;
  Atm *ao[10],*bo[10];
  for(j=0;j<ntt;j++)
  {
   rtt[j]->nhydr=rtt[j]->flag;
  }  

  for(j=0;j<ntt;j++)
  {
    iy=rtt[j]->rotable(ao,4,100);
    for(iu=0;iu<ntt;iu++)
    {
      if(iu==j) continue;
      iz=rtt[iu]->rotable(bo,4,100);
      if(iy!=iz)  continue;
      ig=0;
      for(ix=0;ix<iy;ix++)
      {
        if(abs(ao[ix]->good-bo[ix]->good)>40) {ig=1;break;}
      }
      if(ig==0) rtt[j]->nhydr+=rtt[iu]->flag;
    }
  }
  
  for(j=0;j<ntt;j++)
  {
    for(aaa=rtt[j]->atm;aaa;aaa=aaa->next)aaa->good=1;
  }


  float rew;
  for(iu=0;iu<11;iu++)
  {
    uu=1-iu*0.02;jn=0;
/*
    for(j=0;j<ntt;j++)
    {
       if(rtt[j]==0) continue;
       jn+=rtt[j]->flag;
       if(jn/fe>uu&&j>0) continue;
       rtt[j]->dihedral(fpp[iu],-4);
       rew=1.*rtt[j]->nhydr/mmm;
       fprintf(fpp[iu]," %8.6f\n",rew);
       fflush(fpp[iu]);
       rt[iu]++;
       rrt[iu][i]++;
    }
*/
    for(j=0;j<ntt;j++)
   {
    if(rtt[j]==0) continue;
    jn+=rtt[j]->flag;
    rtt[j]->clashnoring(5,0);
    if(jn*1.0/fe0>uu&&j>0) continue;
    rtt[j]->write(fpp[iu],4);
    fflush(fpp[iu]);
    fprintf(fpp[iu],"HEAD 1 %8.3f %8.3f %8.3f\n",
    rtt[j]->temp[0],rtt[j]->temp[1],rtt[j]->temp[2]);
    fprintf(fpp[iu],"HEAD 2 %8.3f %8.3f %8.3f\n",
    rtt[j]->temp[3],rtt[j]->temp[4],rtt[j]->temp[5]);
    fprintf(fpp[iu],"HEAD 3 %8.3f %8.3f %8.3f\n",
    rtt[j]->temp[6],rtt[j]->temp[7],rtt[j]->temp[8]);
    fprintf(fpp[iu],"STAT %c %8d %8.3f %8.3f\n",
    rtt[j]->name,rtt[j]->rotdegree,rtt[j]->rotenergy,rtt[j]->rotchance);
    //cerr<<rtt[j]->name<<" "<<rtt[j]->rotdegree<<" "<<rtt[j]->rotenergy<<" "<<rtt[j]->rotchance<<endl;
    fflush(fpp[iu]);
    rt[iu]++;
    rrt[iu][i]++;
   }

  }
}
int iu;
for(iu=0;iu<11;iu++)
{
for(i=0;i<30;i++)
{
if(rrt[iu][i])
fprintf(fpp[iu],"%c %d %d\n",char(i+'A'),rrrt[i],rrt[iu][i]);
}
fprintf(fpp[iu],"the total residues: %d\n",rt[iu]);
}
}

void Rotamer::calcrotamer0(Res **rrr, int nnn,int step) {
//Domain  dom;
FILE *fp;
//Build c;
char line[200],line0[200];
//Res *rrr[1000000];
//int nnn;
//Fish cc(&c);
//Pdb *pdb,*pdb0; 
//Chn *chn;
//TRES.flag=100;
//TRES.read("tres");
Res *r,*r1,*r2;
int i,j,k,m,n,kk;
float y,z; 
Atm *a,*a1;
Tatm *tatm;
char c1;
//int step;

//step=atoi(argv[2]);
//fp=fopen(argv[1],"r");
//pdb0=new Pdb;pdb=pdb0;
//read whole pdb
 
for(i=0;i<nnn;i++)
{
  rrr[i]->next=0;
  r=rrr[i];
  m=0;

  for(a=r->atm;a;a=a->next)
  {
    if(a->tatm->rotate==1)
    {
      a->good=a->chi;
      if(a->good<0) a->good=360+a->good;
      if(a->tatm->balance==1)
      {
        if(a->good>180) a->good=a->good-180;
      }
      if(a->good<0){cerr<<"strange!"<<endl;}
      j=a->good;
      n=step;
     
      if(m<2) n=500;
       else if(m>=4&&r->tres->name=='R') {
        if(n<15&&n>5) n=20;
      }
      else if(m>=4&&r->tres->name=='K') {
        if(n<15&&n>5) n=20;
      }

      else if(m>=4&&(r->tres->name=='Q'||r->tres->name=='E'||r->tres->name=='M')) {
        if(n<15&&n>5) n=20;
      }

      else if(a->tatm->name[1]=='O'&&(r->tres->name=='Y'||r->tres->name=='T'||r->tres->name=='S')) {
        if(n<15&&n>5) n=20;
      }
      else
      {
        if(strchr("FWY",r->name))
        if(m<4)n=step;
      }
      a->area=n;
      m++;     
    }
  }
}

//strcpy(c.pdb->name,"resmer");
cerr<<"the number of total residues after kick is: "<<nnn<<endl;

float tem[10];
int temp[10],x,temp1[10];
Res **inres;
int mmm;
int **angle,emp[7];
int *have;
int n1,n2,n3;
int m1;
inres=new Res*[nnn];

angle=new int*[nnn];
for(j=0;j<nnn;j++) angle[j]=new int[7];
 
 
have=new int[nnn];
 
int rt[11],rrt[11][30],rrrt[30];
for(j=0;j<11;j++)rt[j]=0;
for(i=0;i<30;i++)
{
  
  for(j=0;j<11;j++)rrt[j][i]=0;
  rrrt[i]=0;
  mmm=0;
  for(j=0;j<nnn;j++)   
  {
    if((rrr[j]->name-'A')==i) inres[mmm++]=rrr[j];
  }
  if(mmm==0) continue;
  rrrt[i]=mmm;
  for(j=0;j<mmm;j++) for(k=0;k<7;k++) angle[j][k]=0;
  for(k=0;k<7;k++) emp[k]=0;
  cerr<<"0..."<<char(i+'A')<<endl;
  for(j=0;j<mmm;j++)
  {
    k=0;m1=0;
    for(a=inres[j]->atm;a;a=a->next)
    {
      if(a->tatm->rotate==1) 
      {  
        angle[j][k]=a->good;emp[k]=a->area+0.2;
        k++;
        //if(k>=6) break;
      }
    }
    if(m1&&m1!=k) cerr<<"m1 not equal k"<<endl;
    m1=k;
  }
   
  for(j=0;j<mmm;j++)
  {
    if(inres[j]==0) continue;
    //inres[j]->flag+=1;
    
    inres[j]->rotchance=inres[j]->flag;
    n2=0;
    for(n1=0;n1<mmm;n1++)
    {
      if(n1==j) continue;
      if(inres[n1]==0) continue;
      if(abs(angle[j][0]-angle[n1][0])<emp[0]) {have[n2++]=n1;}//inres[j]->flag+=1;}
      else if(360-abs(angle[j][0]-angle[n1][0])<emp[0]) {have[n2++]=n1;}
    } 
    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][1]-angle[have[n1]][1])<emp[1]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][1]-angle[have[n1]][1])<emp[1]) {have[n2++]=have[n1];}
    if(k==2) {
	for(n1=0;n1<n2;n1++){
		//inres[have[n1]]=0;
		//inres[j]->flag+=1;
		inres[j]->rotchance+=inres[have[n1]]->flag;
	} 
	continue;
    } 
    

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][2]-angle[have[n1]][2])<emp[2]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][2]-angle[have[n1]][2])<emp[2]) {have[n2++]=have[n1];}
    if(k==3) {
	for(n1=0;n1<n2;n1++){
		//inres[have[n1]]=0;
		//inres[j]->flag+=1;
		inres[j]->rotchance+=inres[have[n1]]->flag;
	}
	continue;
    }

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][3]-angle[have[n1]][3])<emp[3]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][3]-angle[have[n1]][3])<emp[3]) {have[n2++]=have[n1];}
    if(k==4) {for(n1=0;n1<n2;n1++){inres[j]->rotchance+=inres[have[n1]]->flag; }continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][4]-angle[have[n1]][4])<emp[4]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][4]-angle[have[n1]][4])<emp[4]) {have[n2++]=have[n1];}
    if(k==5) {for(n1=0;n1<n2;n1++){inres[j]->rotchance+=inres[have[n1]]->flag; }continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][5]-angle[have[n1]][5])<emp[5]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][5]-angle[have[n1]][5])<emp[5]) {have[n2++]=have[n1];}
    if(k==6) {for(n1=0;n1<n2;n1++){inres[j]->rotchance+=inres[have[n1]]->flag; }continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][6]-angle[have[n1]][6])<emp[6]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][6]-angle[have[n1]][6])<emp[6]) {have[n2++]=have[n1];}
    if(k==7) {for(n1=0;n1<n2;n1++){inres[j]->rotchance+=inres[have[n1]]->flag;}continue;}
  }

  
 
}

}
void Rotamer::help() {

	cerr<<"crtrotm file.list angle"<<endl;
}
void Rotamer::calcrotamer1(Res **rrr, int nnn,int step) {
//Domain  dom;
FILE *fp;
//Build c;
char line[200],line0[200];
//Res *rrr[1000000];
//int nnn;
//Fish cc(&c);
//Pdb *pdb,*pdb0; 
//Chn *chn;
//TRES.flag=100;
//TRES.read("tres");
Res *r,*r1,*r2;
int i,j,k,m,n,kk;
float y,z; 
Atm *a,*a1;
Tatm *tatm;
char c1;
//int step;

//step=atoi(argv[2]);
//fp=fopen(argv[1],"r");
//pdb0=new Pdb;pdb=pdb0;
//read whole pdb
 
for(i=0;i<nnn;i++)
{
  rrr[i]->next=0;
  r=rrr[i];
  m=0;

  for(a=r->atm;a;a=a->next)
  {
    if(a->tatm->rotate==1)
    {
      a->good=a->chi;
      if(a->good<0) a->good=360+a->good;
      if(a->tatm->balance==1)
      {
        if(a->good>180) a->good=a->good-180;
      }
      if(a->good<0){cerr<<"strange!"<<endl;}
      j=a->good;
      n=step;
      
      if(m<2) n=500;
      
      a->area=n;
      m++;     
    }
  }
}

//strcpy(c.pdb->name,"resmer");
cerr<<"the number of total residues after kick is: "<<nnn<<endl;

float tem[10];
int temp[10],x,temp1[10];
Res **inres;
int mmm;
int **angle,emp[7];
int *have;
int n1,n2,n3;
int m1;
inres=new Res*[nnn];

angle=new int*[nnn];
for(j=0;j<nnn;j++) angle[j]=new int[7];
 
 
have=new int[nnn];

 
int rt[11],rrt[11][30],rrrt[30];
for(j=0;j<11;j++)rt[j]=0;
for(i=0;i<30;i++)
{
  
  for(j=0;j<11;j++)rrt[j][i]=0;
  rrrt[i]=0;
  mmm=0;
  for(j=0;j<nnn;j++)   
  {
    if((rrr[j]->name-'A')==i) inres[mmm++]=rrr[j];
  }
  if(mmm==0) continue;
  rrrt[i]=mmm;
  for(j=0;j<mmm;j++) for(k=0;k<7;k++) angle[j][k]=0;
  for(k=0;k<7;k++) emp[k]=0;
  cerr<<"1...."<<char(i+'A')<<endl;
  for(j=0;j<mmm;j++)
  {
    k=0;m1=0;
    for(a=inres[j]->atm;a;a=a->next)
    {
      if(a->tatm->rotate==1) 
      {  
        angle[j][k]=a->good;emp[k]=a->area+0.2;
        k++;
        //if(k>=6) break;
      }
    }
    if(m1&&m1!=k) cerr<<"m1 not equal k"<<endl;
    m1=k;
  }
  
  for(j=0;j<mmm;j++)
  {
    if(inres[j]==0) continue;
    inres[j]->rotdegree=step;
    n2=0;
    for(n1=0;n1<mmm;n1++)
    {
      if(n1==j) continue;
      if(inres[n1]==0) continue;
      if(abs(angle[j][0]-angle[n1][0])<emp[0]) {have[n2++]=n1;}//inres[j]->flag+=1;}
      else if(360-abs(angle[j][0]-angle[n1][0])<emp[0]) {have[n2++]=n1;}
    } 
    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][1]-angle[have[n1]][1])<emp[1]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][1]-angle[have[n1]][1])<emp[1]) {have[n2++]=have[n1];}
    if(k==2) {
	for(n1=0;n1<n2;n1++){
		inres[have[n1]]=0;
		//inres[j]->flag+=1;
		//tempa[j]+=inres[have[n1]]->flag;
	} 
	continue;
    } 
    

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][2]-angle[have[n1]][2])<emp[2]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][2]-angle[have[n1]][2])<emp[2]) {have[n2++]=have[n1];}
    if(k==3) {
	for(n1=0;n1<n2;n1++){
		inres[have[n1]]=0;
		//inres[j]->flag+=1;
		//tempa[j]+=inres[have[n1]]->flag;
	}
	continue;
    }

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][3]-angle[have[n1]][3])<emp[3]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][3]-angle[have[n1]][3])<emp[3]) {have[n2++]=have[n1];}
    if(k==4) {for(n1=0;n1<n2;n1++){inres[have[n1]]=0;/*inres[j]->flag+=1;*/}continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][4]-angle[have[n1]][4])<emp[4]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][4]-angle[have[n1]][4])<emp[4]) {have[n2++]=have[n1];}
    if(k==5) {for(n1=0;n1<n2;n1++){inres[have[n1]]=0;/*inres[j]->flag+=1;*/}continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][5]-angle[have[n1]][5])<emp[5]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][5]-angle[have[n1]][5])<emp[5]) {have[n2++]=have[n1];}
    if(k==6) {for(n1=0;n1<n2;n1++){inres[have[n1]]=0;/*inres[j]->flag+=1;*/}continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][6]-angle[have[n1]][6])<emp[6]) {have[n2++]=have[n1];}//inres[j]->flag+=1;}
    else if(360-abs(angle[j][6]-angle[have[n1]][6])<emp[6]) {have[n2++]=have[n1];}
    if(k==7) {for(n1=0;n1<n2;n1++){inres[have[n1]]=0;/*inres[j]->flag+=1;*/}continue;}
  }
  
}

}
 
