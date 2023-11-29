#include"../model/source.h"
Rcs RCS(".model");
Tres TRES;
main(int argc,char *argv[])
{
Domain  dom;
FILE *fp,*fp1;
Build c;
char line[200],line0[200];
Res *rrr[1000000];
int nnn;
Fish cc(&c);
Pdb *pdb; 
Chn *chn;
TRES.flag=100;
TRES.read("tres");
Res *r,*r1,*r2;
int i,j,k,m,n,kk;
float y,z; 
Atm *a,*a1;
Tatm *tatm;
char c1;

fp=fopen("lis","r");
pdb=c.pdb;
//read whole pdb
i=0;j=0;
while(fgets(line,200,fp)!=NULL)
{
  if(strlen(line)<10) continue;
  if(i&&j){pdb->next=new Pdb(&c);pdb=pdb->next;}
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
for(pdb=c.pdb;pdb;pdb=pdb->next)
{ 
  for(chn=pdb->chn;chn;chn=chn->next)
  {
  if(chn->res==0) continue;
  i++;
  j+=chn->number;
  cc.header(chn);
  chn->dihedral(0);
  }
}

cerr<<"the number of chains: "<<i<<" "<endl;
cerr<<"the number of total residues: "<<j<<endl;

//kick out the residues not standard or could not be linked
nnn=0;
for(pdb=c.pdb;pdb;pdb=pdb->next)
{
  for(chn=pdb->chn;chn;chn=chn->next)
  {
  for(r=chn->res;r;r=r->next)
  { 
    if(r->number!=r->tres->number)  continue;
    if(r->nemp==0)  continue;
    k=0;
    for(a=r->atm;a;a=a->next)
    {
      if(a->tatm->rotate==1)
      if(a->chi>400) {k=1;break;}
    }
    if(k) continue;
    if(r)rrr[nnn++]=r; 
  }  
  }
} 

for(i=0;i<nnn;i++)
{
  rrr[i]->next=0;
  r=rrr[i];
  m=0;

  for(a=r->atm;a;a=a->next)
  {
    if(a->tatm->rotate==1)
    {
      if(a->chi<0) a->chi=360+a->chi;
      if(a->tatm->balance==1)
      {
        if(a->chi>180) a->chi=a->chi-180;
      }
      if(a->chi<0){cerr<<"strange!"<<endl;}
      j=a->chi;
      n=5;
      if(m<2) n=5;
      else  
      {
        if(strchr("FWYR",r->name))
        if(m<4)n=5;
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

//fp=fopen("/hosts/dodo/home/xiang/rotamer","w");

///fp1=fopen("/hosts/dodo/home/xiang/foot","w");

fp=fopen("rotamer","w");
fp1=fopen("foot","w");
int ii;
int rt,rrt[30],rrrt[30];
rt=0;
for(ii=0;ii<2;ii++)
{
for(i=0;i<30;i++)
{ 
  cerr<<char(i+'A')<<endl;
  rrt[i]=0;rrrt[i]=0;
  mmm=0;
  for(j=0;j<nnn;j++)   
  {
    if((rrr[j]->name-'A')==i) 
    {  
       if(ii==0)inres[mmm++]=rrr[j];
       else
       {
         if(j==0||j==nnn-1) continue; 
         if(rrr[j-1]->name==rrr[j]->name)
         inres[mmm++]=rrr[j];
       }
    }
  }
  rrrt[i]=mmm;
  for(j=0;j<mmm;j++) for(k=0;k<7;k++) angle[j][k]=0;
  for(k=0;k<7;k++) emp[k]=0;

  for(j=0;j<mmm;j++)
  {
    k=0;m1=0;
    for(a=inres[j]->atm;a;a=a->next)
    {
      if(a->tatm->rotate==1) 
      {  
        angle[j][k]=a->chi;emp[k]=a->area+0.2;
        k++;
        if(k>=6) break;
      }
    }
    if(m1&&m1!=k) cerr<<"m1 not equal k"<<endl;
    m1=k;
  }

  for(j=0;j<mmm;j++)
  {
    if(inres[j]==0) continue;
    n2=0;
    for(n1=0;n1<mmm;n1++)
    {
      if(n1==j) continue;
      if(inres[n1]==0) continue;
      if(abs(angle[j][0]-angle[n1][0])<emp[0]) have[n2++]=n1;
    } 
    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][1]-angle[have[n1]][1])<emp[1]) have[n2++]=have[n1];
    if(k==2) {for(n1=0;n1<n2;n1++)inres[have[n1]]=0; continue;} 
    

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][2]-angle[have[n1]][2])<emp[2]) have[n2++]=have[n1];
    if(k==3) {for(n1=0;n1<n2;n1++)inres[have[n1]]=0;continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][3]-angle[have[n1]][3])<emp[3]) have[n2++]=have[n1];
    if(k==4) {for(n1=0;n1<n2;n1++)inres[have[n1]]=0;continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][4]-angle[have[n1]][4])<emp[4]) have[n2++]=have[n1];
    if(k==5) {for(n1=0;n1<n2;n1++)inres[have[n1]]=0;continue;}

    n3=n2;n2=0;
    for(n1=0;n1<n3;n1++)
    if(abs(angle[j][5]-angle[have[n1]][5])<emp[5]) have[n2++]=have[n1];
    if(k==6) {for(n1=0;n1<n2;n1++)inres[have[n1]]=0;continue;}
  }
  for(j=0;j<mmm;j++)
  {
    if(inres[j]==0) continue;

    if(ii==0)inres[j]->dihedral(fp,-2);
    else inres[j]->dihedral(fp1,-2);
    fflush(fp); 
    rt++;
    rrt[i]++;
  }

 } 

for(i=0;i<30;i++)
{
if(rrt[i])
{
if(ii==0)fprintf(fp,"%c %d %d\n",char(i+'A'),rrrt[i],rrt[i]);
else fprintf(fp1,"%c %d %d\n",char(i+'A'),rrrt[i],rrt[i]);
}
}
}
fprintf(fp,"the total residues: %d\n",rt);
exit(0);

}
