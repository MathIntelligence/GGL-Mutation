#include<new.h>
#include"../src/sourceall.h"

Seqdbs::Seqdbs()
{

chn_num=0;
seqchn=0;
}

Seqdbs::~Seqdbs()
{
if(seqchn) delete seqchn;
} 

void Seqdbs::readdbs(char *filnam0,int flag)
{
FILE *fp;
char  *database_directory;
char filnam[256],line[256],**seqprf;
Seqchn *schn,*schn1;
int   i,k,n,h,m[3];

if(filnam0==0) {cerr<<"null file name!\n";return;}
if((fp=fopen(filnam0,"r"))==NULL){
  database_directory=resources["database_directory"];
  sprintf(filnam,"%s/%s",database_directory,filnam0);
  if((fp=fopen(filnam,"r"))==NULL) {cerr<<"no open:"<<line<<endl;return;}
}

seqprf=new char*[3];
for(k=0;k<3;k++) seqprf[k]=new char[3000];
n=0;m[0]=m[1]=m[2]=0;schn=0;schn1=0;

while(fgets(line,256,fp)!=NULL) {
    
  k=strlen(line);
 
  if(k<2||strncmp(line,"    ",4)==0) {n=0;continue;}
  if(line[0]=='!') continue; 
  if(line[0]=='>') 
  {     
      if(schn!=0&&m[0]!=0) addseq(schn,seqprf,m[0],flag);       
      schn1=schn;
      schn=new Seqchn;
      chn_num++;
      if(seqchn==0) seqchn=schn;
      else schn1->next=schn; h=strlen(line);
      schn->name=new char[h];
      strcpy(schn->name,line+1);
      schn->name[h-1]='\0';
      for(i=0;i<h;i++) if(schn->name[i]=='\n') schn->name[i]='\0';
      n=0;m[0]=m[1]=m[2]=0;continue;
  }
  if(flag==0)
  {
    n++;   
    strcpy(seqprf[n-1]+m[n-1],line);   
    m[n-1]+=k-1;    
  }
  else
  {
    strcpy(seqprf[0]+m[0],line);
    m[0]+=k-1;
  }
}

  if(schn!=0&&m[0]!=0)addseq(schn, seqprf,m[0],flag);      
  for(k=0;k<3;k++) delete [] seqprf[k];
  delete [] seqprf;   
}


void Seqdbs::addseq(Seqchn *schn, char **seqprf,int m,int flag)
{
int i,j;
   j=-1;
   for(i=0;i<m;i++)
   {
     if(seqprf[0][i]==' '||seqprf[0][i]=='\0'||seqprf[0][i]=='\n') continue;
     j++;
     seqprf[0][j]= seqprf[0][i];
      seqprf[1][j]= seqprf[1][i];
        seqprf[2][j]= seqprf[2][i];
   }
   m=j+1; 
   seqprf[0][m]='\0';
   seqprf[1][m]='\0';
   seqprf[2][m]='\0';
   
   if(flag==0)schn->seqn=new char[m *3+3];
   else schn->seqn=new char[m+1];
   strncpy(schn->seqn,seqprf[0],m);
   if(flag==0){
   strncpy(schn->seqn+m+1,seqprf[1],m);
   strncpy(schn->seqn+2*m+2,seqprf[2],m);
   schn->seqn[2*m+1]='\0'; schn->seqn[3*m+2]='\0';
   }
   schn->seqn[m]='\0';
}
 
