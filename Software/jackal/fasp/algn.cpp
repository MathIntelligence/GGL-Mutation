#include"sourceclass.h"
extern Sequence sequence;
extern ResourceList resources;
Algn::Algn()
{
  target=0;query=0; score=0;
  result=0;routine=0;flag=-1;
  prfm=0;opncost=gapcost=0;
  outfmt=3;tlen=qlen=0;slen=0;
}

void Algn::seqgive(char *t,char *q,int flg) 
{
int k1,k2,k,i;

destroy();
score=0;opncost=gapcost=0;
target=t;query=q; flag=flg; 
outfmt=3;slen=0;
tlen=strlen(target);
qlen=strlen(query);
prfm=new int[2*tlen+1]; 
for(i=0;i<2*tlen+1;i++) prfm[i]=0;
for(i=0;i<tlen&&flag!=1;i++)
{
  if(target[i+tlen+1]=='h') prfm[i]=1;
  else if(target[i+tlen+1]=='e') prfm[i]=2;
  else prfm[i]=3;
}

k1=0;
for(i=1;i<tlen&&flag!=1;i++)
{
  if(prfm[i]!=prfm[k1])
  {
    k2=i;
    for(k=k1;k<k2;k++) prfm[k+tlen+1]=k1*tlen+k2-1;
    k1=k2;
  }
}

}
 
void Algn::alignment(short *matx,int gap)
{
  float **dmn; 
  int **dmn1;
  int k,i,j,n,k1,k2,k3;
  float d,w1;
 
//  tlen=strlen(target); n=strlen(query);
  
//penalty parameter

  if(matx){opncost=3*matx[9];gapcost=0.2;}
  else {opncost=-2;gapcost=0.05;}

// initialization of alignment variable table 

  dmn=new float*[tlen];  dmn1=new int*[tlen];  
  routine=new int[tlen+qlen+1];//path of alignment 
  for(i=0;i<tlen+qlen+1;i++)routine[i]=-1;

  for(i=0;i<tlen;i++){ dmn [i]=new float[qlen]; dmn1[i]=new int[qlen];}

  for(i=0;i<tlen;i++) for(j=0;j<qlen;j++) {dmn[i][j]=0;dmn1[i][j]=-1;}
   
  scoring(matx,dmn);

//calculation of the maximum of alignment score of alignment
  
  for(i=0;i<tlen;i++)
  {  
     for(j=0;j<qlen;j++)
     {
       if(!i||!j) continue;//bypass the first cell
       d=dmn[i-1][j-1];
       dmn1[i][j]=(i-1)*qlen+j-1;

       for(k=i-2;k>max(-1,i-gap);k--)
       { 
         w1=gapping(dmn1,i,j,k,j-1);
         if(dmn[k][j-1]+w1>d){d=dmn[k][j-1]+w1;dmn1[i][j]=k*qlen+j-1;}
       }      
       for(k=j-2;k>max(-1,j-gap);k--)
       {
         w1=gapping(dmn1,i,j,i-1,k);
         if(dmn[i-1][k]+w1>d) { d=dmn[i-1][k]+w1;dmn1[i][j]=(i-1)*qlen+k;}
       }
       dmn[i][j]+=d; 
   }
 }

//find the maximum of score position
  w1=-10000;
  for(i=0;i<tlen;i++){
    if(w1<dmn[i][qlen-1]){w1=dmn[i][qlen-1];routine[0]=i*qlen+qlen-1;}
  }
  for(i=0;i<qlen;i++){
    if(w1<dmn[tlen-1][i]){w1=dmn[tlen-1][i];routine[0]=(tlen-1)*qlen+i;}
  }

// retrieve the path
// k is the lenght of path
  k=0;
  for(; ;)
  {
   k++;
   i=routine[k-1]/qlen;
   j=routine[k-1]%qlen;
   routine[k]=dmn1[i][j]; 
   if(dmn1[i][j]==-1) break;
  }

  slen=0;
  j=0;k3=0;n=0;
  for(i=k-1;i>-1;i--){
   k1=routine[i]/qlen;
   k2=routine[i]%qlen;
   slen+=1+abs(k1-k2-j);
   j=k1-k2;
   if((query[k2]==target[k1])&&target[k1]!='?')k3++;
   if(flag==0||flag==2)
   {
     if(target[k1+2*tlen+2]=='1') 
     {
       if(strchr("AVLIFWCMYPG",target[k1])) 
       {
          if(strchr("AVLIFWCMYPG",query[k2])) n++;
       }
       else if(strchr("KRHDENQTSPG",target[k1]))
       {
          if(strchr("KRHDENQTSPG",target[k1]))n++;
       }
     }
   }
  }
  i=0;
  for(j=0;j<tlen;j++) if(target[j+2*tlen+2]=='1') i++;
  n=i-n;
  if(i) n=float(n)/i*100;
  else  n=0;   
  slen+=tlen-k1+qlen-k2-2;//overtail
  if(slen)score=float(k3)/slen+n;
  else score=0;

  for(i=0;i<tlen;i++){
   delete [] dmn[i]; delete [] dmn1[i]; 
  }
  delete [] dmn; delete [] dmn1;
}


 
char **Algn::output(FILE *fp,int fff)

{

int i,j,k,k1,k2,k3,m;
int mm[5]; 

if(fff==0)outfmt=5;
else outfmt=3;

//calculation the lenght of alignment sequence

k=-1;
while(routine[++k]!=-1);
k--;

//allocating space for sequence alignment output;

result=new char*[outfmt];

for(i=0;i<outfmt;i++) 
{
  result[i]=new char[slen+1];
  for(j=0;j<slen+1;j++) result[i][j]=' ';
}


// aligned sequence
j=0;
m=0;

for(i=k;i>-1;i--)
{

k1=routine[i]/qlen;
k2=routine[i]%qlen;

if(j>k1-k2)
for(k3=0;k3<j-(k1-k2);k3++){
result[1][m]=query[k1-j+k3];
result[0][m]='-';
m++;
}
else if(j<k1-k2)
for(k3=0;k3<k1-k2-j;k3++){
result[0][m]=target[k2+j+k3];
if(outfmt==5){
result[3][m]=target[k2+j+k3+tlen+1];
result[4][m]=target[k2+j+k3+2*tlen+2];
}
result[1][m]='-';
m++;
}

result[0][m]=target[k1];
if(outfmt==5){result[3][m]=target[k1+tlen+1];result[4][m]=target[k1+2*tlen+2];}
result[1][m]=query[k2];
 
m++;
j=k1-k2;

}

if(k1<strlen(target)-1)
for(i=k1+1;i<strlen(target);i++){
result[0][m]=target[i];
if(outfmt==5) {result[3][m]=target[i+tlen+1];result[4][m]=target[i+2*tlen+2];}
result[1][m]='-';
m++;
}
else if(k2<strlen(query)-1)
for(i=k2+1;i<strlen(query);i++){
result[1][m]=query[i];
result[0][m]='-';
m++;
}

result[0][m]='\0';
result[1][m]='\0';
if(outfmt==5){ result[3][m]='\0'; result[4][m]='\0';}

for(i=0;i<m;i++)
{
  if(result[0][i]==result[1][i]) result[2][i]='*';

  else
  { 
    j=0;
    while(sequence.similar_aa[j++][0]!='\0')
    {
       if(strchr(sequence.similar_aa[j-1],result[0][i])!=NULL&&\
          strchr(sequence.similar_aa[j-1],result[1][i])!=NULL)
       { result[2][i]='.';break;}
    }
    if(result[2][i]!='\0') continue;
    result[2][i]=' ';
  }
}

result[2][m]='\0';

if(!fp) return result;
 
for(i=0;i<outfmt;i++) mm[i]=-1;
fprintf(fp,"\n");

for(k=0;k<slen/60+1;k++) {

for(i=0;i<outfmt;i++)
{ 
   for(j=0;j<60;j++){
    mm[i]++;
    if(mm[i]>=slen) break;
   if(mm[i]%10==0) fprintf(fp," ");
    fprintf(fp,"%c",result[i][mm[i]]);
   }
   fprintf(fp,"\n");
   if(i==outfmt-1) fprintf(fp,"\n");
}
}

cout<<slen*(score-int(score))<<","<<slen<<","<<score<<endl;
return result;

}


Algn::~Algn(){ destroy();}
void Algn::destroy()
{
int i;

if(routine) delete routine; routine=0;
if(prfm) delete [] prfm; prfm=0;
if(result)
{ 
  for(i=0;i<outfmt;i++) delete [] result[i]; 
  delete [] result; 
}
result=0; 
}


float Algn::gapping(int **dmn1,int i,int j,int i1,int j1)
{
float w1;
int k,p1,p2,p3,m,n;
if(flag==1)
{
  if(i-i1>1)
  { 
     p1=0;p2=0;p3=0;//p1:polar;p2:hydrophobic;p3:PG
     
     m=(5-i+i1)/2; if(m<0) m=0; 
     n=0;
     for(k=i1-m;k<=i+m;k++)
     {
        if(k<0||k>=tlen) continue;
        n++;
        if(strchr(sequence.polar_aa,target[k]))p1++;
        if(strchr(sequence.hydro_aa,target[k]))p2++;
        if(target[k]=='P'||target[k]=='G')p3++;
     }
     if(p3+p2==n) {w1=-opncost*1.5-(i-i1+j-j1-2)*gapcost; return w1;}
     if(p3+p1==n) {w1=-opncost*0.5-(i-i1+j-j1-2)*gapcost; return w1;}
     if(p3>0)     {w1=-opncost*0.7-(i-i1+j-j1-2)*gapcost; return w1;}
     w1=-opncost-(i-i1+j-j1-2)*gapcost; return w1;
     
  } 
  else
  {
     p1=0;p2=0;p3=0;//p1:polar;p2:hydrophobic;p3:PG
     m=(5-j+j1)/2; if(m<0) m=0;
     n=0;
     for(k=j1-m;k<=j+m;k++)
     {
        if(k<0||k>=qlen) continue;
        n++;
        if(strchr(sequence.polar_aa,query[k]))p1++;
        if(strchr(sequence.hydro_aa,query[k]))p2++;
        if(query[k]=='P'||query[k]=='G')p3++;
     }
     if(p3+p2==n) {w1=-opncost*1.5-(i-i1+j-j1-2)*gapcost; return w1;}
     if(p3+p1==n) {w1=-opncost*0.5-(i-i1+j-j1-2)*gapcost; return w1;}
     if(p3>0)     {w1=-opncost*0.7-(i-i1+j-j1-2)*gapcost; return w1;}
     w1=-opncost-(i-i1+j-j1-2)*gapcost; return w1;
  }
}
else if(flag==0)
{
  if(prfm[i]==prfm[i1])
  {
    if(prfm[i]!=3) w1=-opncost*1.5-(i-i1+j-j1-2)*gapcost;  
    else w1=-opncost*0.7-(i-i1+j-j1-2)*gapcost; 
  }
  else
  {
    if(prfm[i]==3||prfm[i1]==3) w1=-opncost*0.7-(i-i1+j-j1-2)*gapcost;
    else w1=-opncost*1.5-(i-i1+j-j1-2)*gapcost;
  } 
  return w1;
}
else if(flag==2)
{
 return -opncost*1.5-(i-i1+j-j1-2)*gapcost;  
}
else return -opncost*1.5-(i-i1+j-j1-2)*gapcost;
}


void Algn::scoring(short *matx, float **dmn)
{

int k,i,j,k1,k2;

  if(matx==0&&(flag==1||flag==0)) 
  {cerr<<"please specify the scoring matrice\n";exit(0);}
// initialization of the alignment score

  k=strlen(sequence.amino_acid_order);

  for(i=0;i<tlen;i++)
  {   
    if(flag==1||flag==0)
    {
       for(k1=0;k1<k;k1++) 
       if(target[i]==sequence.amino_acid_order[k1]) break;
       for(j=0;j<qlen;j++)
       {
        for(k2=0;k2<k;k2++)
        if(query[j]==sequence.amino_acid_order[k2]) break;
        if(k1==k||k2==k) {dmn[i][j]=0; continue;}//non-exist
        if(k1<k2) dmn[i][j]=matx[k2*(k2+1)/2+k1]; 
        else dmn[i][j]=matx[k1*(k1+1)/2+k2]; 
       }
     
    }
    else if(flag==2)
    {
       if(strchr(sequence.hydro_aa,target[i])) k1=-1;
       else if(strchr(sequence.polar_aa,target[i])) k1=1;
       else k1=0;
       if(k=='G'||k=='P') k1=0;

       for(j=0;j<qlen;j++)
       {  
         if(strchr(sequence.hydro_aa,query[j])) k2=-1;
         else if(strchr(sequence.polar_aa,query[j])) k2=1;
         else k2=0;
         if(k=='G'||k=='P') k2=0; //GP neutral
         dmn[i][j]=k1*k2;
         if(query[j]==target[i]&&target[i]!='?') dmn[i][j]=1.5;
       }
    }
  }
}
