#include"sourceclass.h"
#include"matrices.h"
//extern Sequence sequence;
//extern ResourceList resources;

Prfalgn::Prfalgn(Protein *ptr)
{
protein=ptr; 
target=protein->chain->get_chain_seq();
query=0;
routine=0;
result=0;
matrix=0;
score=0;    
important=0;
}

Prfalgn::Prfalgn(char *s)
{

protein=0;
target=s;
query=0;
routine=0;
result=0;
matrix=0;
score=0;

}

void Prfalgn::matrice()
{
int i,j,k,m,n,l;
float x;

matrix=new float*[Amino_Acid_Num];
for(i=0;i<Amino_Acid_Num;i++) {
matrix[i]=new float[Amino_Acid_Num];
for(j=0;j<Amino_Acid_Num;j++)matrix[i][j]=0;
}

l=strlen(sequence.amino_acid_order);
k=0;

for(i=0;i<l;i++)
{
 m=sequence.amino_acid_order[i]-'A';
 for(j=0;j<=i;j++)
 {
   n=sequence.amino_acid_order[j]-'A';
   x=matrices[k++];
   matrix[m][n]=x;
   matrix[n][m]=x;
 }
}   
x=-Opncost-Gapcost;
for(i=0;i<27;i++) {matrix[26][i]=x;matrix[i][26]=x;}
matrix[27]=0;

}




void Prfalgn::prfmatrice (char *filnam)
{

int i,j,k;
float maxtr[Amino_Acid_Num][Amino_Acid_Num]; 
char *working_directory;
char line[256];
FILE *fp;
char a,b;

matrice();


for(i=0;i<27;i++)
for(j=0;j<27;j++)
maxtr[i][j]=matrix[i][j];

(*this).~Prfalgn();

important=new int[strlen(target)+1];
for(i=0;i<strlen(target)+1;i++) important[i]=0;

matrix=new float*[strlen(target)+1];
for(i=0;i<strlen(target)+1;i++) matrix[i]=new float[Amino_Acid_Num];

working_directory=resources["working_directory"];
sprintf(line,"%s/%s",working_directory,filnam);
if((fp=fopen(line,"r"))==NULL){cerr<<"no open:"<<line<<endl;}

i=0;
while(fgets(line,256,fp))
{
if(line[0]=='!') continue; 
 j=0;

 k=strlen(line);

 for(; ;){
    i=atoi(line+j)-protein->chain->begin; 
    important[i]=1;
//    flag[i++]=atoi(line+j);
    while(line[j++]!=','&&j<=k-1);
    if(j==k) break;
 }

}

for(i=0;i<strlen(target);i++)
{
  a=target[i];
  if(a>'Z'||a<'A') a=26;
  else a=a-'A';
  
  for(j=0;j<strlen(sequence.amino_acid_order);j++)
  {
    b=sequence.amino_acid_order[j];      
    b=b-'A'; 
    if(important[i]) matrix[i][b]=maxtr[a][b];
    else {
          if(maxtr[a][b]>0)matrix[i][b]=maxtr[a][b];
          else matrix[i][b]=0;
    }

  }

}

  matrix[strlen(target)]=0;

}

void Prfalgn::align(char *seq) 
{

float **dmn,**pmn,**qmn;
char a,b;
int m,n,i,j,flag; 
int k,k1,k2 ;
float d,w1;

query=seq;

m=strlen(target);

n=strlen(query);

// initialization of alignment variable table 
dmn=new float*[m];//the three matrices for calculation
pmn=new float*[m];
qmn=new float*[m];
routine=new int[m+n+1];//path of alignment 
for(i=0;i<m+n+1;i++)routine[i]=-1;

for(i=0;i<m;i++){
  dmn[i]=new float[n];
  pmn[i]=new float[n];
  qmn[i]=new float[n];
}

for(i=0;i<m;i++)
for(j=0;j<n;j++){
pmn[i][j]=-10000;
qmn[i][j]=-10000;
}


i=-1;
while(matrix[++i]);
if(i==27)flag=0;
else flag=1;
  

// initialization of the alignment score

for(i=0;i<m;i++)
{
  if(flag) k=i;
  else
  { 
    a=target[i];
    if(a>'Z'||a<'A') k=26;
    else k=a-'A';
  }
   for(j=0;j<n;j++)
  {
    b=query[j];
    if(b>'Z'||b<'A') k1=26;
    else k1=b-'A';
    dmn[i][j]=matrix[k][k1];
  }
}




//calculation of the maximum of alignment score of alignment
w1=-Opncost-Gapcost;
for(i=0;i<m;i++)
{
  for(j=0;j<n;j++)
  {
    if(!i||!j) continue;
    d=dmn[i][j];
    dmn[i][j]=dmn[i-1][j-1];
    if(i>1)pmn[i][j]=max((dmn[i-2][j-1]+w1),(pmn[i-1][j]-Gapcost));
    if(j>1)qmn[i][j]=max((dmn[i-1][j-2]+w1),(qmn[i][j-1]-Gapcost));
    dmn[i][j]=max(dmn[i][j],pmn[i][j]);
    dmn[i][j]=max(dmn[i][j],qmn[i][j])+d;
  }
}

//find the maximum of score position
w1=-10000;
for(i=0;i<m;i++){
if(w1<dmn[i][n-1]){w1=dmn[i][n-1];routine[0]=i*n+n-1;}
}
for(i=0;i<n;i++){
if(w1<dmn[m-1][i]){w1=dmn[m-1][i];routine[0]=(m-1)*n+i;}
}


// retrieve the path
// k is the lenght of path

k=0;

re200:

k1=routine[k]/n;
if(flag) i=k1;
else{
  a=target[k1];
  if(a>'Z'||a<'A') i=26;
  else i=a-'A';
}

k2=routine[k]%n;
b=query[k2];
if(b>'Z'||b<'A') j=26;
else j=b-'A';
dmn[k1][k2]-=matrix[i][j];


if(abs(dmn[k1-1][k2-1]-dmn[k1][k2])<0.001)
{
routine[++k]=(k1-1)*n+k2-1;
goto re100;
}


for(i=k1-2;i>-1;i--)
{
w1=-Opncost-Gapcost*(k1-1-i);
if(abs(dmn[i][k2-1]+w1-dmn[k1][k2])<0.001){
routine[++k]=i*n+k2-1;
goto re100;
}
}
for(i=k2-2;i>-1;i--)
{
w1=-Opncost-Gapcost*(k2-1-i);
if(abs(dmn[k1-1][i]+w1-dmn[k1][k2])<0.001){
routine[++k]=(k1-1)*n+i;
goto re100;
}

}

re100:

if(routine[k]/n>0&&routine[k]%n>0) goto re200;

for(i=0;i<m;i++){
delete [] dmn[i];
delete [] pmn[i];
delete [] qmn[i];
}
delete [] dmn;
delete [] pmn;
delete [] qmn;

}


 
char **Prfalgn::output(FILE *fp)

{
int i,j,k,k1,k2,k3,m,n,m2;

//calculation the lenght of alignment sequence

k=-1;
while(routine[++k]!=-1);
k--;
k3=0;

n=strlen(query);
m=strlen(target);
j=0;
for(i=k;i>-1;i--){
k1=routine[i]/n;
k2=routine[i]%n;
k3+=1+abs(k1-k2-j);
j=k1-k2;
}
k3+=m-k1+n-k2-2;//overtail



//allocating space for sequence alignment output;


result=new char*[4];

for(i=0;i<4;i++) 
{
  result[i]=new char[k3+1];
  for(j=0;j<k3+1;j++) result[i][j]='\0';
}


// aligned sequence
j=0;
m=0;

for(i=k;i>-1;i--)
{

k1=routine[i]/n;
k2=routine[i]%n;

if(j>k1-k2)
for(k3=0;k3<j-(k1-k2);k3++){
result[1][m]=query[k1-j+k3];
result[0][m]='-';
result[3][m]=' ';
m++;
}

else if(j<k1-k2)
for(k3=0;k3<k1-k2-j;k3++){
result[0][m]=target[k2+j+k3];
if(important[k2+j+k3]) result[3][m]='?';
else result[3][m]=' ';
result[1][m]='-';
m++;
}

result[0][m]=target[k1];
if(important[k1]) result[3][m]='?';
else result[3][m]=' ';
result[1][m]=query[k2];
m++;
j=k1-k2;

}

if(k1<strlen(target)-1)
for(i=k1+1;i<strlen(target);i++){
result[0][m]=target[i];
if(important[i]) result[3][m]='?';
else result[3][m]=' ';
result[1][m]='-';
m++;
}
else if(k2<strlen(query)-1)
for(i=k2+1;i<strlen(query);i++){
result[1][m]=query[i];
result[0][m]='-';
result[3][m]=' ';
m++;
}


//delete [] routine;
result[0][m]='\0';
result[1][m]='\0';
result[3][m]='\0';

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


//s=get_chain_seq();
k3=strlen(result[0]); //sequence length
 

//for(k=0;k<residue_num/Seq_Output_Len_ALL+1;k++) {

i=-1;m=-1;n=-1;m2=-1;
for(k=0;k<k3/Seq_Output_Len_ALL+1;k++) {
  
   for(j=0;j<Seq_Output_Len_ALL;j++){
    i++;
    if(i>=k3) break;
    if(i%Seq_Output_Len_ind==0) fprintf(fp,"%4d ",i+1);
    fprintf(fp,"%c",result[0][i]);
   }

   fprintf(fp,"\n");

   for(j=0;j<Seq_Output_Len_ALL;j++){
    m++;
    if(m>=k3) break;
    if(m%Seq_Output_Len_ind==0)  fprintf(fp,"     " );
    fprintf(fp,"%c",result[1][m]);
   }
   fprintf(fp,"\n");
  
   for(j=0;j<Seq_Output_Len_ALL;j++){
    n++;
    if(n>=k3) break;
    if(n%Seq_Output_Len_ind==0) fprintf(fp,"     " );
    fprintf(fp,"%c",result[2][n]);
   }

   fprintf(fp,"\n");
   

   for(j=0;j<Seq_Output_Len_ALL;j++){
    m2++;
    if(m2>=k3) break;
    if(m2%Seq_Output_Len_ind==0) fprintf(fp,"     " );
    fprintf(fp,"%c",result[3][m2]);
   }

   fprintf(fp,"\n\n\n"); 
   
}
 
return result;

}


Prfalgn::~Prfalgn()
{
int i,j;

if(routine)delete routine;

if(result){
  for(i=0;i<4;i++) delete [] result[i];
  delete [] result;
}
j=0;
while(matrix[j++]);
for(i=0;i<j;i++) delete matrix[i];
delete [] matrix;
if(important) delete [] important;
}
