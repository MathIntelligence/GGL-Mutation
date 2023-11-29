#include"../src/sourceall.h"
Seqchn::Seqchn()
{
name=0;
seqn=0;
next=0;
flag=1;
}

Seqchn::~Seqchn()
{
if(seqn) delete [] seqn;
if(name) delete [] name;
}


void Seqchn::write(FILE *fp,int flg)
{
  
int n,m,k,j,i;
char *s;
int number;

i=-1;
m=-1;
n=-1;
flag=flg;
s=seqn;
number=strlen(s);
fprintf(fp,">%s",name);
fprintf(fp,"\n\n");

for(k=0;k<number/Seq_Output_Len_ALL+1;k++) {
   for(j=0;j<Seq_Output_Len_ALL;j++){
    i++;
    if(i>=number) break;
    if(i%Seq_Output_Len_ind==0) fprintf(fp," ");
    fprintf(fp,"%c",s[i]);
   }
   fprintf(fp,"\n");
   if(flag==0)
   {
     for(j=0;j<Seq_Output_Len_ALL;j++){
       m++;
       if(m>=number) break;
       if(m%Seq_Output_Len_ind==0) fprintf(fp," ");
       fprintf(fp,"%c",s[1+number+m]);
     }
     fprintf(fp,"\n");
     for(j=0;j<Seq_Output_Len_ALL;j++){
       n++;
       if(n>=number) break;
       if(n%Seq_Output_Len_ind==0) fprintf(fp," ");
       fprintf(fp,"%c",s[2+2*number+n]);
     }
     fprintf(fp,"\n\n");if(i>=number-1)  return; 
  }

}

}

