#include "source.h"
#include <time.h>

Strhandler::Strhandler()
{
}

Strhandler::~Strhandler()
{
}

void Strhandler::init()
{
}
 
char *Strhandler::lastindexof(char *str,char *tok) {

  if(tok==0) return 0;

  int n=strlen(tok);

  char *s=str;
  char *t=0;

  while(s) {

     t=strstr(s,tok);
     if(t) {
        t+=n;
  	s=t;
	continue;
     }
     else {
	return s;
     }
  }
  return 0;
}

void Strhandler::clearendemptyspace(char *line) {

	if(line==0) return;

	char *s=strdup(line);

	s=clearendchar(s," ");

	if(s) strcpy(line,s);
	if(s) delete [] s;
}

void Strhandler::ungets(char *s, FILE *fp)
{
  char *t;

  t=s;
  while(t) 
  {
    ungetc(*t,fp);
    t++;
  }
}
int Strhandler::getasciinum(char *s)
{
  char *p;
  int n;
  p=s;n=0;
  while(s&&*p)
  {
    n+=*p;
    p++;
  }
  return n;
}

int Strhandler::notuntilstr(char *file,char *token)
{
  char *f;
  int m,tot;

  if(file==0||token==0) return 0;
 
  m=strlen(token);
  f=file;tot=0;
  while(strncmp(f,token,m)==0)
  {
    f+=m;tot+=m; 
  }
  return tot;
}

int Strhandler::notuntilchar(char *file,char *token)
{
  char *f;
  int tot;

  if(file==0||token==0) return 0;
  
  f=file;tot=0;
  while(strchr(token,*f))
  {
    f++;tot++;
  }
  return tot;
}

int Strhandler::notuntilstrlast(char *file,char *token)
{
  char *f;
  int m,tot;

  if(file==0) return 0;
  f=file;tot=strlen(file);
  m=strlen(token); 
  if(token==0) return tot;
  while(strncmplast(f,token,m)==0)
  {
    f-=m;tot-=m;
  }
  return tot;
}

int Strhandler::notuntilcharlast(char *file,char *token)
{
  char *f;
  int tot;

  if(file==0) return 0;
 
  f=file;tot=strlen(f);
  if(token==0) return tot;

  while(strchr(token,f[tot-1]))
  {
    tot--;
  }
  return tot;
}


char **Strhandler::myallocat(char **name,int maxm)
{
char  **lnks;
int   i,num;

num=gettotnum(name);
if(num>=maxm) return  name;
lnks=opnarray(maxm);
for(i=0;i<num;i++)
{
  lnks[i]=name[i];name[i]=0;
}
name=strdel(name);
name=lnks;lnks=0;
return name;
}     


char *Strhandler::myallocat(char *name,int maxm)
{
char  *lnks;
int   num;

if(name) num=strlen(name);
else num=0;

if(num>=maxm) return name;

lnks=opnstr(maxm);

if(name) strcpy(lnks,name);
name=strdel(name);

name=lnks;lnks=0;
return name;
}


int Strhandler::gettotnum(char **s)
{
   int n;
   n=0; while(s&&s[n])n++;
   return n;
}

int Strhandler::gettotnum(int **s)
{
   int n;
   n=0; while(s&&s[n])n++;
   return n;
}

int Strhandler::gettotnum(float **s)
{
   int n;
   n=0; while(s&&s[n])n++;
   return n;
}

char *Strhandler::straddup(char *a,char *b,char *c)
{
  a=straddup(a,b);
  a=straddup(a,c);
  return a;
}
char *Strhandler::straddup(char *a,char *add)
{
  char *b;

  if(add==0) return a;
  if(a==0) return opnstr(add);

  b=opnstr(strlen(a)+strlen(add)+10);

  sprintf(b,"%s%s",a,add);
  a=strdel(a); 

  a=b;

  return b;
}

char *Strhandler::assemblefilepath(char **tt)
{
   int i,n,m;
   char *s;

   n=gettotnum(tt);
   m=0; for(i=0;i<n;i++) m+=strlen(tt[i])+10;
   s=opnstr(m+10);
   for(i=0;i<n;i++)
   {
      strcat(s,tt[i]);
      if(i!=n-1) strcat(s,"/");
   }
   return s;
}

char *Strhandler::assemblefilepath(char **tt,char *v)
{
   int i,n;
   char *s;

   n=gettotnum(tt);

   s=opnstr(1000);
   strcat(s,v);
   
   for(i=0;i<n;i++)
   {
      strcat(s,tt[i]);
      if(i!=n-1) strcat(s,"/");
   }
   return s;
}


char **Strhandler::addchararray(char **newf,char **oldf)
{
   int n,i;

   char **tt;

   n=0;
   tt=opnarray(gettotnum(newf)+gettotnum(oldf)+100);

   i=0;
   while(newf&&newf[i])
   {
      tt[n++]=opnstr(newf[i++]);
   }  
   i=0;
   while(oldf&&oldf[i])
   {
      tt[n++]=opnstr(oldf[i++]);
   }  
   return tt;
}

char **Strhandler::attacharray(char **newf,char **oldf)
{
   int n,i,m;


   n=gettotnum(oldf);
   m=gettotnum(newf);
   
   for(i=0;i<n;i++)
   {
      newf[m+i]=opnstr(oldf[i]);
   }

   return newf;
}


char **Strhandler::addchartoarray(char **newf,char *oldf)
{
   int n,i;
   char line[100000];
   char **tt;

   if(oldf==0) oldf=opnstr("");
   if(newf==0) return 0;
   tt=opnarray(gettotnum(newf)+100);

   n=0; i=0;
   while(newf&&newf[i])
   {
      sprintf(line,"%s%s",oldf,newf[i]);
      tt[n++]=opnstr(line);i++;
   }  
   newf=strdel(newf); newf=tt;
   return newf;
}


char **Strhandler::delchartoarray(char **newf,char *oldf)
{
   int n,i;
   char **tt;
   char *s;

   if(oldf==0) oldf=opnstr("");
   if(gettotnum(newf)==0) return 0;

   tt=opnarray(gettotnum(newf)+100);

   n=0; i=0;
   while(newf&&newf[i])
   {
      s=strstr(newf[i],oldf);
      if(s) s+=strlen(oldf);
      else  s=newf[i];
      tt[n++]=opnstr(s);i++;
   } 
   newf=strdel(newf);newf=tt;
   return tt;
}


char *Strhandler::opnfile(char *f)
{
   char *s;
   FILE *fp; 
   if(f==0) return 0;
   fp=fopen(f,"r");
   if(fp==0) return 0;
   s=opnfile(fp);
   fclose(fp);
   return s;
}

int Strhandler::writeout(char **ss,FILE *fp)
{
  int n;

  if(ss==0) return 0;

  n=0;
  while(ss[n])
  {
    fprintf(fp,"%s<br>\n",ss[n]);
    n++;
  }
  return n;
}

char **Strhandler::shrink(char **s,int n)
{
   int i,m;
   if(s==0) return 0;
   
   m=0;
   for(i=0;i<n;i++)
   {
      if(s[i]==0) continue;
      s[m++]=s[i];
   }
   s[m]=0;
   return s;
}

char **Strhandler::opnfilebyline(char *f)
{
   char **s;
   FILE *fp;
   if(f==0) return 0;
   fp=fopen(f,"r");
   if(fp==0) return 0;
   s=opnfilebyline(fp);
   fclose(fp);
   return s;
}

char **Strhandler::opnfilebyline(char *f,int a,int b)
{
   char line[10000],**t;
   int m,n,i,len;
   FILE *fp;

   if(a>b) return 0;
   if(f==0) return 0;
   fp=fopen(f,"r");
   if(fp==0) return 0;

   t=new char*[10000]; len=10000;
   for(i=0;i<len;i++) t[i]=0;

   n=-1;m=0;
   while(fgets(line,10000,fp)!=NULL)
   {
    n++;
    if(n<a) continue;
    else if(n>b) break;
     
    clearendchar(line,"\n\t\r ");
    if(line[0]=='\0') continue;
    if(m>len-10)
    {
      t=myallocat(t,len+10000);
      //t=(char **)realloc(t,len+10000);
      len=len+10000;
    }
    t[m++]=strdup(line);
    if(m>b-a+50) break;
   }
   fclose(fp);
   return t;
}


char **Strhandler::opnfilebylinesimple(char *f)
{
   char line[10000],**t;
   int i,m,len;
   FILE *fp;

   if(f==0) return 0;
   fp=fopen(f,"r");
   if(fp==0) {
	cerr<<"could not open file ..."<<endl;
	exit(0);
   }

   t=new char*[10000]; len=10000;
   for(i=0;i<len;i++) t[i]=0;

   m=0;
   while(fgets(line,10000,fp)!=NULL)
   {
     if(m>len-10)
     {
       t=myallocat(t,len+10000);
       len=len+10000;
     }
     t[m++]=strdup(line);
   }
   t[m]=0;
   fclose(fp);
   return t;
}

char **Strhandler::opnfilebyline(FILE *fp)
{
  char line[10000],*s[50000],**t;
  int n,i;

  if(fp==0) return 0;

  for(i=0;i<50000;i++) s[i]=0;
  n=0;
  while(fgets(line,10000,fp)!=NULL)
  {
    clearendchar(line,"\n\t\r ");
    if(line[0]=='\0') continue;
    s[n++]=strdup(line); 
    if(n==50000)
    {
       cerr<<"too many file lines.ignore more than 50000\n";
       break;
    }
  }
  t=new char*[n+1];
  for(i=0;i<n;i++) t[i]=s[i];
  t[n]=0;
  return t;
}

char *Strhandler::opnfile(FILE *fp)
{
  char line[10000],*s;
  int n,m,i; 

  if(fp==0) return 0;

  m=20000;
  s=new char[20000];
  s[0]='\0';n=0;
  while(fgets(line,10000,fp)!=NULL)
  {
    i=strlen(line);
    if(i+n>m)
    {
       m+=10000;
       s=(char *)realloc(s,m);
    } 
    strcat(s,line);
    //strcat(s,"\n");
    n+=i;
  }
  s=(char *)realloc(s,n+10);
  return s;
} 

char *Strhandler::clearchar(char *tmp,char *s)
{
  int i,n,j;
  if(tmp==0||s==0) return tmp;
  n=strlen(tmp);
  if(n==0) return 0;
  if(strlen(s)==0) return tmp;
  j=0;
  for(i=0;i<n;i++)
  {
    if(strchr(s,tmp[i])) continue;
    tmp[j++]=tmp[i];
  }
  tmp[j]='\0';
  return tmp;
}

char *Strhandler::clearfirstchar(char *tmp,char *t)
{
  int i,n,j,k;
  if(tmp==0||t==0) return tmp;
  n=strlen(tmp);
  if(n==0) return 0;
  if(strlen(t)==0) return tmp;
  j=0;k=0;
  for(i=0;i<n;i++)
  {
    if(k==1) {tmp[j++]=tmp[i]; continue;}
    if(strchr(t,tmp[i])) continue; 
    k=1; 
    tmp[j++]=tmp[i];    
  }
  tmp[j]='\0';
  return tmp;
}

char *Strhandler::clearlastchar(char *tmp,char *t)
{
  int n,i;
  if(tmp==0||t==0) return tmp;
  n=strlen(tmp);
  if(n==0) return 0;
  if(strlen(t)==0) return tmp;
  for(i=n-1;i>=0;i--) 
  {
    if(strchr(t,tmp[i])==0) break;
    tmp[i]='\0';
  }
  return tmp; 
}

char *Strhandler::clearendchar(char *tmp,char *s)
{
  tmp=clearlastchar(tmp,s);
  tmp=clearfirstchar(tmp,s); 
  return tmp;
}

char *Strhandler::lowercase(char *tmp)
{
  int i,n;
  if(tmp==0) return 0;
  n=strlen(tmp);
  if(n==0) return 0;
  for(i=0;i<n;i++) 
  {
     if(tmp[i]>='A'&&tmp[i]<='Z') tmp[i]+=32;
  }
  return tmp;
}

char *Strhandler::uppercase(char *tmp)
{
  int i,n;
  if(tmp==0) return 0;
  n=strlen(tmp);
  if(n==0) return 0;
  for(i=0;i<n;i++) 
  {
     if(tmp[i]>='a'&&tmp[i]<='z') tmp[i]-=32;
  }
  return tmp;
}

char *Strhandler::opnstr(char *tmp)
{
  if(tmp==0) return 0;
  if(strlen(tmp)==0) return 0;
  return strdup(tmp);
}

int *Strhandler::opnint(int *tmp,int n)
{
  int i,*t;
  if(tmp==0) return 0;
  if(n==0)   return 0;
  t=new int[n]; 
  for(i=0;i<n;i++) t[i]=tmp[i];
  return t;
}

int *Strhandler::opnint(int n)
{
  int i,*t;
  if(n==0)   return 0;
  t=new int[n];
  for(i=0;i<n;i++) t[i]=0;
  return t;
}

char *Strhandler::opnstr(int n)
{
  char *t;
  int i;
  if(n==0)   return 0;
  t=new char[n];
  for(i=0;i<n;i++) t[i]='\0';
  return t;
}

char **Strhandler::opnstr(char **tmp)
{
  char **s;
  int i,n;
  
  if(tmp==0) return 0;
  n=0; while(tmp[n])n++;
  if(n==0) return 0;
  s=new char*[n+10];
  for(i=0;i<n;i++)
  {
    s[i]=opnstr(tmp[i]);
  }
  s[n]=0;
  return s;
}


char *Strhandler::changestr(char *tmp,char *a, char *b)
{
  return changestr(tmp,a,b,1000000);
}

char *Strhandler::changestr(char *tmp,char *a, char *b,int total)
{
   return changestrloop(tmp,a,b,total,-1);
}

char *Strhandler::changestrloop(char *tmp,char *a, char *b,int total,int loopf)
{
  int i,j,n,tot,k;
  char *s,*f;
  char *tmp0;
  char **temp;

  re300:
  
  if(tmp==0||a==0||b==0) return  tmp;

  if(strstr(tmp,a)==0) return tmp;

  i=strlen(a);j=strlen(b);n=strlen(tmp);

  if(i==0||j==0||n==0) return tmp;

  if(tmp[0]=='\0'||a[0]=='\0'||b[0]=='\0') return tmp;

  temp=opnarray(n+100);

  tot=0; f=tmp;
  while(1)
  {
    s=strstr(f,a);
    if(s)
    {
       temp[tot++]=s;
       f=s+i;
       if(tot>=total) break; 
    }
    else break;
  }
   

  tmp0=new char[n+tot*j+200];
  tmp0[0]='\0';

  f=tmp;
  for(k=0;k<tot;k++)
  {
    s=temp[k];
    *s='\0';
    strcat(tmp0,f);
    strcat(tmp0,b);
    f=s+i;
  }
  strcat(tmp0,f);
  if(tmp) tmp=strdel(tmp);
  if(temp) delete [] temp;
  tmp=tmp0;

  if(loopf!=-1&&tot>0) goto re300;
  
  re200:
  return tmp0;
}




char *Strhandler::strdel(char *s)
{
  if(s) {delete [] s;s=0;}
  return 0;
}

float *Strhandler::floatdel(float *s) {
  if(s) {delete [] s; s=0;}
  return 0;
}
int *Strhandler::intdel(int *s)
{

  if(s) {delete [] s; s=0;}
  return 0;
}

char **Strhandler::strdel(char **garage)
{
  int i,n;

  if(garage)
  {
    n=gettotnum(garage);
    for(i=0;i<n;i++)
    {
       delete [] garage[i]; 
       garage[i]=0;
    }
    delete [] garage;
    garage=0;
  }
  return 0;
}

int **Strhandler::intdel(int **garage)
{
  int i,n;

  if(garage)
  {
    n=gettotnum(garage);
    for(i=0;i<n;i++)
    {
       delete [] garage[i];
       garage[i]=0;
    }
    delete [] garage;
    garage=0;
  }
  return 0;
}

float **Strhandler::floatdel(float **garage)
{
  int i,n;

  if(garage)
  {
    n=gettotnum(garage);
    for(i=0;i<n;i++)
    {
       delete [] garage[i];
       garage[i]=0;
    }
    delete [] garage;
    garage=0;
  }
  return 0;
}


char *Strhandler::strcasestr(char *s,char *t)
{
  int i,j,ii,i0;
  int n,m;
  char a,b;
 
  if(s==0||t==0) return 0;
  n=strlen(s);m=strlen(t);
  if(n==0||m==0) return 0;
  for(i=0;i<n;i++)
  {
    ii=0;
    for(j=0,i0=i;j<m&&i0<n;j++,i0++)
    {
      if(t[j]>='A'&&t[j]<='Z') b=t[j]+32;
      else b=t[j];
      if(s[i0]>='A'&&s[i0]<='Z') a=s[i0]+32;
      else a=s[i0];
      if(a!=b) {ii=1;break;}
    }
    if(ii==0&&j==m) 
       return s+i;
  }
  return 0;
}

int *Strhandler::pairitem(char *str,char *tg)
{
  int n,i,k,j,m;
  int *pairlist;
  int len,len0;

  //tag length
  len=strlen(tg);
  if(len==0) return 0;
  //pair maximum
  len0=strlen(str);
  if(len0==0) return 0;
  pairlist=new int[2*len0];
  for(i=0;i<2*len0;i++) pairlist[i]=-1;

  n=0;m=0;
  for(i=0;i<len0;i++)
  {
     if(strncmp(str+i,tg,len)==0)
     {
        if(i)
        {
           pairlist[n*2]=m;
           pairlist[n*2+1]=i-1;
           n++;
        }
        m=i+len; i=m-1;
     }
  }
  pairlist[n*2]=m;pairlist[n*2+1]=len0-1; n++;

  m=0;
  for(i=0;i<n;i++)
  {
    k=0;
    for(j=pairlist[2*i];j<=pairlist[i*2+1];j++)
    {
       if(str[j]==' '||str[j]=='\n'||str[j]=='\r'||str[j]=='\t') continue;
       k=1;break;
    }
    if(k==0) continue;
    pairlist[2*m]=pairlist[2*i];
    pairlist[2*m+1]=pairlist[2*i+1];
    m++;
  }
  for(i=m;i<n;i++) {pairlist[2*i]=-1;pairlist[2*i+1]=-1;}
  return pairlist;
}



char *Strhandler::addflag(char *flag,char *s)
{
  int i,j;
  char *tmp; 
  if(s==0) return flag;
  if(flag==0) 
  {
     return opnstr(s);
  }
  i=strlen(flag);j=strlen(s);
  if(i==0) return opnstr(s);
  if(j==0) return flag;
  i=i+j+10;
  tmp=new char[i];
  strcpy(tmp,flag);
  strcat(tmp,s);
  flag=strdel(flag); 
  flag=tmp;
  return flag;
}

char *Strhandler::cutflag(char *flag,char *s)
{
  int n,i,j;
  char *t,*f;
  if(s==0) return flag;
  if(flag==0) return 0;
  i=strlen(flag);j=strlen(s);
  if(i==0) return 0;
  if(j==0) return flag;
  t=flag;
  n=strlen(s);
  while((f=strstr(t,s)))
  {
    for(i=0;i<n;i++) *(f+i)=' ';
  }
  flag=clearendchar(flag," ");
  return flag;
}

int Strhandler::islowerchar(char *s)
{
  int i,n;
  if(s==0) return -1;
  n=strlen(s);
  if(n==0) return -1;
  for(i=0;i<n;i++)
  {
     if(s[i]=='\0') break;
     if(s[i]=='\r'||s[i]=='\t'||s[i]=='\n') continue;
     if(s[i]>='a'&&s[i]<='z') continue;
     return 0;
  }
  return 1;
}

int Strhandler::isupperchar(char *s)
{
  int i,n;
  if(s==0) return -1;
  n=strlen(s);
  if(n==0) return -1;
  for(i=0;i<n;i++)
  {
     if(s[i]=='\0') break;
     if(s[i]=='\r'||s[i]=='\t'||s[i]=='\n') continue;
     if(s[i]>='A'&&s[i]<='Z') continue;
     return 0;
  }
  return 1;
}

int Strhandler::ischar(char *s)
{
  int i,n;
  if(s==0) return -1;
  n=strlen(s);
  if(n==0) return -1;

  for(i=0;i<n;i++)
  {
     if(s[i]=='\0') break;
     if(s[i]=='\r'||s[i]=='\t'||s[i]=='\n') continue;
     if(s[i]>='A'&&s[i]<='Z') continue;
     if(s[i]>='a'&&s[i]<='z') continue;
     return 0;
  }
  return 1;
}

int Strhandler::isfloat(char *s)
{
  int i,n;
  if(s==0) return -1;
  n=strlen(s);
  if(n==0) return -1;


  for(i=0;i<n;i++)
  {
     if(s[i]=='\0') break;
     if(s[i]>='0'&&s[i]<='9') continue;
     if(s[i]=='.') continue;
     if(s[i]=='\r'||s[i]=='\t'||s[i]=='\n') continue;
     return 0;
  }
  return 1;
}

int Strhandler::isinteger(char *s)
{
  int i,n;
  if(s==0) return -1;
  n=strlen(s);
  if(n==0) return -1;

  for(i=0;i<n;i++)
  {
     if(s[i]=='\0') break;
     if(s[i]=='\r'||s[i]=='\t'||s[i]=='\n') continue;
     if(s[i]>='0'&&s[i]<='9') continue;
     return 0;
  }
  return 1;
}

int Strhandler::isempty(char *file,char *s)
{
  int i,n;
  if(file==0) return -1;
  n=strlen(file);
  if(n==0) return -1;
  for(i=0;i<n;i++)
  {
    if(strchr(s,file[i])) continue;
    return 0;
  }
  return 1;
}

int Strhandler::isstdchar(char *f)
{
   int i;
   if(f==0) return 0;
   i=0;
   while(f[i])
   {
     if(isstdchar(f[i])==0) return 0;
     i++;
   }
   return 1;
}
int Strhandler::isstdchar(char f)
{
   if(f>='a'&&f<='z') return 1;
   if(f>='A'&&f<='Z') return 1;
   if(f>='0'&&f<='9') return 1;     
   return 0;
}
int Strhandler::isgoodchar(char f)
{
   if(f>='a'&&f<='z') return 1;
   if(f>='A'&&f<='Z') return 1;
   if(f>='0'&&f<='9') return 1;   
   if(f=='`') return 1;
   if(f=='~') return 1;
   if(f=='!'||f=='@'||f=='#'||f=='$'||f=='%'||f=='^') return 1;
   if(f=='&'||f=='*'||f=='('||f==')'||f=='_'||f=='-') return 1;
   if(f=='['||f==']'||f=='{'||f=='}') return 1;
   if(f=='+'||f=='='||f=='|'||f=='\\'||f==':'||f==';') return 1;
   if(f=='\''||f=='"'||f=='<'||f==',') return 1;
   if(f=='>'||f=='.'||f=='?'||f=='/') return 1;
   return 0;   
}

int Strhandler::isemailendchar(char f)
{
    if(f>='a'&&f<='z') return 0;
    if(f>='A'&&f<='Z') return 0;
    if(f>='0'&&f<='9') return 0;
    if(f=='_'||f=='-'||f=='+') return 0;
    if(f=='.') return 0;
    if(f=='%') return 0;
    return 1;
}

char *Strhandler::getpart(char *file,int a,int b)
{
   int k,n,i;
   char *s;
   if(b-a<0) return 0;
   if(file==0) return 0;
   n=strlen(file);
   if(b>n) b=n;
   if(a>=n) return 0;
   s=new char[b-a+20];
   if(n==0) return 0;
   k=0;
   for(i=0;i<n;i++)
   {
     if(i>=a&&i<=b) s[k++]=file[i];
   }
   s[k]='\0';
   return s;
}

char *Strhandler::getstr(int n)
{
  char *s;
  s=new char[32];
  s[0]='\0';
  sprintf(s,"%i",n);
  return s;
}

char *Strhandler::getstr(float n)
{
  char *s;
  s=new char[32];
  s[0]='\0';
  sprintf(s,"%f",n);
  return s;
}
 
int Strhandler::strncasecmplast(char *s,char *s0,int n)
{
  int i1,i2,i;
  char a,b;
  i1=strlen(s);i2=strlen(s0);
  if(i1-n<0||i2-n<0) return 1;
  for(i=0;i<n;i++)
  {
    a=s[i1-i-1];b=s0[i2-i-1];
    if(a>='A'&&a<='Z')a=a+32;
    if(b>='A'&&b<='Z')b=b+32;
    if(a!=b) return 1;
  }
  return 0;
}

int Strhandler::mystrncasecmp(char *s,char *s0,int n)
{
  int i1,i2,i;
  char a,b;
  if(s==0&&s0) return -1;
  else if(s&&s0==0) return 1;
  i1=strlen(s);i2=strlen(s0);
  if(i1-n<0||i2-n<0) return 1;
  for(i=0;i<n;i++)
  {
    a=s[i];b=s0[i];
    if(a>='A'&&a<='Z')a=a+32;
    if(b>='A'&&b<='Z')b=b+32;
    if(a!=b) return 1;
  }
  return 0;
}

int Strhandler::strncmplast(char *s,char *s0,int n)
{
  int i1,i2,i;
  char a,b;
  i1=strlen(s);i2=strlen(s0);
  if(i1-n<0||i2-n<0) return 1;
  for(i=0;i<n;i++)
  {
    a=s[i1-i-1];b=s0[i2-i-1];
    if(a!=b) return 1;
  }
  return 0;
}

char **Strhandler::pairbytoken(char *str,char *tg)
{
  int n,i;
  char *s,t,*tmp;
  char **pairlist;
  int m,len;

  if(str==0||tg==0) return 0;
  len=strlen(tg);
  m=strlen(str)+1;
  pairlist=new char*[m];
  for(i=0;i<m;i++) pairlist[i]=0;

  n=0; tmp=str;
  while((s=strstr(tmp,tg)))
  {
    t=*s;
    *s='\0';
    pairlist[n++]=opnstr(tmp);
    *s=t;
    tmp=s+len;if(tmp[0]=='\0') break;
  }
  if(tmp[0]!='\0') pairlist[n++]=opnstr(tmp);
  n=0;
  for(i=0;i<m;i++)
  {
    if(pairlist[i]==0) continue;
    pairlist[n++]=pairlist[i];
  }
  for(i=n;i<m;i++) pairlist[i]=0;
  return pairlist;
}

char **Strhandler::pairbytokensimple(char *str,char *tg)
{
  int n,i;
  char *s,t,*tmp;
  char **pairlist;
  int m,len;

  if(str==0||tg==0) return 0;
  len=strlen(tg);
  m=strlen(str)+1;
  pairlist=new char*[m];
  for(i=0;i<m;i++) pairlist[i]=0;

  n=0; tmp=str;
  while((s=strstr(tmp,tg)))
  {
    t=*s;
    *s='\0';
    pairlist[n++]=strdup(tmp);
    *s=t;
    tmp=s+len;if(tmp[0]=='\0') break;
  }
  if(tmp[0]!='\0') pairlist[n++]=strdup(tmp);
  n=0;
  for(i=0;i<m;i++)
  {
    if(pairlist[i]==0) continue;
    pairlist[n++]=pairlist[i];
  }
  for(i=n;i<m;i++) pairlist[i]=0;
  return pairlist;
}



int Strhandler::stlen(char *s)
{
  if(s==0) return 0;
  else     return strlen(s);
}

char *Strhandler::getbetween(char *file,char *a,char *b)
{
  char *s,*t;
  if(file==0) return 0;
  if(a==0||b==0) return file;
  s=strstr(file,a);
  if(s==0) { file=strdel(file); return 0;}
  s+=strlen(a);
  t=strstr(s,b);
  if(t==0) { file=strdel(file); return 0;}
  *t='\0';
  t=strdup(s);
  file=strdel(file);
  return t;
}

char *Strhandler::getbetweenwithend(char *file,char *a,char *b)
{
  char *s,*t;
  if(file==0) return 0;
  if(a==0||b==0) return file;
  s=strstr(file,a);
  if(s==0) { file=strdel(file); return 0;}
  s+=strlen(a);
  t=strstr(s,b);
  if(t==0) 
  { 
     t=opnstr(s);
     file=strdel(file);
     file=t;return file;
  }
  *t='\0';
  t=strdup(s);
  file=strdel(file);
  return t;
}
char *Strhandler::getbetweeninclude(char *file,char *a,char *b)
{
  char *s,*t;
  if(file==0) return 0;
  if(a==0||b==0) return file;
  s=strstr(file,a);
  if(s==0) { file=strdel(file); return 0;}
  t=strstr(s+strlen(a),b);
  if(t==0) { file=strdel(file); return 0;}
  t+=strlen(b);
  *t='\0';
  t=strdup(s);
  file=strdel(file);
  return t;
}
int Strhandler::getposition(char *file,char *tag)
{
  int n,i;

  if(file==0||tag==0) return -1;

  n=strlen(file);
  if(n==0) return -1;

  for(i=0;i<n;i++)
  { 
      if(file+i==tag) return i;
  }
  return -1;
}

char *Strhandler::insert(char *file,int n,char *tmp)
{
  int m,m0;
  char *t;
  if(file==0||tmp==0) return file;
  m=strlen(file); m0=strlen(tmp);
  if(m0==0) return file;

  if(m<=n) return strcat(file,tmp);
  else if(n<0)
  {
     t=opnstr(m+m0+10);
     strcat(t,tmp);
     strcat(t,file);      
     file=strdel(file); 
     file=t;
     return file;
  }

  file=(char *)realloc(file,m+m0+10);
  t=opnstr(file+n); 
  file[n]='\0';
  strcat(file,tmp);
  if(t) strcat(file,t);
  t=strdel(t);
  return file;
}
 

char *Strhandler::insert(char *page,char *token,char *tmp)
{
  char *t;
  int n,i;
 
  if(page==0||token==0||tmp==0) return page;

  t=strstr(page,token);
  if(t==0) return page;
  t+=strlen(token);
  n=strlen(page);
  for(i=0;i<n;i++)
  {
     if(page+i==t) break;
  }
  if(i<n)
  {
     page=insert(page,i,tmp);
  }
  return page;
}

char *Strhandler::insertbefore(char *page,char *token,char *tmp)
{
  char *t;
  int n,i;
 
  if(page==0||token==0||tmp==0) return page;

  t=strstr(page,token);
  if(t==0) return page;
  //t+=strlen(token);
  n=strlen(page);
  for(i=0;i<n;i++)
  {
     if(page+i==t) break;
  }
  if(i<n)
  {
     page=insert(page,i-1,tmp);
  }
  return page;
}

char ** Strhandler::getdelfilebyline(char *file,char *str)
{
  int i,n,m;
  char **tt;

  tt=opnfilebyline(file);

  n=0;while(tt&&tt[n]) n++;
  if(n==0) return 0;

  m=0;
  for(i=0;i<n;i++)
  {
    if(strcmp(tt[i],str)==0)
    {
       tt[i]=strdel(tt[i]); continue;
    }
    tt[m++]=tt[i];
  }
  tt[m]=0;
  return tt;
}

int Strhandler::delfilebyline(char *file,char *str)
{
  int i,n,m;
  char **tt;

  tt=opnfilebyline(file);
 
  n=0;while(tt&&tt[n]) n++;
  if(n==0) return -1;

  m=0;
  for(i=0;i<n;i++) 
  {
    if(strcmp(tt[i],str)==0) 
    {
       tt[i]=strdel(tt[i]); continue;
    }
    tt[m++]=tt[i];
  }
  tt[m]=0;
  prnfilebyline(file,tt);
  tt=strdel(tt);
  return m;
}

int Strhandler::prnfile(char *file,char *str)
{
  FILE *fp;
  int i;


  fp=fopen(file,"w");

  if(fp==0) return 0;

  if(str==0) 
  {
    fclose(fp); 
    return 0;
  }

  i=fprintf(fp,"%s",str);

  fclose(fp);
  return i;
}

int Strhandler::prnfilebyline(char *file,char **str)
{
  FILE *fp;
  int i,n,j;

  j=0;
  fp=fopen(file,"w");
  if(fp==0) return 0;
  if(str==0) 
  {
   fclose(fp); 
   return 0;
  }
  n=0;while(str[n])n++;
  if(n==0) return n;
  for(i=0;i<n;i++)
  {
    j+=fprintf(fp,"%s\n",str[i]);
    //cout<<i<<" <br>"<<endl;
  } 
  fclose(fp);
  return j;
}

int  Strhandler::prnfilebyline(FILE *fp,char **str)
{
  int i,n,j;

  j=0;
  if(fp==0) return 0;
  if(str==0) return 0;
  n=0;while(str[n])n++;
  if(n==0) return n;
  for(i=0;i<n;i++)
  {
    j+=fprintf(fp,"%s\n",str[i]);
  }
  fflush(fp);
  return j;
}


int Strhandler::findid(char *s0,char *s)
{
  int j,m;
  char *t;
  
  if(s0==0||s==0) return -1;

  t=strstr(s0,s);
  if(t==0) return -1;

  m=strlen(s0);
  for(j=0;j<m;j++)
  {
     if(s0+j==t) return j;
  }
  return -1; 
}

int Strhandler::findid(char **s0,char *s)
{
  int n;
 
  if(s0==0||s==0) return -1;
  n=0;
  while(s0[n])
  {
    if(strstr(s0[n],s)==0) {n++;continue;}
    return n;
  }
  return -1;
}

char *Strhandler::cuttextbetween(char *file,char *first,char *last)
{
  int tot;
  int i,j,k,n,m;
  
  tot=0;
  while(1)
  {
  tot++;
  if(tot>100) return file;
  if(first==0||last==0||file==0) return file;  //does not exist

  i=findid(file,first);

  if(i==-1) return file;  //does not exist
  
  n=strlen(first);
  
  j=findid(file+i+n,last); 

  if(j==-1) return file;  //does not exist

  j+=i+n;

  n=strlen(last);

  j+=n-1;
  
  n=strlen(file);
  
  m=0;
  for(k=0;k<n;k++)
  {
    if(k>=i&&k<=j) continue;  //kick out
    file[m++]=file[k]; 
  }
  file[m]='\0';

  }

  //return file;
}


char *Strhandler::cuttextinclude(char *file, char *token)
{
  int i,j,k,n,m,tot;
  int all;

  all=0;
  while(1)
  {
     all++;
     if(all>100) return file;
     if(file==0||token==0) return file;

     m=findid(file,token);
  
     if(m==-1) return file; //does not exist token
  
     n=strlen(token);
  
     tot=strlen(file);
     
     for(i=m+n;i<tot;i++)
     {
        if(file[i]==' ') break;
     }
  
     for(j=m-1;j>=0;j--)
     {
        if(file[j]==' ') break;
     }
  
     m=0;
     for(k=0;k<tot;k++)
     {
        if(k>j&&k<i) continue;  //kick out
        file[m++]=file[k];
     }
     file[m]='\0';
  }
  //return file;
}

int Strhandler::pushfile(char *file,char *name)
{   
  char **rec,**rec0;
  int n,i;
  if(name==0||name[0]=='\0') return 0;
  if(file==0) return 0;

  rec=opnfilebyline(file);
  n=0; while (rec&&rec[n]) n++;   
  rec0=opnarray(n+10);  
  rec0[0]=opnstr(name);
  for(i=0;i<n;i++)
  {
    rec0[i+1]=rec[i]; rec[i]=0;
  }
  rec0[n+1]=0;
  prnfilebyline(file,rec0); 
  rec=strdel(rec); rec0=strdel(rec0);
  return 1;

}


int Strhandler::pushfile(char *file,char **name)
{   
  char **rec,**rec0;
  int n,i,n0;
  if(name==0||name[0]=='\0') return 0;
  if(file==0) return 0;

  rec=opnfilebyline(file);

  n=0; while (rec&&rec[n]) n++;   
  n0=0; while(name&&name[n0])n0++;

  rec0=opnarray(n+n0+10);  

  for(i=0;i<n0;i++) rec0[i]=opnstr(name[i]);  
  for(i=0;i<n;i++) { rec0[i+n0]=rec[i]; rec[i]=0; }
  rec0=shrink(rec0,n+n0+3);
  prnfilebyline(file,rec0); 
  rec=strdel(rec); rec0=strdel(rec0);
  return 1;

}


int Strhandler::attachfilebylineifnoexist(char *name,char *file)
{
  char **tt;
  int n;
  tt=opnfilebyline(name);
  n=ifexist(tt,file);
  if(n>-1) return 0;
  
  n=attachfilebyline(name,file);
  
  return 1;
}

int Strhandler::ifexist(char **s,char *line)
{
  int n,i;

  n=0;while(s&&s[n])n++;
  if(n==0) return -1; 
  if(line==0) return -1;
  for(i=0;i<n;i++)
  {
    if(strcmp(s[i],line)==0) return i;
  }
  return -1;
}

int Strhandler::attachfile(char *name,char *file)
{
  FILE *fp;

  if(name==0||name[0]=='\0') return 0;
  if(file==0) return 0;

  fp=fopen(name,"a");
  if(fp==0) return 0;
  fprintf(fp,"%s",file);
  fclose(fp);
  return 1;
}

int Strhandler::attachfilebyline(char *name,char *file)
{
  FILE *fp;

  if(name==0||name[0]=='\0') return 0;
  if(file==0) return 0;

  fp=fopen(name,"a");
  if(fp==0) return 0;
  fprintf(fp,"%s\n",file);
  fclose(fp);
  return 1;
}



int Strhandler::attachfile(char *name,char **file)
{
  FILE *fp;
  int n;
  if(name==0||name[0]=='\0') return 0;
  if(file==0) return 0;

  fp=fopen(name,"a");
  if(fp==0) return 0;
  n=0;
  while(file[n]) 
  {
    fprintf(fp,"%s\n",file[n]);
    n++;
  }
  fclose(fp);
  return 1;
}

char *Strhandler::strstrlast(char *s0,char *token)
{
   char *s,*t;
   int n;
   if(s0==0) return 0;
   if(token==0) return 0;
   n=strlen(token);
   s=s0;t=0;
   while((s=strstr(s,token))!=NULL)
   {
     t=s;s+=n;
   }
   return t;
}

int Strhandler::viewablechar(char *s)
{
   int j,i,n;

   i=0;n=0;

   while(s[i]!='\0')
   {
      j=isgoodchar(s[i]);
      if(j==1) n+=2;
      else     n++; 
      i++;
   }
   return n;
}

char *Strhandler::translate(char **ctoe,char *cat)
{
   Strhandler cc;
   char *filetmp,**name;

   int i,n;
    
   if(ctoe==0) return 0;
   n=0;while(ctoe[n])n++;
   filetmp=0;
   for(i=0;i<n;i++)
   {
     name=cc.pairbytoken(ctoe[i],":");
     if(name==0||name[0]==0||name[1]==0) continue;
     if(strcmp(name[0],cat)==0)
     {
        filetmp=cc.opnstr(name[1]);
        filetmp=cc.clearendchar(filetmp,"/n/r/t ");
        name=cc.strdel(name);
        if(strcmp(filetmp,"<!>")==0)
        {
           filetmp=cc.strdel(filetmp);
        }
        return filetmp;
     }
     name=cc.strdel(name);
   }
   return 0;
}

char **Strhandler::opnarray(int n)
{
   int i;

   char **s;

   s=new char*[n];

   for(i=0;i<n;i++) s[i]=0;
   return s;
}

int **Strhandler::opnarrayint(int n)
{
   int i;

   int **s;

   s=new int*[n];

   for(i=0;i<n;i++) s[i]=0;
   return s;
}


char *Strhandler::gethttptime(long int n)
{
   char name[100];
   char **tt,*s,*t;
   struct tm *tme;

   tme=gmtime(&n); 
   s=asctime(tme);
   t=opnstr(s);
   t=clearendchar(t,"\r\t\n ");
   tt=pairbytoken(t," ");
   sprintf(name,"%s, %s %s %s %s GMT",tt[0],tt[2],tt[1],tt[4],tt[3]);
   t=strdel(t);tt=strdel(tt);
   t=opnstr(name);
   return t;
}

char *Strhandler::gettimeymd(int n)
{
   struct tm *s;
   time_t tloc;
   int y,m,d;
   char *t;


   tloc=time(0)+n;
   s = gmtime (&tloc);
   
   if(s==0) return 0;
   d=s->tm_mday;
   m=s->tm_mon+1;
   y=s->tm_year+1900;
 
   t=new char[100];

   sprintf(t,"/%i/%i/%i",y,m,d);
   return t;
}

char *Strhandler::gettimeymd_s(int n)
{
   struct tm *s;
   time_t tloc;
   int y,m,d;
   char *t;


   tloc=time(0)+n;
   s = gmtime (&tloc);
  
   if(s==0) return 0;
   d=s->tm_mday;
   m=s->tm_mon+1;
   y=s->tm_year+1900;

   t=new char[100];

   sprintf(t,"%i_%i_%i",y,m,d);
   return t;
}

int Strhandler::istimeright(int day)
{
   struct tm *s;
   time_t tloc;
   int d;

   tloc=time(0);
   s = gmtime (&tloc);

   if(s==0) return 0;
   d=s->tm_mday;

   if(d==day)         return 1;
   else               return 0;
}

int Strhandler::istimeright(int mon,int day)
{
   struct tm *s;
   time_t tloc;
   int m,d;

   tloc=time(0);
   s = gmtime (&tloc);
 
   if(s==0) return 0;
   d=s->tm_mday;
   m=s->tm_mon+1;

   if(m==mon&&d==day) return 1;
   else               return 0;
}

int Strhandler::istimeright(int year,int mon,int day)
{
   struct tm *s;
   time_t tloc;
   int y,m,d;

   tloc=time(0);
   s = gmtime (&tloc);
 
   if(s==0) return 0;
   d=s->tm_mday;
   m=s->tm_mon+1;
   y=s->tm_year+1900;
 
   if(m==mon&&d==day&&y==year) return 1;
   else                        return 0;
}

char *Strhandler::getfile(char *file)
{
   char *s;
   char *t;

   if(file==0) return 0;
      
   s=strstrlast(file,"/");
   if(s==0) 
   {
     t=opnstr(file);
     return t;
   }
   else
   {
     t=opnstr(s+1);
     return t;
   }
}

char *Strhandler::getpath(char *file)
{
   char *s,*t;

   if(file==0) return 0;
   s=strstrlast(file,"/");
   if(s==0)
   {
     return 0;
   }
   else
   {
     *s='\0';
     t=opnstr(file);
     *s='/';  
     return t;
   }
}

int Strhandler::filefound(char *file)
{
   FILE *fp;
   fp=fopen(file,"r");
   if(fp==0) return 0;
   else      
   {
      fclose(fp);
      return 1; 
   }
}
void Strhandler::setzero(char *s,int n)
{
  int i;

  for(i=0;i<n;i++)
  {
     s[i]='\0';
  }

}

char ** Strhandler::clearendchar(char **t,char *s) {

	int n=gettotnum(t);

	int ii=0;

	for(int i=0;i<n;i++) {
		t[i]=clearendchar(t[i],s);
		if(t[i]) t[ii++]=t[i];
	}

	t[ii]=0;

	return t;
}

char *Strhandler::getsegment(char *s,int start,int end) {

	if(s==0) return 0;

	int n=strlen(s);

	if(n<start) return 0;

	if(n<end) return strdup(s+start);

	char *ss=new char[n];

	for(int i=start;i<end;i++) {

		ss[i-start]=s[i];
	}

	ss[end-start]='\0';

	return ss;

}

void Strhandler::clearbadascii(char *s){

	if(s==0) return;

	int n=strlen(s);

	int i;
	int j=0;
	for(i=0;i<n;i++) {
		char c=s[i];
		if(c>='a'&&c<='z') {
			s[j++]=c;
			continue;
		}
		else if(c>='A'&&c<='Z') {
			s[j++]=c;
                        continue;
		}
		else if(c>='0'&&c<='9') {
			s[j++]=c;
                        continue;
		}
		else if(strchr("-.*",c)) {
			s[j++]=c;
                        continue;
		}
		else if(c=='\'') {
			s[j++]='*';
                        continue;
		}
	}
	s[j]='\0';
}
