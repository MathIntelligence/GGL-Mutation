#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "source.h"
#define SKIP " \t\n:;,"
 
Rcs::Rcs(char *filename)

{
int i,j,k;
char line[256],*pointer,*placeholder;
char *workingdirectory;
FILE *fp;
char rpath[20][1000],rname[20][100];

number=0;name=0;path=0;

// Look in currentdirectory first.  Otherwise HOME directory.

for(i=0;i<20;i++)
{
for(j=0;j<1000;j++) rpath[i][j]='\0';
for(j=0;j<100;j++) rname[i][j]='\0';
}

if((fp=fopen(filename,"r"))==NULL) {
	char *t=getenv("JACKALDIR");
	if((workingdirectory=getenv("JACKALDIR"))!=NULL) {
		sprintf(line,"%s",workingdirectory);
	}
	else {
		printf("Could not find resource file: %s\n",filename); 
		cerr<<"check the environment path:JACKALDIR"<<endl;
        	cerr<<"something wrong in specifying jackal.dir"<<endl;
		exit(0); 
	}
}
  
if(fp==NULL) if((fp=fopen(line,"r"))==NULL)
{ printf("Could not find resource file: %s\n",filename); exit(0); }
while(fgets(line,256,fp)!=NULL)
{
  placeholder=line;  
  while(strchr(SKIP,*placeholder)) placeholder++; //Skip whitespace.
  pointer=placeholder;
  if(strchr(placeholder,':'))
  {  
     number++;
     while(*placeholder!=':') placeholder++;  
     *placeholder='\0';
     strcpy(rname[number-1],pointer);    
     pointer=++placeholder;
     while(strchr(SKIP,*pointer)) pointer++; 
  }
  placeholder=pointer;
  while(!strchr(SKIP,*placeholder)) placeholder++;    
  *placeholder='\0';
  if(*pointer=='\0') continue;
  placeholder=rpath[number-1];
  if(*placeholder=='\0') strcpy(placeholder,pointer);
  else
  {
    while(*placeholder!='\0'||*(placeholder+1)!='\0')  placeholder++;
    strcpy(placeholder+1,pointer);
  }

}
for(i=0;i<number;i++) {path=new char*[number];name=new char*[number];}
for(i=0;i<number;i++)
{
j=strlen(rname[i]);
name[i]=new char[j+1];
strcpy(name[i],rname[i]);
name[i][j]='\0';
pointer=rpath[i];
placeholder=pointer;
j=1;
while(*placeholder!='\0'||*(placeholder+1)!='\0'){placeholder++; j++;}
j=j+2;
path[i]=new char[j];
for(k=0;k<j;k++)path[i][k]=rpath[i][k];
path[i][j-1]='\0';
}

char *mlib=RCS["key"];
if(mlib==0||strstr(mlib,"XIANG2002")==0) {
cerr<<endl<<"the key in jackal.dir is not correct"<<endl;
cerr<<endl<<"key:"<<mlib<<endl<<endl;
cerr<<"contact the author at: zx11@columbia.edu"<<endl;
cerr<<endl;
cerr<<"please tell us who you are and the intention of using this software"<<endl;
cerr<<"thank you for your interests in Jackal\n";
cerr<<endl;
exit(0);
}
else {

	cerr<<endl;
	cerr<<"If you detect any bugs, please download the newest version at: "<<endl;
	cerr<<"http://trantor.bioc.columbia.edu/~xiang/jackal"<<endl;
	cerr<<"We also appreciate your reporting the bug to the author"<<endl;
	cerr<<endl;

}

}

Rcs::~Rcs()
{
  int i;
  for(i=0;i<number;i++) {delete [] name[i]; delete [] path[i];}
  delete [] name; delete [] path;
}

char *Rcs::operator[](char *s)
{
int i,j;
j=strlen(s);
for(i=0;i<number;i++)
if(strncmp(name[i],s,j)==0) return path[i];
return 0;
}
