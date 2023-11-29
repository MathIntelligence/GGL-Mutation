#include"../model/source.h"
Rcs RCS(".model");
Tres TRES;
main(int argc,char *argv[])
{
Domain  dom;
FILE *fp;
Build c;
char line[200],line0[200];
Res *rrr[1000000];
int nnn;
Fish cc(&c);
Pdb *pdb; 
TRES.flag=100;
TRES.read("tres");
Res *r,*r1,*r2;
int i,j,k,m,n,kk;
float x,y,z; 
Atm *a,*a1;

c.pdb->read("1myf_mod",'1');

cc.header(c.pdb->chn);

i=0;

r1=0;
for(r=c.pdb->chn->res;r;r=r->next)
{
if(r->name=='L')
{
  r->id=i++;
  r->id0=r->id;
  if(r1==0) r1=r;
  else {dom.link(r1,r,r->id);r1=r;}
  r->write(stdout);
}
}


exit(0);
}
