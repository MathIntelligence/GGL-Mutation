#include"../model/source.h"
Rcs RCS(".model");
Tres TRES;
main(int argc,char *argv[])
{
Build c;
Res *res; 
char line[200];
char chn;
TRES.flag=100;
TRES.read("tres");

if(argc<=3)chn='1';
else 
{
chn=argv[3][0];
if(chn=='0') chn='1';
}
strcpy(line,argv[1]);
c.pdb->read(line,chn);
printf("\n%s_%c\n",line,chn);;
for(res=c.pdb->chn->res;res;res=res->next)
{
   if(argv[2][0]=='0') 
   {
     res->dihedral(stdout);
   }
   else if(res->name==argv[2][0]) 
   {
     res->dihedral(stdout);
   }
}
exit(0);
}
