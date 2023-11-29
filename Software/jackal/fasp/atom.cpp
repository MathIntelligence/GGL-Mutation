#include "stdafx.h"
#include <ctype.h>
#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include "constants.h"
#include "vector3d.h"
#include "color.h"
#include "atom.h"
#include "structure.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


using namespace Troll;



Atom::Atom(Residue *r,char *n,float b) : residue(r), name(n,r), bfactor(b)

{

charge=0.0;
radius=0.0;

}


int Atom::Is(char *p)

{
int i;

if(residue->description==NULL) throw(TE_UnknownAtom,this);
if(residue->description->atom.find(name)==residue->description->atom.end())
  throw(TE_UnknownAtom,this);

AtomDescription *ad=residue->description->atom[name];

for(i=0;i<ad->properties.size();i++)
  if(strcmp(ad->properties[i],p)==0) return 1;

return 0;

}


float Atom::Area(unsigned long flag)

{
AtomDescription *ad;

switch(flag)
  {
  case A_Current:
  return area;

  case A_Buried:
  if(residue->description->atom.find(name)==residue->description->atom.end()) return 0.0;
  ad=residue->description->atom[name];
//  return ad->aref-area;

  case A_Extended:
  if(residue->description->atom.find(name)==residue->description->atom.end()) return 0.0;
  ad=residue->description->atom[name];
//  return ad->aref;

  case A_PercentBuried:
  if(residue->description->atom.find(name)==residue->description->atom.end()) return 0.0;
  ad=residue->description->atom[name];
//  return (1.0-area/ad->aref)*100.0;

  default:
  return 0.0;
  }

}


int AtomName::operator==(char *n)

{
ResidueDescription *rd=residue->description;
std::map<char*,char*,ltstr>::iterator aa;

if(!strcmp(name,n)) return 1;
for(aa=rd->alias.begin();aa!=rd->alias.end();aa++)
  if( !strcmp((*aa).first,n) && !strcmp((*aa).second,name) ) return 1;

return 0;

}

int AtomName::operator==(AtomName& n)

{

return (*this)==n.name;

}


int AtomName::operator=(char *n)

{
int i;
char *tn=n;

while( *tn==' ') *tn++;

strcpy(name,tn);

i=strlen(name)-1;
tn=&name[i];

while( *tn==' ' ) tn--;
*(++tn)='\0';

return 1;

}


AtomName::AtomName(char *n,Residue *r)

{
int i=0;
char t[5];
char *tn=n;

(*this)=n;

residue=r;

}




::ostream& Troll::operator<<(::ostream& os,Atom& at)

{
os << "ATOM  "
   << setw(5) << at.param 
   << " "
   << at.name
   << " "
   << setw(3) << at.residue->name
   << " "
   << at.residue->chain->name[0]
   << at.residue->rid
   << "   "
   << setprecision(3) << setw(8) << at.pos.x
   << setw(8) << at.pos.y
   << setw(8) << at.pos.z
   << setprecision(2) << setw(6) << 1.00
   << setw(6) << at.bfactor
   << "      " << at.residue->chain->name
   << endl;

return os;

}



::ostream& Troll::operator<<(ostream& os,AtomName& an)

{
int l;
char aname[5];

os.setf(ios::fixed);

l=strlen(an.name);
if(l==4) strcpy(aname,an.name);
else if(l==3)
  {
  if((an.name[l-1]>='0' && an.name[l-1]<='9') || an.name[l-1]=='*')
    sprintf(aname," %s",an.name);
  else sprintf(aname,"%s ",an.name);
  }
else if(l==2) sprintf(aname," %s ",an.name);
else if(l==1) sprintf(aname," %s  ",an.name);
else sprintf(aname,"ERR ");

os << aname;

return os;

}



