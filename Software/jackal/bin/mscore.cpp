#include"../fasp/head/source.h"
Rcs RCS("jackal.dir");
Tres TRES;
main(int argc,char *argv[])
{
Pdb  pdb;
TRES.flag=102;
TRES.read("tres");
TRES.setdipolecharge();
TRES.setnewengcoeff(0,0.5);
TRES.readmore("back_small_rotamer","side_small_rotamer");
//TRES.set3backrotamer("backbone");
TRES.smt=1;
TRES.prepareengtable();
TRES.readdistpopular("bound_nearnext1168.dist","internal");
TRES.readdistpopular("bound_hbondall1168.dist","hbond");
pdb.read(argv[1],'1');
pdb.chn->header();
pdb.chn->addhatoms();
pdb.chn->configure();
pdb.chn->write("a");
Disc disc;
disc.pdb=&pdb;
disc.ready();
//disc.setup(pdb.chn->res,10);
Res *r=pdb.chn->isres0(0);
disc.setrange();
disc.setup(pdb.chn->res,1000,1);
disc.cutoff=10;
disc.setdipoleascharge();
disc.setrange();
//disc.write(stdout);
//disc.logg=1;
strcpy(disc.force,"udD");
float d=disc.clash(pdb.chn->res,10000);
int k=0;
for(k=0;k<200;k++) {
if(disc.energy[k]==0) continue;
cerr<<char(k)<<":"<<disc.energy[k]<<endl;
}
cerr<<"the total energy is: "<<d<<endl;
}

