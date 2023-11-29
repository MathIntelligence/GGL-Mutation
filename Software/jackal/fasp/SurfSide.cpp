#include"source.h"

SurfSide::SurfSide()
{
}

SurfSide::~SurfSide()
{
}

void SurfSide::readpdb(char *file) {

        Strhandler te;
        FILE *fp;

        fp=fopen(file,"r");

        if(fp==0) {
                cerr<<"could not open file:"<<file<<endl;
                return;
        }

        char line[1000];

	int nnn=0;
        while(fgets(line,1000,fp)!=NULL) {

                char c='1';

                //if(c=='0') c='1';

                line[4]='\0';

                Pdb pdb;

                te.lowercase(line);

                pdb.read(line,c);

                cerr<<line<<" "<<c<<endl;
                Res *r,*rr;

                if(pdb.chn==0) continue;
		
		pdb.chn->header();
		pdb.setflg(1);
		pdb.chn->surface(10,1.4);
		pdb.chn->buildhbond();
		pdb.chn->setdsspstr();

                for(r=pdb.chn->res;r;r=r->next) {
                         
			nnn++;
			if(r->ambgt==0) continue;
			if(strchr("GA",r->name)) continue;
			cout<<line<<" "<<nnn<<" "<<r->name<<r->id<<":"<<r->sec<<" "<<1-r->bury(4,100);
			cout<<" "<<r->ambgtrmsd(0)<<" "<<r->ambgtrmsd(1)<<" "<<r->ambgtrmsd(2);
			cout<<" "<<r->ambgtrmsd(3)<<" "<<r->ambgtrmsd(4)<<endl;				 	                         
                }
        }

}
