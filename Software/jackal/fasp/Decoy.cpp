#include"source.h"

Decoy::Decoy()
{
 native=0;
 decoy=0;
 algn=1;
}

Decoy::~Decoy()
{
  if(decoy) delete decoy;
}

void Decoy::randomskip(int n) {
	
	srandom(n);
}

void Decoy::createdecoyfly(int angle,int number,Chn *nat) {
 
	Chn *t;
	Stralg alg;
	int n=0;
        alg.flag=2;
	while(1) {

		t=randmz(native,angle);	
                cout<<"the rmsd:"<<alg.superimpose(t,nat,1)<<" "<<alg.superimpose(t,native,1)<<endl;
		delete t;t=0;
		if(n++>number) break;
	}
	 
}

void Decoy::createdecoy(int angle,int number) {

        Chn *t;
        Stralg alg;
        int n=0;
        alg.flag=2;
        Chn *s=0;
        for(s=decoy;s&&s->next;s=s->next);
	
        while(1) {

                t=randmz(native,angle);
                if(algn) cout<<"the rmsd:"<<alg.superimpose(t,native,1)<<endl;
		if(s==0) {
                        decoy=t;
                        s=t;
                }
                else {
                        s->next=t;
                        s=s->next;
                }
                if(n++>number) break;
        }
}


Chn *Decoy::randmz(Chn *chn,int angle) {

	Res *r;
	Atm *a;
	Rotate rot;

	Chn *s=new Chn(chn);
	s->configure();
	if(angle<=0) return s;
	for(r=s->res;r;r=r->next) 
	for(a=r->atm;a;a=a->next) {
			
 		int n=random()%angle;	
		
		if(a->tatm->rotate==0) continue;
		if(a->tatm->id>=4) continue;

		rot.rotate(a,n,1);
 	} 
	return s; 
}

void Decoy::write(char *s) {

	Chn *ss;

	char line[1000];
	int n=0;
	for(ss=decoy;ss;ss=ss->next) {
		ss->addresid(1);
		ss->addatmid0(1);
		sprintf(line,"%s_%i.pdb",s,n++);
		ss->write(line);
	}

}
