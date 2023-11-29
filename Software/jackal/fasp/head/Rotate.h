#ifndef _Rotate
#define _Rotate

class Rotate
{
 public:
 Rotate();
 ~Rotate();
 void hook(Res *);
 void hook(Res *from,int);
 void hook(Res *from,Res *to); //hook from cb
 void hook(Res *from,Res *to,int);
 void rotateomega(Atm *atm0, float angle,int n,int flg);
 void hookheader(Res *to);
 void link(Res *from,Res *to,int end,int flg);
 void rotate(Atm *, float angle); //rotate atom in this residue;
 void rotate(Atm *, float angle, int flg); //rotate atom in n residue;
 void rotate(float *,float *,float,float);
 void rotate(Atm *, float angle, int n,int flg); //rotate atom in n residue;
 int *rota;
 int flag;
};
#endif
