#ifndef _ICOSAHEDRON
#define _ICOSAHEDRON


#define NVER 2600
#define NFCE 5200
#define NEDGE 7800

class Icosahedron {
public:

  Icosahedron(int);
  ~Icosahedron();
  int divide(int,int);
  void trimid();
  void face(int);

  int nv,nf,ne;
  int **edg;
  int **edg2;
  int **fc;
  int *nedg;
  int *vncl;
  Vector **fcemd;
  int **fcedg;
  Vector *ver;
} ;

#endif


