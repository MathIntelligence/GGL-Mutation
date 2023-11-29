#ifndef _ObjectBox
#define _ObjectBox

class ObjectBox {

public:
  ObjectBox();
 ~ObjectBox();
  int ne,max;
  int *element;
  void remove(int);
  void add(int);
};
#endif
