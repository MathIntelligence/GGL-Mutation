#ifndef _Contact
#define _Contact

class Atm;
class Contact
{
public:

Contact(int);
Contact();
~Contact();
void setengarray(int);
void setatmarray(int);
Atm  **atm;
float *distance;
float *energy;
float total;
int num;
int flag;
Contact *next;
};

#endif
