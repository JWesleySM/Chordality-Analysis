#include <stdlib.h>

void foo(int a, int b, int c, int d, int e, int r1, int r2, int r3){
  c=r3;
  a=r1;
  b=r2;
  d=0;
  e=a;
  while(e>0){
    d=d+b;
    e=e-1;
  }
  r1=d;
  r3=c;
  return;
}
