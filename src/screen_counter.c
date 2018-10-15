#include <stdio.h>

screen_counter_(i,n,m)

int  *i;
int  *n;
int  *m;

{

float percent;
int   ipercent;

percent= *i;

percent=(1.- percent/ *n)*100.;
ipercent=percent;

if ( *m>0)
  {
  printf(" Doing step %d of %d - %d percent remaining - %d iterations\r", *i, *n, ipercent, *m);
  }

else
  {
  printf(" Doing %d of %d - %d percent remaining \r", *i, *n, ipercent);
  }

fflush(0);

}
