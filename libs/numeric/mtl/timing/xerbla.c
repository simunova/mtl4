//#include "f2c.h"

#include <stdio.h>

typedef long int integer;

int xerbla_(char *srname, integer *info)
{

  printf("** On entry to %6s, parameter number %2i had an illegal value\n",
	 srname, *info);
  return 0;
}
