#include <stdlib.h>

int main(int argc, char **argv)
{
  int i, num;

  for(i=0;i<9;i++)
    printf("Shfit by 1 << %d = %d\n", i, 1<<i);

  num = 1<<0;
  num += 1<<3;
  num += 1<<8;
  printf("num = %d\n", num);
  for(i=0;i<9;i++)
    printf("%d >> %d div 2 = %d\n", num, i, (num >> i) % 2);

  return 0;
}
