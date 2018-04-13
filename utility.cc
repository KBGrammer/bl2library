#include <sys/stat.h> 
#include "utility.h"

long file_size(char * filename)
{
  FILE *pfile;
  struct stat fstat;
  stat(filename,&fstat);
  return fstat.st_size;
}
