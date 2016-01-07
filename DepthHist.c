#include <stdio.h>
#include "bam.h"

int
main(int argc, char** argv)
{
  tamFile samf;
  samf = sam_open("fr.sam");  
  if(!samf){
    exit(EXIT_FAILURE);
  }
  sam_close(samf);  
}
