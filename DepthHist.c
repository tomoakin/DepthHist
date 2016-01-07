#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>

int
main(int argc, char** argv)
{
  int n_target;
  htsFile *htsf;
  bam_hdr_t * header_p;
  htsf = hts_open("fr.sam", "r");
  if(!htsf){
    fputs("sam/bam file open failed\n", stderr);
    exit(EXIT_FAILURE);
  }
  header_p = sam_hdr_read(htsf);
  if(!header_p){
    fputs("no header\n", stderr);
    exit(EXIT_FAILURE);
  }
  sam_close(htsf);  
}
