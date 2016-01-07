#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>

int
main(int argc, char** argv)
{
  int n_target;
  int i;
  htsFile *htsf;
  bam_hdr_t * header_p;
  int **depth_buffer;
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
  fprintf(stdout, "n_targets: %d\n", header_p->n_targets);
  depth_buffer = malloc(sizeof(int*) * header_p->n_targets);
  if(!depth_buffer){fputs("memory allocation failed", stderr);exit(EXIT_FAILURE);}
  for(i=0;i<header_p->n_targets;i++){
    fprintf(stdout, "%s\t%d\n", header_p->target_name[i], header_p->target_len[i]);
    depth_buffer[i] = calloc(header_p->target_len[i], sizeof(int));
    if(!depth_buffer[i]){fputs("memory allocation failed", stderr);exit(EXIT_FAILURE);}
  }
  sam_close(htsf);
}
