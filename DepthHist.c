#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <htslib/sam.h>

int64_t read_sam_and_fill_depth_buffer(htsFile*, bam_hdr_t*, int**);
void write_depths_as_wig(FILE*, bam_hdr_t*, int**);
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
  read_sam_and_fill_depth_buffer(htsf, header_p, depth_buffer);
  write_depths_as_wig(stdout, header_p, depth_buffer);
  sam_close(htsf);
}
int64_t
read_sam_and_fill_depth_buffer(htsFile*htsf, bam_hdr_t*header_p, int**depth_buffer)
{
  bam1_t *r1;
  bam1_t *r2;
  int retv1, retv2;
  r1=bam_init1();
  r2=bam_init1();
  retv1 = sam_read1(htsf, header_p, r1);
  fprintf(stdout, "%d:%d:%d:%d:%d\n", r1->core.tid, r1->core.pos, r1-> core.qual, r1->core.mtid, r1->core.mpos);
  retv2 = sam_read1(htsf, header_p, r2);
  fprintf(stdout, "%d:%d:%d:%d:%d\n", r2->core.tid, r2->core.pos, r2-> core.qual, r2->core.mtid, r2->core.mpos);
  return 0;
}
void
write_depths_as_wig(FILE*stdout, bam_hdr_t*header_p, int**depth_buffer)
{
}

