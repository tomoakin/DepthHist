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
/*  fprintf(stdout, "n_targets: %d\n", header_p->n_targets); */
  depth_buffer = malloc(sizeof(int*) * header_p->n_targets);
  if(!depth_buffer){fputs("memory allocation failed", stderr);exit(EXIT_FAILURE);}
  for(i=0;i<header_p->n_targets;i++){
/*    fprintf(stdout, "%s\t%d\n", header_p->target_name[i], header_p->target_len[i]); */
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
  bam1_t *t;
  int threshold = 30;
  int count = 0;
  int retv1, retv2;
  r1=bam_init1();
  r2=bam_init1();
  retv1 = sam_read1(htsf, header_p, r1);
  while(retv1 == 0){
    count ++;
/*    fprintf(stdout, "%s:%d:%d:%d:%d:%d\n", bam_get_qname(r1), r1->core.tid, r1->core.pos, r1-> core.qual, r1->core.mtid, r1->core.mpos);*/
    if(r1->core.qual < threshold) { /* more condition may come until a good mapping is found */
      retv1 = sam_read1(htsf, header_p, r1);
      continue;
    }
    retv2 = sam_read1(htsf, header_p, r2);
    while(retv2 == 0){ /* seeking matching record that pairs r1 */
      count ++;
/*      fprintf(stdout, "%s:%d:%d:%d:%d:%d\n", bam_get_qname(r2), r2->core.tid, r2->core.pos, r2-> core.qual, r2->core.mtid, r2->core.mpos);*/
      if(strcmp(bam_get_qname(r1),bam_get_qname(r2))!=0){
        /* r2 is a new read. Thus, we need to make r1 point to this and start to seek again */
       t=r2; r1=r2; r2=t;
       break;
      }
      /* now we have two records of single fragment */
      if( r2-> core.qual >= threshold &&r1->core.mtid == r2->core.tid && r1->core.mpos == r2->core.pos){
        int i,f,t;
        if(r1->core.pos < r2->core.pos){
          f = r1->core.pos;
          t = r2->core.pos + bam_cigar2rlen(r2->core.n_cigar, bam_get_cigar(r2));
        }else{
          f = r2->core.pos;
          t = r1->core.pos + bam_cigar2rlen(r1->core.n_cigar, bam_get_cigar(r1));
        }
/*        fprintf(stderr,"fill target %d: from %d to %d\n", r1->core.tid, f, t); */
        for(i=f; i<t; i++){
          depth_buffer[r1->core.tid][i] += 1;
        }
        retv1 = sam_read1(htsf, header_p, r1);
        break;
      }
      retv2 = sam_read1(htsf, header_p, r2);
    }
  }
  return count;
}
void
write_depths_as_wig(FILE*out, bam_hdr_t*header_p, int**depth_buffer)
{
  int i,j;
  fprintf(out, "track type=wiggle_0\n");
  for(i=0;i<header_p->n_targets;i++){
    fprintf(out, "fixedStep chrom=%s start=1%d\n", header_p->target_name[i]);
    for(j=0; j < header_p-> target_len[i]; j++){
      fprintf(out, "%d\n", depth_buffer[i][j]);
    }
  }
}

