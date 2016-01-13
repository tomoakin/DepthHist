#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <htslib/sam.h>

typedef struct _pair_params{
  int min_valid_mapq;
  int min_proper_insert;
  int max_proper_insert;
} pair_params;
typedef struct _region_params{
  int depth_threshold;
  int non_reporting_margin;
} region_params;

int64_t read_sam_and_fill_depth_buffer(htsFile*, bam_hdr_t*, int**, pair_params);
void write_depths_as_wig(FILE*, bam_hdr_t*, int**, region_params);
void usage();
int
main(int argc, char** argv)
{
  int n_target;
  int i;
  pair_params pair_p = {40, 10000, 40000};
  region_params region_p = {5, 11000};
  htsFile *htsf;
  bam_hdr_t * header_p;
  int **depth_buffer;
  if(argc != 2){
    usage();
    exit(EXIT_FAILURE);
  }
  htsf = hts_open(argv[1], "r");
  if(!htsf){
    fputs("sam file open failed\n", stderr);
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
  read_sam_and_fill_depth_buffer(htsf, header_p, depth_buffer, pair_p);
  write_depths_as_wig(stdout, header_p, depth_buffer, region_p);
  sam_close(htsf);
}
int64_t
read_sam_and_fill_depth_buffer(htsFile*htsf, bam_hdr_t*header_p, int**depth_buffer, pair_params param)
{
  bam1_t *r1;
  bam1_t *r2;
  bam1_t *t;
  int count = 0;
  int retv1, retv2;
  r1=bam_init1();
  r2=bam_init1();
  retv1 = sam_read1(htsf, header_p, r1);
  while(retv1 == 0){
    count ++;
#if 0
    fprintf(stderr, "r1 %s\t%d:%d:%d:%d:%d\n", bam_get_qname(r1), r1->core.tid, r1->core.pos, r1-> core.qual, r1->core.mtid, r1->core.mpos);
#endif 
    if(r1->core.qual < param.min_valid_mapq) { /* more condition may come until a good mapping is found */
      retv1 = sam_read1(htsf, header_p, r1);
      continue;
    }
    retv2 = sam_read1(htsf, header_p, r2);
    while(retv2 == 0){ /* seeking matching record that pairs r1 */
      count ++;
#if 0
      fprintf(stderr, "r2 %s\t%d:%d:%d:%d:%d\n", bam_get_qname(r2), r2->core.tid, r2->core.pos, r2-> core.qual, r2->core.mtid, r2->core.mpos);
#endif 
      if(strcmp(bam_get_qname(r1),bam_get_qname(r2))!=0){
        /* r2 is a new read. Thus, we need to make r1 point to this and start to seek again */
        t=r2; r2=r1; r1=t;
        break;
      }
      /* now we have two records of single fragment */
      if( r2-> core.qual >= param.min_valid_mapq && r1->core.qual >= param.min_valid_mapq && 
          r1-> core.tid == r2->core.tid &&
          r1-> core.mtid == r2->core.tid && r1->core.mpos == r2->core.pos &&
          r1-> core.isize >= param.min_proper_insert && r1 -> core.isize <= param.max_proper_insert){
        int i,f,t;
        if(r1->core.pos < r2->core.pos){
          f = r1->core.pos;
          t = r2->core.pos + bam_cigar2rlen(r2->core.n_cigar, bam_get_cigar(r2));
        }else{
          f = r2->core.pos;
          t = r1->core.pos + bam_cigar2rlen(r1->core.n_cigar, bam_get_cigar(r1));
        }
#if 0
        fprintf(stderr,"fill target %d: from %d to %d\n", r1->core.tid, f, t);
#endif 
        for(i=f; i<t; i++){
          depth_buffer[r1->core.tid][i] += 1;
        }
        retv1 = sam_read1(htsf, header_p, r1);
        if(retv1 == 0)
#if 0
          fprintf(stderr, "r1 %s\t%d:%d:%d:%d:%d\n", bam_get_qname(r1), r1->core.tid, r1->core.pos, r1-> core.qual, r1->core.mtid, r1->core.mpos);
#endif 
        break;
      }
      retv2 = sam_read1(htsf, header_p, r2);
    }
    retv1 = retv2;
  }
  return count;
}
void
write_depths_as_wig(FILE*out, bam_hdr_t*header_p, int**depth_buffer, region_params param)
{
  int i,j;
  fprintf(out, "track type=wiggle_0\n");
  for(i=0;i<header_p->n_targets;i++){
    fprintf(out, "fixedStep chrom=%s start=1 step=1\n", header_p->target_name[i]);
    for(j=0; j < header_p-> target_len[i]; j++){
      fprintf(out, "%d\n", depth_buffer[i][j]);
      if(j> param.non_reporting_margin && 
         j< header_p-> target_len[i] - param.non_reporting_margin &&
         depth_buffer[i][j] < param.depth_threshold)
        fprintf(stderr, "%s %i %i\n", header_p->target_name[i], j, depth_buffer[i][j]);
    }
  }
}
void usage()
{
  fputs("DepthHist [-d depth_threshold] [-n non_reporting_margin] [-m min_mapq] [-i min_insert] [-a max_insert] [-s sam_file] > wigfile", stderr);
}

