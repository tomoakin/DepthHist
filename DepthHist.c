#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <htslib/sam.h>

int getopt(int argc, char * const argv[],
                  const char *optstring);
extern char *optarg;
extern int optind, opterr, optopt;

const char * const track_type_string = "track type=wiggle_0\n";


typedef struct _pair_params{
  int min_valid_mapq;
  int min_proper_insert;
  int max_proper_insert;
} pair_params;
typedef struct _region_params{
  int depth_threshold;
  int non_reporting_margin;
} region_params;

int read_wig_and_add_to_depth(const char *wigfilename, bam_hdr_t* header_p, int **depth_buffer);
int64_t read_sam_and_fill_depth_buffer(htsFile*, bam_hdr_t*, int**, pair_params);
void write_depths_as_wig(FILE*, bam_hdr_t*, int**, region_params, FILE*);
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
  int c;
  FILE *out=stdout;
  FILE *log=NULL;
  char*endptr;

  while ((c = getopt(argc, argv, "d:n:m:i:a:s:o:l:")) != -1){
    switch (c){
    case 'd':
      region_p.depth_threshold = strtol(optarg,&endptr,0);
      break;
    case 'n':
      region_p.non_reporting_margin = strtol(optarg,&endptr,0);
      break;
    case 'm':
      pair_p.min_valid_mapq = strtol(optarg,&endptr,0);
      break;
    case 'i':
      pair_p.min_proper_insert = strtol(optarg,&endptr,0);
      break;
    case 'a':
      pair_p.max_proper_insert = strtol(optarg,&endptr,0);
      break;
    case 's':
      htsf = hts_open(optarg, "r");
      if(!htsf){
        fputs("sam file open failed\n", stderr);
        exit(EXIT_FAILURE);
      }
      break;
    case 'o':
      out = fopen(optarg, "w");
      if(!out){
        fputs("output file open failed\n", stderr);
        exit(EXIT_FAILURE);
      }
      break;
    case 'l':
      log = fopen(optarg, "w");
      if(!out){
        fputs("log file open failed\n", stderr);
        exit(EXIT_FAILURE);
      }
      break;
    default:
      usage();
    }
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
  { /* for each file in command line option, parse as wig and sum into depth_buffer */
    int i;
    for(i = optind; i< argc; i++)
      read_wig_and_add_to_depth(argv[i], header_p, depth_buffer);
  }
  write_depths_as_wig(out, header_p, depth_buffer, region_p, log);
  sam_close(htsf);
}

int
is_proper_pair(bam1_t *r1, bam1_t *r2, pair_params param)
{
  if( r2-> core.qual >= param.min_valid_mapq && 
      r1-> core.qual >= param.min_valid_mapq && 
      r1-> core.tid == r2->core.tid &&
      r1-> core.mtid == r2->core.tid && 
      r1-> core.mpos == r2->core.pos &&
      r1-> core.isize >= param.min_proper_insert && 
      r1 -> core.isize <= param.max_proper_insert) 
    return 1;
  return 0;
}

void
fill_depth_buffer(bam1_t *r1, bam1_t *r2, int**depth_buffer)
{
  int i,f,t;
  f = r1->core.pos;
  t = r2->core.pos + bam_cigar2rlen(r2->core.n_cigar, bam_get_cigar(r2));
  for(i=f; i<t; i++){
    depth_buffer[r1->core.tid][i] += 1;
  }
}

int
read_wig_and_add_to_depth(const char *wigfilename, bam_hdr_t* header_p, int **depth_buffer)
{
  int retv = 0;
  int linebufsize=4096; /*1 page; */
  char* linebuf,*retp;
  FILE*wigfh;
  char * scaff_name;
  int start; int step;
  linebuf = malloc(linebufsize);
  if(!linebuf){fputs("memory allocation failed", stderr);exit(EXIT_FAILURE);}
  scaff_name = malloc(linebufsize);
  if(!scaff_name){fputs("memory allocation failed", stderr);exit(EXIT_FAILURE);}
  wigfh=fopen(wigfilename,"r");
  if(!wigfh){fprintf(stderr, "wig file open failed: %s\n", wigfilename);perror(NULL);exit(EXIT_FAILURE);}
  retp=fgets(linebuf,linebufsize,wigfh);
  if(!retp){fprintf(stderr, "wig file could not read: %s\n", wigfilename);perror(NULL);exit(EXIT_FAILURE);}
  if(strcmp(linebuf,track_type_string,linebufsize)){
    fprintf(stderr, "The first line do not match exactly: %s\n", wigfilename);
    fprintf(stderr, "in file : %s\n", linebuf);
    fprintf(stderr, "expected: %s\n", track_type_string);
    exit(EXIT_FAILURE);
  }
/*  track type=wiggle_0 */
/*fixedStep chrom=scaffold_blahblah start=1 step=1 */
  while(retp=fgets(linebuf,linebufsize,wigfh)){
    int i,ch, scaff_id;
    sscanf(linebuf, "fixedStep chrom=%s start=%i step=%i", scaff_name, &start, &step);
    scaff_id=bam_name2id(header_p, scaff_name);
    for(i = start-1; i < header_p->target_len[scaff_id]; i+= step){
      ch = getc(wigfh); ungetc(ch,wigfh); 
      if(!isdigit(ch))break;
      retp=fgets(linebuf,linebufsize,wigfh);
      if(!retp){fprintf(stderr, "wig file readerr: %s\n", wigfilename);perror(NULL);exit(EXIT_FAILURE);}
      depth_buffer[scaff_id][i] += atoi(linebuf);
    }
  }
  free(linebuf);
  free(scaff_name);
  return retv;
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
  while(retv1 >= 0){
    count ++;
    if(r1->core.qual < param.min_valid_mapq) { /* more condition may come until a good mapping is found */
      retv1 = sam_read1(htsf, header_p, r1);
      continue;
    }
    retv2 = sam_read1(htsf, header_p, r2);
    while(retv2 >= 0){ /* seeking matching record that pairs r1 */
      count ++;
      if(strcmp(bam_get_qname(r1), bam_get_qname(r2))!=0){
        /* r2 is a new read. Thus, we need to make r1 point to this and start to seek again */
        t=r2; r2=r1; r1=t;
        break;
      }
      /* now we have two records of single fragment */
      if(is_proper_pair(r1,r2, param)){
        if(r1->core.pos < r2->core.pos){
          fill_depth_buffer(r1,r2,depth_buffer);
        }else{
          fill_depth_buffer(r2,r1,depth_buffer);
        }
        retv1 = sam_read1(htsf, header_p, r1);
        if(retv1 == 0)
        break;
      }
      retv2 = sam_read1(htsf, header_p, r2);
    }
    retv1 = retv2;
  }
  return count;
}
void
write_depths_as_wig(FILE*out, bam_hdr_t*header_p, int**depth_buffer, region_params param, FILE* log)
{
  int i,j;
  fprintf(out, track_type_string);
  for(i=0;i<header_p->n_targets;i++){
    fprintf(out, "fixedStep chrom=%s start=1 step=1\n", header_p->target_name[i]);
    for(j=0; j < header_p-> target_len[i]; j++){
      fprintf(out, "%d\n", depth_buffer[i][j]);
      if(log && j> param.non_reporting_margin && 
         j< header_p-> target_len[i] - param.non_reporting_margin &&
         depth_buffer[i][j] < param.depth_threshold)
        fprintf(log, "%s %i %i\n", header_p->target_name[i], j, depth_buffer[i][j]);
    }
  }
}
void usage()
{
  fputs("DepthHist [-d depth_threshold] [-n non_reporting_margin] [-m min_mapq] [-i min_insert] [-a max_insert] [-s sam_file] [-o output_wigfile] [-l low_depth_points_file] [wig files]", stderr);
}

