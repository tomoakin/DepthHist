# DepthHist
Calculate clone depth histogram in wig format from SAM/BAM

Mate-pair information is used to connect contigs into scaffold in genome assembly.
When there are a lot of mate-pair data and there is few mates spanning a certain
region that could be a indicator of misassembly.

"View As Pairs" capability of  IGV [http://www.broadinstitute.org/igv/AlignmentData]
is a good way to view it. However, to scan a large genome, I prefer to have a 
histogram like data. 
A depth plot, sometimes called coverage plot, showing read depth of each position
is a hint, but here I want to have a plot of clone depth that include the unsequenced
region.

# Requirements
This program uses HTSlib to parse SAM format data.
In ihe input sam file, paired read should be in a cluster of
reads with the same name. Therefore, the input should be
just the output of an alignment program like bwa, or
name sorted.

# Usage
DepthHist [-d depth_threshold] [-n non_reporting_margin] [-m min_mapq] [-i min_insert] [-a max_insert] [-s sam_file] [-o output]

# Example work flow
    bwa mem -p mp.fq  >  mp.sam
    DepthHist -d 3 -n 10000 -m 40 -i 10000 -a 40000 -s mp.sam -o mp.wig 2> mp.low_cov_points
    ruby range_compress.rb mp.low_cov_points > mp.low_cov_regions

