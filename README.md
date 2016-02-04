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
    DepthHist [-d depth_threshold] \
        [-n non_reporting_margin] \
        [-m min_mapq] \
        [-i min_insert] [-a max_insert] \
        [-s sam_file] [-o output] [wig files]

This program reads a sam_file and find pairs of reads having both at least min_mapq.
Calculate the span of the pair reads and if the span size is between min_insert and max_insert
accept as a good pair and increment the depth at each positions as depth.
The output is a standard uncompressed text representation of [wig format](https://genome.ucsc.edu/goldenpath/help/wiggle.html), which 
can be converted to a more efficient binary format using wigToBigWig (https://genome.ucsc.edu/goldenpath/help/bigWig.html).
The output defaults to stdout if not specified.
While reporting the depths to wig format output, points with depth lower than
depth_threshold will be written to stderr.
The low depth points can be gathered as low depth regions with
range_compress.rb.

# Example work flow
The following example assumes bash grammer, to redirect.

    bwa index ref.fa
    bwa mem -t 32 -p ref.fa mp.fq > mp.sam
    java -Xmx10G -XX:ParallelGCThreads=3 -jar path/to/picard.jar \
        SortSam I=mp.sam O=mp.bam \
        SO=coordinate CREATE_INDEX=true MAX_RECORDS_IN_RAM=5000000
    java -Xmx10G -XX:ParallelGCThreads=3 -jar path/to/picard.jar \
        CollectInsertSizeMetrics I=mp.bam \
        H=mp.hist O=mp.metrics HISTOGRAM_WIDTH=40000 MINIMUM_PCT=0.3

read mp.metrics and choose appropriate parameter (should be automated, but not yet)

    DepthHist -d 3 -n 10000 -m 40 -i 10000 -a 40000 -s mp.sam -o mp.wig -l mp.low_depth_points
    paste <(fatt name ref.fa) <(fatt len ref.fa) > ref.sizes
    wigToBigWig mp.wig ref.sizes mp.bw
    ruby range_compress.rb mp.low_depth_points > mp.low_depth_regions

    DepthHist -d 0 -n 7000 -m 40 -i 7000 -a 14000 -s mp2.sam -o mp2.wig -l mp2.low_depth_points
    samtools view -HS mp2.sam > samhead
    DepthHist -d 3 -n 7000 -s samhead -o mpc.wig -l mpc.low_depth_points mp.wig mp2.wig

Different library may have different valid size range. The depth from different libraries processed
with different parameters may be summed after counting each library, rather than
processing at once.

# BUILD
    wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
    tar jxvf samtools-1.3.tar.bz2
    make -C samtools-1.3/htslib-1.3
    make
