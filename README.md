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

Because the mapping data comes as SAM or BAM file and we need to
process a large volume of data, I try 
to implement in C using samtools/HTSlib.
