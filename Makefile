DepthHist: DepthHist.c
	cc -g -O3 -I./samtools-1.3/htslib-1.3 DepthHist.c ./samtools-1.3/htslib-1.3/libhts.a  -lpthread -lz -o $@
pUC18.fa.sa: pUC18.fa
	bwa index pUC18.fa
fr.sam: fr.fa pUC18.fa.sa
	bwa mem -p pUC18.fa fr.fa > fr.sam
clean:
	rm DepthHist pUC18.fa.sa fr.sam 
