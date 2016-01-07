pUC18.fa.sa: pUC18.fa
	bwa index pUC18.fa

fr.sam: fr.fa pUC18.fa.sa
	bwa mem -p pUC18.fa fr.fa > fr.sam
