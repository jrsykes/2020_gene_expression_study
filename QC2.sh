SPECIES=$1
LAYOUT=$2

if [ $LAYOUT == 'PAIRED' ]
then
	TRIMMED_LIBS=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/PAIRED/male/*.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/PAIRED/female/*.fq) ; do readlink -f $file; done | paste -sd " " - )
elif [ $LAYOUT == 'SINGLE' ]
then
	TRIMMED_LIBS=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/SINGLE/male/*.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/SINGLE/female/*.fq \
		/projects/sykesj/analyses/$SPECIES/trimmomatic/PAIRED/male/*.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/PAIRED/female/*.fq) ; do readlink -f $file; done | paste -sd " " - )
fi



/home/sykesj/software/FastQC/fastqc --outdir /projects/sykesj/analyses/$SPECIES/fastqc2 $TRIMMED_LIBS \
	&& multiqc /projects/sykesj/analyses/$SPECIES/fastqc/ -o /projects/sykesj/analyses/$SPECIES/fastqc/ && rm -f /projects/sykesj/analyses/$SPECIES/fastqc/SRR* \
	&& multiqc /projects/sykesj/analyses/$SPECIES/fastqc2/ -o /projects/sykesj/analyses/$SPECIES/fastqc2/ && rm -f /projects/sykesj/analyses/$SPECIES/fastqc2/SRR*







