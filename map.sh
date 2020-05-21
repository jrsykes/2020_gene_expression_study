#!/bin/bash
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=40gb
#SBATCH --ntasks=6
#SBATCH --output=/projects/sykesj/StdOut/R-%x.%j-map.out
#SBATCH --error=/projects/sykesj/StdOut/R-%x.%j-map.err


SPECIES=$1
SRR=$2
SEX=$3
LAYOUT=$4
DUEL_LAYOUT=$5
WD=$6



if [ "$DUEL_LAYOUT" == YES ]
then
	PAIRED_BUSCO_SCORE=$(sed '8q;d' "$WD"/analyses/"$SPECIES"/busco/BUSCO_out_"$SPECIES"_PAIRED.txt | awk -F[CS] '{print $2}' | sed 's/[^0-9]*//g')
	SINGLE_BUSCO_SCORE=$(sed '8q;d' "$WD"/analyses/"$SPECIES"/busco/BUSCO_out_"$SPECIES"_SINGLE.txt | awk -F[CS] '{print $2}' | sed 's/[^0-9]*//g')

	if [ "$PAIRED_BUSCO_SCORE" -ge "$SINGLE_BUSCO_SCORE" ]
	then
		BEST_TRANS_IDX=PAIRED_"$SPECIES".idx
		touch "$WD"/analyses/"$SPECIES"/kallisto/mapped_to_PAIRED_idx
	
	elif [ "$PAIRED_BUSCO_SCORE" -lt "$SINGLE_BUSCO_SCORE" ]
	then
		BEST_TRANS_IDX=SINGLE_"$SPECIES".idx
		touch "$WD"/analyses/"$SPECIES"/kallisto/mapped_to_SINGLE_idx
	
	else
		touch "$WD"/analyses/"$SPECIES"/kallisto/"$SRR"_FailedToMap
		exit 1
	fi
fi



kallisto_map () {

	mkdir "$WD"/analyses/"$SPECIES"/kallisto/"$SRR"
	#mkdir /scratch/"$WD"/map_"$SRR"

	if [ $LAYOUT == 'PAIRED' ]
	then
		if [ $DUEL_LAYOUT == 'YES' ]
		then
			TRANS_IDX="$BEST_TRANS_IDX"
		else
			TRANS_IDX=PAIRED_"$SPECIES".idx
		fi

		kallisto quant -t 16 -i "$WD"/analyses/"$SPECIES"/kallisto/"$TRANS_IDX" -o "$WD"/analyses/"$SPECIES"/kallisto/"$SRR" \
			-b 100 "$WD"/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/"$SEX"/"$SRR"\_1.fq "$WD"/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/"$SEX"/"$SRR"\_2.fq


	elif [ $LAYOUT == 'SINGLE' ]
	then
		if [ $DUEL_LAYOUT == 'YES' ]
		then
			TRANS_IDX="$BEST_TRANS_IDX"
		else
			TRANS_IDX=SINGLE_"$SPECIES".idx
		fi

		READ_LENGTH=$(awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("%f\n",m);}'  "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/$SEX/$SRR\_s.fq)
		SD=$(awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("%f\n",sq/n-m*m);}' "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/$SEX/$SRR\_s.fq)

		kallisto quant -t 16 -i "$WD"/analyses/"$SPECIES"/kallisto/"$TRANS_IDX" -o "$WD"/analyses/"$SPECIES"/kallisto/"$SRR" -b 100 \
			--single -l "$READ_LENGTH" -s "$SD" "$WD"/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/"$SEX"/"$SRR"\_s.fq

	fi

	#rsync -a /scratch/"$WD"/map_"$SRR"/"$SRR" "$WD"/analyses/"$SPECIES"/kallisto && rm -rf /scratch/"$WD"/map_"$SRR"



	#ln -s "$WD"/analyses/"$SPECIES"/kallisto/"$SRR" "$WD"/analyses/"$SPECIES"/kallisto/kal_results/kal_files/
	#rm -rf "$WD"/analyses/"$SPECIES"/kallisto/"$SRR"/"$SRR"

}

kallisto_map $SPECIES $SRR $SEX $LAYOUT 