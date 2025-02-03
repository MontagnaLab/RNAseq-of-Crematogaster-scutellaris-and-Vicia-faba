# script to identify contaminat contigs from nr+blastp result
# $1 - blastp input 
# $2 - contaminants contigs output
# $3 - target taxonomy from which blast hits must derive (e.g. Pancrustacea or Viridiplantae)
# $4 - proportion of blast hits must derive from target taxonomy (e.g. 0.8 or 0.5)

for i in $(awk -F "\t" '{print $1}' $1 | sort -u);

        do

	a=$(grep $i $1 | awk -F "\t" '{print $(NF-1)" "$NF}' | sort -g | head -500 | awk '{print $NF}');

        echo $a | sed "s/ /\n/g" | sed "s/;/\n/g" > $1"tmp.taxon.lst";

 	cat $1"tmp.taxon.lst"

        taxonkit lineage $1"tmp.taxon.lst"  --data-dir dbs/taxdump > $1"tmp.lineage.lst" ;

	cat $1"tmp.lineage.lst"

        tot=$(wc -l $1"tmp.lineage.lst" | awk '{print $1}');

        bound=$(echo "$tot * $4" | bc | awk -F "." '{print $1}')

	if ! [[ $(grep -c $3 $1"tmp.lineage.lst") -ge bound ]] 
	        then echo $i >> $2;
        fi;

done

if [[ ! -f "$2" ]]; then touch $2; fi
