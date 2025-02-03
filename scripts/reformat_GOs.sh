# script to reformat PANNZER2 and EggNOG-mapper output for TopGO
# $1 - PANNZER2 predictions.lst input file
# $2 - EggNOG-mapper annotations.tsv input file
# $3 - output file

if [[ -f $3 ]]; then echo "annotation already exists"; exit; fi

for i in $(awk '{print $1}' $1 | sort -u);

        do

	GOs=$(grep $i $1 | awk '{print $3}');
        refGOs=$(echo $GOs | sed "s/ /, GO:/g");
        if ! [ -z "$refGOs" ];
                then echo -e "$i\t$refGOs" | sed "s/\t/\tGO:/" | sed 's/_i[0-9]*\.p[0-9]*//g' >> $3;
                else echo -e $i  >> $3;
        fi;

done;

for i in $(awk '{print $1}' $2 | grep TRINITY | awk 'BEGIN{FS=OFS="_"}NF{NF-=1};1' | sort -u); 
	
        do

        GOs=$(grep $i $2 | awk -F "\t" '{print $10}' | sed 's/,/, /g'); 
        if ! [ "$GOs" == - ]; 
                then if grep -q $i $3;
                        then 
			ln=$(grep -n $i $3 | awk -F ":" '{print $1}')
			exGOs=$(grep $i $3 | awk -F "\t" '{print $2}');
			nwGOs=$(echo $GOs $exGOs | tr ' ' '\n' | sort | uniq | tr '\n' ' ' | sed -e 's/[[:space:]]*$//');
			sed -i "${ln}d" $3;
			echo -e "$i\t$nwGOs" >> $3;
                        else 
			echo -e "$i\t$GOs" >> $3;
                fi;
	fi;

done;

