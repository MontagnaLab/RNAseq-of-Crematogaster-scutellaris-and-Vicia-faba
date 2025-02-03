cd comparative_genomics

echo -e "transcript\tstratum\tmodule"

for transcript in $(grep ">" ../annotations/crema/crema.Trinity.fasta.transdecoder.pep | awk '{print $1}' | awk -F "." '{print $1}' | awk -F "_" 'BEGIN { OFS = "_" } NF{NF-=1};1' | tr -d ">");
	
	do

		module=$(grep $transcript ../abundances/crema/crema_WGCNA_gene/$1_module* | awk -F "_" '{print $5}')

		if [ -z "${module}" ]; then module=0; fi

		line=$(grep $transcript crema/Orthogroups.txt)

			if echo $line | grep -q -e GCA_000507165 -e GCA_000376725; then 
				echo -e "$transcript\tPAL\t$module";

			elif echo $line | grep -q -e GCA_003018175 -e GCF_000696155 -e GCF_021461395; then
                                echo -e "$transcript\tNEO\t$module";

			elif echo $line | grep -q -e GCF_000001215 -e GCF_000005575 -e GCF_014905235 -e GCF_009731565 -e GCF_000002335 -e GCF_000390285 -e GCF_000696795 -e GCF_005508785; then
                                echo -e "$transcript\tHOL\t$module";

                        elif echo $line | grep -q -e GCF_910589235 -eGCF_000612105 -e GCF_021155765 -e GCF_015476425 -e GCF_003254395 -e GCF_009193385; then
                                echo -e "$transcript\tHYM\t$module";

                        elif echo $line | grep -q -e GCF_003651465 -e GCA_001045655 -e GCF_010583005 -e GCF_000217595 -e GCF_003672135 -e GCF_002006095; then
                                echo -e "$transcript\tFOR\t$module";

                        elif echo $line | grep -q -e CF_000143395 -e GCF_016802725 -e GCF_000956235 -e GCF_000949405 -e GCF_003070985; then
                                echo -e "$transcript\tMYR\t$module";
			else 
				echo -e "$transcript\tspecies-specific\t$module";
			fi

done
