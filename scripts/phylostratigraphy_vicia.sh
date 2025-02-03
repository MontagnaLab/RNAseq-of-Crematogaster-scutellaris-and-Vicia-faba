cd comparative_genomics

echo -e "transcript\tstratum\texpression"

for transcript in $(grep ">" ../annotations/vicia/vicia.Trinity.fasta.transdecoder.cds | awk -F "." '{print $1}' | sed "s/_i[0-9]*//" | tr -d ">");
	
	do

                expression=$(grep $transcript ../abundances/vicia/vicia_deseq_gene/*_genes.lst | awk -F "_" '{print $4}')

                if [ -z "${expression}" ]; then expression="-"; fi

		line=$(grep $transcript vicia/Orthogroups.txt)

			if echo $line | grep -q -e GCA_019776745; then 
				echo -e "$transcript\tSPE\t$expression";

			elif echo $line | grep -q -e GCF_001865875 -e GCA_019845555 -e GCF_001605985 -e GCF_016808335 -e GCF_000313855; then
                                echo -e "$transcript\tANG\t$expression";

			elif echo $line | grep -q -e GCF_000003745 -e GCF_004153795 -e GCF_000188115 -e GCA_003112345 -e GCA_008009225 -e GCF_001879475 -e GCA_022085475; then
                                echo -e "$transcript\tEUD\t$expression";

                        elif echo $line | grep -q -e GCF_013401445 -e GCF_000001735 -e GCF_006381635 -e GCF_000493195 -e GCF_005239225 -e GCF_002906115 -e GCF_000346465; then
                                echo -e "$transcript\tROS\t$expression";

                        elif echo $line | grep -q -e GCF_004799145 -e GCA_014851425; then
                                echo -e "$transcript\tFAB\t$expression";

                        elif echo $line | grep -q -e GCF_003473485 -e GCF_000331145 -e GCF_000004515 -e GCF_020283565 -e GCF_001865875 -e GCF_000499845 -e GCA_003370565 -e GCA_004329165 -e GCF_003086295; then
                                echo -e "$transcript\tPAP\t$expression";
			else 
				echo -e "$transcript\tspecies-specific\t$expression";
			fi

done
