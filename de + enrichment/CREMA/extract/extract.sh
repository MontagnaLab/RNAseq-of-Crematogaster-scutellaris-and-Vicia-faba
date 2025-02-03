echo transcript AD_AvsB AD_AvsC AD_AvsD AD_BvsC AD_BvsD AD_CvsD CT_AvsB CT_AvsC CT_AvsD CT_BvsC CT_BvsD CT_CvsD

for i in $(blastp -query $1 -db crema.Trinity.fasta.transdecoder.ol.pep -evalue $2 -qcov_hsp_perc $3 -outfmt "6 qseqid sseqid evalue qcovs" | awk '{print $2}' | awk -F "_i" '{print $1}' | sort -u);

	do

#
	
	if grep -q $i ../AD/AD_AvsB_UP_genes_lst.txt; then AD_AvsB=UP
	elif grep -q $i ../AD/AD_AvsB_DN_genes_lst.txt; then AD_AvsB=DN		
	else AD_AvsB=NA
	fi
	
	if grep -q $i ../AD/AD_AvsC_UP_genes_lst.txt; then AD_AvsC=UP
	elif grep -q $i ../AD/AD_AvsC_DN_genes_lst.txt; then AD_AvsC=DN		
	else AD_AvsC=NA
	fi
	
	if grep -q $i ../AD/AD_AvsD_UP_genes_lst.txt; then AD_AvsD=UP
	elif grep -q $i ../AD/AD_AvsD_DN_genes_lst.txt; then AD_AvsD=DN		
	else AD_AvsD=NA
	fi

#	

	if grep -q $i ../AD/AD_BvsC_UP_genes_lst.txt; then AD_BvsC=UP
	elif grep -q $i ../AD/AD_BvsC_DN_genes_lst.txt; then AD_BvsC=DN		
	else AD_BvsC=NA
	fi
	
	if grep -q $i ../AD/AD_BvsD_UP_genes_lst.txt; then AD_BvsD=UP
	elif grep -q $i ../AD/AD_BvsD_DN_genes_lst.txt; then AD_BvsD=DN		
	else AD_BvsD=NA
	fi
		
	if grep -q $i ../AD/AD_CvsD_UP_genes_lst.txt; then AD_CvsD=UP
	elif grep -q $i ../AD/AD_CvsD_DN_genes_lst.txt; then AD_CvsD=DN		
	else AD_CvsD=NA
	fi
	
#	
	
	if grep -q $i ../CT/CT_AvsB_UP_genes_lst.txt; then CT_AvsB=UP
	elif grep -q $i ../CT/CT_AvsB_DN_genes_lst.txt; then CT_AvsB=DN		
	else CT_AvsB=NA
	fi
	
	if grep -q $i ../CT/CT_AvsC_UP_genes_lst.txt; then CT_AvsC=UP
	elif grep -q $i ../CT/CT_AvsC_DN_genes_lst.txt; then CT_AvsC=DN		
	else CT_AvsC=NA
	fi
	
	if grep -q $i ../CT/CT_AvsD_UP_genes_lst.txt; then CT_AvsD=UP
	elif grep -q $i ../CT/CT_AvsD_DN_genes_lst.txt; then CT_AvsD=DN		
	else CT_AvsD=NA
	fi
	
#	

	if grep -q $i ../CT/CT_BvsC_UP_genes_lst.txt; then CT_BvsC=UP
	elif grep -q $i ../CT/CT_BvsC_DN_genes_lst.txt; then CT_BvsC=DN		
	else CT_BvsC=NA
	fi
	
	if grep -q $i ../CT/CT_BvsD_UP_genes_lst.txt; then CT_BvsD=UP
	elif grep -q $i ../CT/CT_BvsD_DN_genes_lst.txt; then CT_BvsD=DN		
	else CT_BvsD=NA
	fi
		
	if grep -q $i ../CT/CT_CvsD_UP_genes_lst.txt; then CT_CvsD=UP
	elif grep -q $i ../CT/CT_CvsD_DN_genes_lst.txt; then CT_CvsD=DN		
	else CT_CvsD=NA
	fi	

#	
	
echo $i $AD_AvsB $AD_AvsC $AD_AvsD $AD_BvsC $AD_BvsD $AD_CvsD $CT_AvsB $CT_AvsC $CT_AvsD $CT_BvsC $CT_BvsD $CT_CvsD
	
done