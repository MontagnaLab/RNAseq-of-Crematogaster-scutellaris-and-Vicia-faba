echo -e "hit\tevalue\tqcover\ttranscript\tlogFC_AB_AD\tpadj_AB_AD\tlogFC_AC_AD\tpadj_AC_AD\tlogFC_AD_AD\tpadj_AD_AD\tlogFC_AB_CT\tpadj_AB_CT\tlogFC_AC_CT\tpadj_AC_CT\tlogFC_AD_CT\tpadj_AD_CT\tAD\tCT"

for i in $(awk -F "_i" '{print $1}' $1 | grep "TRI" | tr -d ">");

do 

	hit=$(grep $i $1 | awk '{print $2}'); if [[ -z $hit ]]; then hit="NA"; fi 
	evl=$(grep $i $1 | awk '{print $3}'); if [[ -z $evl ]]; then evl="NA"; fi	
	qcv=$(grep $i $1 | awk '{print $4}'); if [[ -z $qcv ]]; then qcv="NA"; fi

	expAB_AD=$(grep $i abundances/crema/crema_deseq2_gene/RSEM_crema.filtered.gene.counts.matrix.A_AD_vs_B_AD.DESeq2.DE_results | awk '{print $7"\t"$11}');
	if [[ -z $expAB_AD ]] ; then expAB_AD="NA\tNA";fi
	expAC_AD=$(grep $i abundances/crema/crema_deseq2_gene/RSEM_crema.filtered.gene.counts.matrix.A_AD_vs_C_AD.DESeq2.DE_results | awk '{print $7"\t"$11}');
	if [[ -z $expAC_AD ]]; then expAC_AD="NA\tNA";fi
	expAD_AD=$(grep $i abundances/crema/crema_deseq2_gene/RSEM_crema.filtered.gene.counts.matrix.A_AD_vs_D_AD.DESeq2.DE_results | awk '{print $7"\t"$11}'); 
	if [[ -z $expAD_AD ]]; then expAD_AD="NA\tNA";fi
	mod_AD=$(grep $i abundances/crema/crema_WGCNA_gene/AD_module_*_genes.lst | awk -F "_" '{print $5}')
	if [[ -z $mod_AD ]]; then mod_AD="NA";fi

	expAB_CT=$(grep $i abundances/crema/crema_deseq2_gene/RSEM_crema.filtered.gene.counts.matrix.A_CT_vs_B_CT.DESeq2.DE_results | awk '{print $7"\t"$11}');
	if [[ -z $expAB_CT ]]; then expAB_CT="NA\tNA";fi
	expAC_CT=$(grep $i abundances/crema/crema_deseq2_gene/RSEM_crema.filtered.gene.counts.matrix.A_CT_vs_C_CT.DESeq2.DE_results | awk '{print $7"\t"$11}');
	if [[ -z $expAC_CT ]]; then expAC_CT="NA\tNA";fi
	expAD_CT=$(grep $i abundances/crema/crema_deseq2_gene/RSEM_crema.filtered.gene.counts.matrix.A_CT_vs_D_CT.DESeq2.DE_results | awk '{print $7"\t"$11}'); 
	if [[ -z $expAD_CT ]]; then expAD_CT="NA\tNA";fi
	mod_CT=$(grep $i abundances/crema/crema_WGCNA_gene/CT_module_*_genes.lst | awk -F "_" '{print $5}')
        if [[ -z $mod_CT ]]; then mod_CT="NA";fi
	
	echo -e "$hit\t$evl\t$qcv\t$i\t$expAB_AD\t$expAC_AD\t$expAD_AD\t$expAB_CT\t$expAC_CT\t$expAD_CT\t$mod_AD\t$mod_CT";
	
done
