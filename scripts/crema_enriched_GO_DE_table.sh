awk -F "\t" '{print $1"\t"$2}' enrichment/GO_crema/crema_A_*BP_enrichment.tsv | sort -u | tail -n +2 > tot_GOs.tmp

echo "GOterm,AB_CT_UP,AB_CT_DN,AC_CT_UP,AC_CT_DN,AD_CT_UP,AD_CT_DN,AB_AD_UP,AB_AD_DN,AC_AD_UP,AC_AD_DN,AD_AD_UP,AD_AD_DN,description";

while read line;

	do

	j=$(echo $line | awk '{print $1}')
	d=$(echo $line | awk '{$1=""; print $0}')
#AB_CT_UP
	if grep -q $j enrichment/GO_crema/crema_A_CT_vs_B_CT_UP_BP_enrichment.tsv;
		then AB_CT_UP=$(grep $j enrichment/GO_crema/crema_A_CT_vs_B_CT_UP_BP_enrichment.tsv | awk -F "\t" '{print $6}');
		else AB_CT_UP="NA";
	fi
#AB_CT_DN
        if grep -q $j enrichment/GO_crema/crema_A_CT_vs_B_CT_DN_BP_enrichment.tsv;
                then AB_CT_DN=$(grep $j enrichment/GO_crema/crema_A_CT_vs_B_CT_DN_BP_enrichment.tsv | awk -F "\t" '{print $6}');
                else AB_CT_DN="NA";
        fi
#AC_CT_UP
        if grep -q $j enrichment/GO_crema/crema_A_CT_vs_C_CT_UP_BP_enrichment.tsv;
                then AC_CT_UP=$(grep $j enrichment/GO_crema/crema_A_CT_vs_C_CT_UP_BP_enrichment.tsv | awk -F "\t" '{print $6}');
                else AC_CT_UP="NA";
        fi
#AC_CT_DN
        if grep -q $j enrichment/GO_crema/crema_A_CT_vs_C_CT_DN_BP_enrichment.tsv;
                then AC_CT_DN=$(grep $j enrichment/GO_crema/crema_A_CT_vs_C_CT_DN_BP_enrichment.tsv | awk -F "\t" '{print $6}');
                else AC_CT_DN="NA";
        fi
#AD_CT_UP
        if grep -q $j enrichment/GO_crema/crema_A_CT_vs_D_CT_UP_BP_enrichment.tsv;
                then AD_CT_UP=$(grep $j enrichment/GO_crema/crema_A_CT_vs_D_CT_UP_BP_enrichment.tsv | awk -F "\t" '{print $6}');
                else AD_CT_UP="NA";
        fi
#AD_CT_DN
        if grep -q $j enrichment/GO_crema/crema_A_CT_vs_D_CT_DN_BP_enrichment.tsv;
                then AD_CT_DN=$(grep $j enrichment/GO_crema/crema_A_CT_vs_D_CT_DN_BP_enrichment.tsv | awk -F "\t" '{print $6}');
                else AD_CT_DN="NA";
        fi

##################################################

#AB_AD_UP
	if grep -q $j enrichment/GO_crema/crema_A_AD_vs_B_AD_UP_BP_enrichment.tsv;
                then AB_AD_UP=$(grep $j enrichment/GO_crema/crema_A_AD_vs_B_AD_UP_BP_enrichment.tsv | awk -F "\t" '{print $6}');
                else AB_AD_UP="NA";
        fi
#AB_AD_DN
	if grep -q $j enrichment/GO_crema/crema_A_AD_vs_B_AD_DN_BP_enrichment.tsv;
                then AB_AD_DN=$(grep $j enrichment/GO_crema/crema_A_AD_vs_B_AD_DN_BP_enrichment.tsv | awk -F "\t" '{print $6}');
                else AB_AD_DN="NA";
        fi
#AC_AD_UP
	if grep -q $j enrichment/GO_crema/crema_A_AD_vs_C_AD_UP_BP_enrichment.tsv;
                then AC_AD_UP=$(grep $j enrichment/GO_crema/crema_A_AD_vs_C_AD_UP_BP_enrichment.tsv | awk -F "\t" '{print $6}');
                else AC_AD_UP="NA";
        fi
#AC_AD_DN
	if grep -q $j enrichment/GO_crema/crema_A_AD_vs_C_AD_DN_BP_enrichment.tsv;
                then AC_AD_DN=$(grep $j enrichment/GO_crema/crema_A_AD_vs_C_AD_DN_BP_enrichment.tsv | awk -F "\t" '{print $6}');
                else AC_AD_DN="NA";
        fi
#AD_AD_UP
	if grep -q $j enrichment/GO_crema/crema_A_AD_vs_D_AD_UP_BP_enrichment.tsv;
                then AD_AD_UP=$(grep $j enrichment/GO_crema/crema_A_AD_vs_D_AD_UP_BP_enrichment.tsv | awk -F "\t" '{print $6}');
                else AD_AD_UP="NA";
        fi
#AD_AD_DN
	if grep -q $j enrichment/GO_crema/crema_A_AD_vs_D_AD_DN_BP_enrichment.tsv;
                then AD_AD_DN=$(grep $j enrichment/GO_crema/crema_A_AD_vs_D_AD_DN_BP_enrichment.tsv | awk -F "\t" '{print $6}');
                else AD_AD_DN="NA";
        fi

##################################################

	echo "$j,$AB_CT_UP,$AB_CT_DN,$AC_CT_UP,$AC_CT_DN,$AD_CT_UP,$AD_CT_DN,$AB_AD_UP,$AB_AD_DN,$AC_AD_UP,$AC_AD_DN,$AD_AD_UP,$AD_AD_DN,$d";

done < tot_GOs.tmp

rm tot_GOs.tmp
