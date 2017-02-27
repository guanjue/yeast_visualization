script_bin='/Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/'

python $script_bin'gff3_expand_sort_all.py' -f /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/gene_group/gff3_file/Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_uniq.gff3 -u 0_2000 -d 0_2000 -r T -o Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA.bed
cat Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA.bed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort.bed

bedtools closest -a Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort.bed -b h2a_h3_merged.midpoint.bed -iu -D a -t first -io > Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort_plus1.1bp.bed
bedtools closest -a Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort.bed -b h2a_h3_merged.midpoint.bed -id -D a -t last -io > Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort_minus1.1bp.bed

paste Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort_minus1.1bp.bed Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort_plus1.1bp.bed | awk -F '\t' -v OFS='\t' '{print $1,$2-500,$3+500,$4,$5,$6,$7,$8,$9,$10,$17,$18,$19,$20}' > Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort_1000bp_pm.txt
paste Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort_minus1.1bp.bed Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort_plus1.1bp.bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$17,$18,$19,$20,$20-$10}' | sort -u | sort -k15,15n > Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSScenter_1bp_nfr_sort.txt
paste Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort_minus1.1bp.bed Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_TSS_sort_plus1.1bp.bed | awk -F '\t' -v OFS='\t' '{print $1,int(($8+$18)/2)-500,int(($8+$18)/2)+501,$4,$5,$6,$7,$8,$9,$10,$17,$18,$19,$20,$20-$10}' | sort -u | sort -k15,15n > Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort.txt

cat Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}'> Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort.bed









paste Xu_2009_ORF-Ts_V64_minus1.1bp.bed Xu_2009_ORF-Ts_V64_plus1.1bp.bed | awk -F '\t' -v OFS='\t' '{print $1,$2-500,$3+500,$4,$5,$6,$7,$8,$9,$10,$17,$18,$19,$20,$20-$10}' | sort -k14,14n > Xu_2009_ORF-Ts_V64_1000bp_pm_sort_TSSplus1.bed
paste Xu_2009_ORF-Ts_V64_minus1.1bp.bed Xu_2009_ORF-Ts_V64_plus1.1bp.bed | awk -F '\t' -v OFS='\t' '{print $1,$18-500,$19+500,$4,$20,$6}' | sort -k5,5n > Xu_2009_ORF-Ts_V64_1000bp_pm_sort_plus1_distance2TSS.bed

mv Xu_2009_ORF-Ts_V64_plus1.1bp.bed Xu_2009_ORF-Ts_V64_plus1.1bp.3intersect.bed
mv Xu_2009_ORF-Ts_V64_minus1.1bp.bed Xu_2009_ORF-Ts_V64_minus1.1bp.3intersect.bed

cat Xu_2009_ORF-Ts_V64_plus1.1bp.3intersect.bed | awk -F '\t' -v OFS='\t' '{print $7,$8-500,$9+500,$4,$10,$6}' | sort -k5,5n > Xu_2009_ORF-Ts_V64_plus1.1bp.3intersect_sort.bed
cat Xu_2009_ORF-Ts_V64_minus1.1bp.3intersect.bed | awk -F '\t' -v OFS='\t' '{print $7,$8-500,$9+500,$4,$10,$6}' | sort -k5,5n > Xu_2009_ORF-Ts_V64_minus1.1bp.3intersect_sort.bed

sort -k4,4 Xu_2009_ORF-Ts_V64_plus1.1bp.3intersect_sort.bed > Xu_2009_ORF-Ts_V64_plus1.1bp.3intersect_name_sort.bed
sort -k4,4 Xu_2009_ORF-Ts_V64_minus1.1bp.3intersect_sort.bed > Xu_2009_ORF-Ts_V64_minus1.1bp.3intersect_name_sort.bed

paste Xu_2009_ORF-Ts_V64_minus1.1bp.3intersect_name_sort.bed Xu_2009_ORF-Ts_V64_plus1.1bp.3intersect_name_sort.bed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $0,$8-$2; else print $0,$2-$8}' | sort -k13,13n | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > Xu_2009_ORF-Ts_V64_minus1plus1.1bp.3intersect_nfr_sort.txt
paste Xu_2009_ORF-Ts_V64_minus1.1bp.3intersect_name_sort.bed Xu_2009_ORF-Ts_V64_plus1.1bp.3intersect_name_sort.bed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $0,$8-$2; else print $0,$2-$8}' | sort -k13,13n | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > Xu_2009_ORF-Ts_V64_minus1plus1.1bp.3intersect_nfr_sort.tmp.txt

cat Xu_2009_ORF-Ts_V64_minus1plus1.1bp.3intersect_nfr_sort.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' > Xu_2009_ORF-Ts_V64_minus1.1bp.3intersect_nfr_sort.bed
cat Xu_2009_ORF-Ts_V64_minus1plus1.1bp.3intersect_nfr_sort.txt | awk -F '\t' -v OFS='\t' '{print $7,$8,$9,$10,$11,$12}' > Xu_2009_ORF-Ts_V64_plus1.1bp.3intersect_nfr_sort.bed







