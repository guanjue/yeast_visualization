#!/bin/bash

script_bin='/Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/'
while read info; do
	id=$(echo $info | awk -F ' ' -v OFS='\t' '{print $1}' )
	bam_file=$(echo $info | awk -F ' ' -v OFS='\t' '{print $2}' )
	echo $id'_experiment.events'
	echo $bam_file'.bam'

	rm -r $id'_gene_TSS_TES'
	mkdir $id'_gene_TSS_TES'
	cp $id'_experiment.events' $id'_gene_TSS_TES'
	cp $bam_file'.bam' $id'_gene_TSS_TES'
	cd $id'_gene_TSS_TES'
	#####################################
	echo 'get multiGPS pks'
	tail -n+8 $id'_experiment.events' | awk -F '\t' -v OFS='\t' '{print $1}' | awk -F ':' -v OFS='\t' '{print $1,$2,$2+1}' | sort -k1,1 -k2,2n > $id'_experiment.sort.midpoint'
	echo 'expand bed file and sort'
	python $script_bin'gff3_expand_sort_all.py' -f /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/gene_group/gff3_file/Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_uniq.gff3 -u 500_2000 -d 100_2000 -r T -o Xu_2009_ORF-Ts_V64_sort.bed
	sort -u Xu_2009_ORF-Ts_V64_sort.bed  | sort -k1,1 -k2,2n > Xu_2009_ORF-Ts_V64_cosort.bed
	echo 'get bound genes'
	bedtools intersect -a Xu_2009_ORF-Ts_V64_cosort.bed -b $id'_experiment.sort.midpoint' -wa | sort -u > $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed'
	python $script_bin'bed_expand_sort.py' -f $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed' -u -500 -d -100 -r F -o $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed'	
	echo extract bound genes
	cat $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'_Xu_2009_ORF-Ts_V64.genelist.txt'
	echo 'expand bed file and sort'
	python $script_bin'bed_expandmidpoint_sort.py' -f $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' -u 2000 -d 2000 -r F -o $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.bed'
	echo 'gene group split'
	python $script_bin'gene_group_split.py' -t $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.bed' -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F

	#echo 'generate bam.bai'
	samtools index $bam_file'.bam'
	echo 'get heatmap table'
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.bed.genegroup' -s 6 -n 5 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	echo 'generate heatmap'
	python $script_bin'gene_group_split_notation.py' -t $id'_Xu_2009_ORF-Ts_V64_'$bam_file'_read1_combined.tabular' -a 0 -g $script_bin'gene_group_list.txt' -b 0 -i T -n 5
	Rscript $script_bin'heatmap.R' $id'_Xu_2009_ORF-Ts_V64_'$bam_file'_read1_combined.tabular.split_line' $id'_Xu_2009_ORF-Ts_V64_'$bam_file'_read1_combined.tabular.colorsplit' $id'_Xu_2009_ORF-Ts_V64_'$bam_file'_read1_combined' darkorange2 0.97
	cp $id'_Xu_2009_ORF-Ts_V64_'$bam_file'_read1_combined.png' $id'_image_10.png'
	cp $id'_Xu_2009_ORF-Ts_V64_'$bam_file'_read1_combined.color_split.png' $id'_image_11.png'
	cp $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.bed.genegroup' $id'_genemidpoint.bed'

	### for all genes
	python $script_bin'gff3_expand_sort_all.py' -f /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/gene_group/gff3_file/Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_uniq.gff3 -u 0_2000 -d 0_2000 -r T -o Xu_2009_ORF-Ts_V64_sort_OD.bed
	sort -u Xu_2009_ORF-Ts_V64_sort_OD.bed > Xu_2009_ORF-Ts_V64_sort_uniq.bed
	python $script_bin'bed_expandmidpoint_sort.py' -f Xu_2009_ORF-Ts_V64_sort_uniq.bed -u 2000 -d 2000 -r F -o Xu_2009_ORF-Ts_V64_sort_midpoint_expand.bed
	echo 'gene group split'
	python $script_bin'gene_group_split.py' -t Xu_2009_ORF-Ts_V64_sort_midpoint_expand.bed -a 3 -g $script_bin'gene_group_list_9.txt' -b 0 -i F
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c Xu_2009_ORF-Ts_V64_sort_midpoint_expand.bed.genegroup -s 6 -n 1 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	python $script_bin'bin_row_col.py' -i 'Xu_2009_ORF-Ts_V64_sort_midpoint_expand_'$bam_file'_read1_combined.tabular' -r 25 -c 5 -o 'Xu_2009_ORF-Ts_V64_sort_midpoint_expand_'$bam_file'_read1_combined_binned.tabular'
	python $script_bin'gene_group_split_notation.py' -t 'Xu_2009_ORF-Ts_V64_sort_midpoint_expand_'$bam_file'_read1_combined_binned.tabular' -a 0 -g $script_bin'gene_group_list_9.txt' -b 0 -i T -n 5
	Rscript $script_bin'heatmap.R' 'Xu_2009_ORF-Ts_V64_sort_midpoint_expand_'$bam_file'_read1_combined_binned.tabular.split_line' 'Xu_2009_ORF-Ts_V64_sort_midpoint_expand_'$bam_file'_read1_combined_binned.tabular.colorsplit' 'Xu_2009_ORF-Ts_V64_midpointexpand_'$bam_file'_read1_combined' black 0.97
	cp 'Xu_2009_ORF-Ts_V64_midpointexpand_'$bam_file'_read1_combined.png' $id'_image_18.png'
	cp 'Xu_2009_ORF-Ts_V64_midpointexpand_'$bam_file'_read1_combined.color_split.png' $id'_image_19.png'
	cp Xu_2009_ORF-Ts_V64_sort_midpoint_expand.bed.genegroup $id'_genemidpoint_all.bed'

	#java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$id'_Xu_2009_ORF-Ts_V64_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $id'_Xu_2009_ORF-Ts_V64.gene_len.png' -a 0 -c 20 -f png

	#####################################
	echo 'bound genes TSS'
	cat $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > $id'_Xu_2009_ORF-Ts_V64_TSS.bed'
	echo 'sort by dist to TSS'
	#cat $id'_experiment.sort.midpoint' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,".",".","+"}' > $id'_experiment.midpoint.pn.txt'
	#cat $id'_experiment.sort.midpoint' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,".",".","-"}' >> $id'_experiment.midpoint.pn.txt'
	#sort -k1,1 -k2,2n $id'_experiment.midpoint.pn.txt' > $id'_experiment.sort.midpoint.pn.txt'
	bedtools closest -a $id'_Xu_2009_ORF-Ts_V64_TSS.bed' -b $id'_experiment.sort.midpoint' -D a -t first > $id'_Xu_2009_ORF-Ts_V64_TSS.bound.txt'
	sort -k10,10nr $id'_Xu_2009_ORF-Ts_V64_TSS.bound.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-750,$2+251,$4,$5,$6; else print $1,$2-250,$2+751,$4,$5,$6}' > $id'_Xu_2009_ORF-Ts_V64_TSS.bound.sort.bed'
	echo 'gene group split'
	python $script_bin'gene_group_split.py' -t $id'_Xu_2009_ORF-Ts_V64_TSS.bound.sort.bed' -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F
	echo 'get heatmap table'
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $id'_Xu_2009_ORF-Ts_V64_TSS.bound.sort.bed.genegroup' -s 6 -n 1 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	echo 'generate heatmap'
	python $script_bin'gene_group_split_notation.py' -t $id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular' -a 0 -g $script_bin'gene_group_list.txt' -b 0 -i T -n 5
	Rscript $script_bin'heatmap.R' $id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular.split_line' $id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular.colorsplit' $id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined' aquamarine4 0.97
	#java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $id'_Xu_2009_ORF-Ts_V64.TSS_ref.png' -a 0 -c 20 -f png
	cp $id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.png' $id'_image_13.png'
	cp $id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.color_split.png' $id'_image_13a.png'
	cp $id'_Xu_2009_ORF-Ts_V64_TSS.bound.sort.bed.genegroup' $id'_TSS.bed'

	#####################################
	echo 'bound genes TES'
	cat $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$3-1,$3,$4,$5,$6; else print $1,$2,$2+1,$4,$5,$6}' | sort -k1,1 -k2,2n > $id'_Xu_2009_ORF-Ts_V64_TES.bed'
	echo 'sort by dist to TES'
	#cat $id'_experiment.sort.midpoint' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,".",".","+"}' > $id'_experiment.midpoint.pn.txt'
	#cat $id'_experiment.sort.midpoint' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,".",".","-"}' >> $id'_experiment.midpoint.pn.txt'
	#sort -k1,1 -k2,2n $id'_experiment.sort.midpoint.pn.txt' > $id'_experiment.sort.midpoint.pn.txt'
	bedtools closest -a $id'_Xu_2009_ORF-Ts_V64_TES.bed' -b $id'_experiment.sort.midpoint' -D a -t first > $id'_Xu_2009_ORF-Ts_V64_TES.bound.txt'
	sort -k10,10n $id'_Xu_2009_ORF-Ts_V64_TES.bound.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+501,$4,$5,$6; else print $1,$2-500,$2+501,$4,$5,$6}' > $id'_Xu_2009_ORF-Ts_V64_TES.bound.sort.bed'
	echo 'gene group split'
	python $script_bin'gene_group_split.py' -t $id'_Xu_2009_ORF-Ts_V64_TES.bound.sort.bed' -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F
	echo 'get heatmap table'
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $id'_Xu_2009_ORF-Ts_V64_TES.bound.sort.bed.genegroup' -s 6 -n 1 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	echo 'generate heatmap'
	python $script_bin'gene_group_split_notation.py' -t $id'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.tabular' -a 0 -g $script_bin'gene_group_list.txt' -b 0 -i T -n 5
	Rscript $script_bin'heatmap.R' $id'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.tabular.split_line' $id'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.tabular.colorsplit' $id'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined' brown3 0.97
	#java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$id'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $id'_Xu_2009_ORF-Ts_V64.TSS_ref.png' -a 0 -c 20 -f png
	cp $id'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.png' $id'_image_15.png'
	cp $id'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.color_split.png' $id'_image_15a.png'
	cp $id'_Xu_2009_ORF-Ts_V64_TES.bound.sort.bed.genegroup' $id'_TES.bed'

	#####################################
	##################
	### nucleosome part
	echo 'bound nfr'
	#python $script_bin'extract_bound.py' -b $id'_Xu_2009_ORF-Ts_V64.genelist.txt' -t /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/Yeast_plus_one_sacCer3.gff -i 4 -o $id'_Yeast_plus_one_sacCer3.bed'
	#echo 'sort by dist to TSS'
	#bedtools closest -a $id'_Yeast_plus_one_sacCer3.bed' -b $id'_experiment.sort.midpoint' -D a > $id'_Yeast_plus_one_sacCer3.bound.txt'
	#sort -k10,10nr $id'_Yeast_plus_one_sacCer3.bound.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+501,$4,$5,$6; else print $1,$2-500,$2+501,$4,$5,$6}' > $id'_Yeast_plus_one_sacCer3.bound.sort.bed'
	#echo 'gene group split'
	#python $script_bin'gene_group_split.py' -t $id'_Yeast_plus_one_sacCer3.bound.sort.bed' -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F
	#echo 'get heatmap table'
	python $script_bin'extract_bound.py' -b $id'_Xu_2009_ORF-Ts_V64.genelist.txt' -t /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/gene_group/nuclesome/Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort.bed -i 4 -o $id'_Yeast_nfrmidpoint_sacCer3.bed'
	python $script_bin'gene_group_split.py' -t $id'_Yeast_nfrmidpoint_sacCer3.bed' -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F

	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $id'_Yeast_nfrmidpoint_sacCer3.bed.genegroup' -s 6 -n 1 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	echo 'generate heatmap'
	python $script_bin'gene_group_split_notation.py' -t $id'_Yeast_nfrmidpoint_sacCer3_'$bam_file'_read1_combined.tabular' -a 0 -g $script_bin'gene_group_list.txt' -b 0 -i T -n 5
	Rscript $script_bin'heatmap.R' $id'_Yeast_nfrmidpoint_sacCer3_'$bam_file'_read1_combined.tabular.split_line' $id'_Yeast_nfrmidpoint_sacCer3_'$bam_file'_read1_combined.tabular.colorsplit' $id'_Yeast_nfrmidpoint_sacCer3_'$bam_file'_read1_combined' cornflowerblue 0.97
	cp $id'_Yeast_nfrmidpoint_sacCer3_'$bam_file'_read1_combined.png' $id'_image_17.png'
	cp $id'_Yeast_nfrmidpoint_sacCer3_'$bam_file'_read1_combined.color_split.png' $id'_image_17a.png'
	cp $id'_Yeast_nfrmidpoint_sacCer3.bed.genegroup' $id'_nfrmidpoint.bed'

	#java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $id'_Xu_2009_ORF-Ts_V64.TSS_ref.png' -a 0 -c 20 -f png



	python $script_bin'gene_group_split.py' -t /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/gene_group/nuclesome/Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort.bed -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/gene_group/nuclesome/Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort.bed.genegroup -s 6 -n 1 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	python $script_bin'bin_row_col.py' -i 'Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort_'$bam_file'_read1_combined.tabular' -r 5 -c 5 -o 'Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort_'$bam_file'_read1_combined_binned.tabular'
	python $script_bin'gene_group_split_notation.py' -t 'Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort_'$bam_file'_read1_combined_binned.tabular' -a 0 -g $script_bin'gene_group_list.txt' -b 0 -i T -n 5
	Rscript $script_bin'heatmap.R' 'Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort_'$bam_file'_read1_combined_binned.tabular.split_line' 'Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort_'$bam_file'_read1_combined_binned.tabular.colorsplit' 'Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA_NFRcenter_1001bp_nfr_sort_'$bam_file'_read1_combined' black 0.97



	cd ..
done < id_list.txt

