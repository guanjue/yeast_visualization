script_bin='/Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/'

python $script_bin'gff3_expand_sort.py' -f Xu_2009_CUTs_V64.gff3 -u 0 -d 0 -r T -o CUTs.bed
python $script_bin'gff3_expand_sort.py' -f van_Dijk_2011_XUTs_V64.gff3 -u 0 -d 0 -r T -o XUTs.bed
python $script_bin'gff3_expand_sort.py' -f Xu_2009_SUTs_V64.gff3 -u 0 -d 0 -r T -o SUTs.bed
python $script_bin'gff3_expand_sort.py' -f saccharomyces_cerevisiae_Y_prime_element.gff -u 0 -d 0 -r T -o Y_prime_element.bed
python $script_bin'gff3_expand_sort.py' -f saccharomyces_cerevisiae_X_element.gff -u 0 -d 0 -r T -o X_element.bed
python $script_bin'gff3_expand_sort.py' -f saccharomyces_cerevisiae_tRNA_gene_modified.gff -u 0 -d 0 -r T -o tRNA_gene_modified.bed

cat CUTs.bed | awk -F '\t' -v OFS='\t' '{print $4}' > CUTs.genename.txt
cat SUTs.bed | awk -F '\t' -v OFS='\t' '{print $4}' > SUTs.genename.txt
cat XUTs.bed | awk -F '\t' -v OFS='\t' '{print $4}' > XUTs.genename.txt
cat Y_prime_element.bed | awk -F '\t' -v OFS='\t' '{print $4}' > Y_prime_element.genename.txt
cat X_element.bed | awk -F '\t' -v OFS='\t' '{print $4}' > X_element.genename.txt
cat tRNA_gene_modified.bed | awk -F '\t' -v OFS='\t' '{print $4}' > tRNA_gene_modified.genename.txt

tail -n+43 Xu_2009_ORF-Ts_V64.gff3 | awk -F '\t' -v OFS='\t' '{print $9}' | awk -F ';' -v OFS='\t' '{print $1,$2}' | awk -F '=' -v OFS='\t' '{print $2,$3}' | awk -F '\t' -v OFS='\t' '{print $1,$3}' | awk -F '%' -v OFS='\t' '{print $1}'> Xu_2009_ORF-Ts_V64.id2name


python name2id.py -g 'SAGA_TSS_Xu_2009_ORF_Ts_V64.gff' -s 'Xu_2009_ORF-Ts_V64.id2name' -o 'SAGA.genename.txt'
python name2id.py -g 'RP_137_genes_TSS_Xu_2009.gff' -s 'Xu_2009_ORF-Ts_V64.id2name' -o 'RP.genename.txt'
#python name2id.py -g 'TFIID_TSS_Xu_2009_ORF_Ts_V64.gff' -s 'Xu_2009_ORF-Ts_V64.id2name' -o 'TFIID.genename.txt'

python getTFIID.py -a Xu_2009_ORF-Ts_V64.id2name -r RP.genename.txt -s SAGA.genename.txt -t TFIID.genename.txt

mkdir genename
cp *.genename.txt genename/