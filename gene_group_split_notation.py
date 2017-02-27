def gene_group_split(target_filename,target_id_column,gene_group_list,gene_group_id_column,header,split_line_num):
	gene_group0=[]
	#gene_group_id={1:'RP',2:'SAGA',3:'TFIID',4:'CUTs',5:'SUTs',6:'XUTs',7:'tRNA',8:'X_elements'}
	gene_group_file0=open(gene_group_list,'r')
	for records in gene_group_file0:
		gene_group0.append(records.split()[0])
	#gene_group0=['RP_137_genes_Xu_2009','SAGA_Xu_2009_ORF_Ts_V64','TFIID_Xu_2009_ORF_Ts_V64']
	target=open(target_filename,'r')
	target0=[]
	for records in target:
		target0.append(records.split())

	target0_genegroup = open(target_filename+'.split_line','w')
	target0_genegroup_colorsplit = open(target_filename+'.colorsplit','w')

	if header=='T':
		## write header
		for rec in target0[0]:
			target0_genegroup.write(rec+'\t')
		target0_genegroup.write('\n')
	j=0
	target0_reduce=[]
	for gene_group in gene_group0:
		j=j+1
		for h in range(0,split_line_num):
			target0_genegroup.write('split_line'+'\t''split_line'+'\t')
			for i in range(2,len(target0[0])):
				target0_genegroup.write('0'+'\t')
			target0_genegroup.write('\n')
			target0_genegroup_colorsplit.write(str(j)+'\n')
		### read gene group gene names
		gene_group_file = open(gene_group,'r')
		gene_group1={}
		for records in gene_group_file:
			gene_group1[records.split()[gene_group_id_column]]=' '
		gene_group_file.close()

		## write data
		if header=='T':
			start_row=1
		else:
			start_row=0
		k=0
		for records in target0[0:]:

			if records[target_id_column] in gene_group1:
				k=k+1
				for rec in records:
					target0_genegroup.write(rec+'\t')
				target0_genegroup.write('\n')
				target0_genegroup_colorsplit.write(str(j)+'\n')
			else:
				target0_reduce.append(records)
		target0=target0_reduce
		target0_reduce=[]
		print(gene_group)
		print(k)
	target0_genegroup.close()
	target0_genegroup_colorsplit.close()


############################################################################
### python gene_group_split.py -t Xu_2009_ORF-Ts_V64_sort_pk_plus1.bed -a 3 -g gene_group_list.txt -b 0 -i F
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h:t:a:g:b:i:n:")
	except getopt.GetoptError:
		print 'python gene_group_split.py -t <target_filename> -a <target_id_column> -g <gene_group_list> -b <gene_group_id_column> -i <header T or F> -n <split line number>'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python gene_group_split.py -t <target_filename> -g <gene_group_list> -g <gene_group_list> -b <gene_group_id_column> -i <header T or F> -n <split line number>'
			sys.exit()
		elif opt=="-t":
			target_filename=str(arg.strip())
		elif opt=="-a":
			target_id_column=int(arg.strip())
		elif opt=="-g":
			gene_group_list=str(arg.strip())
		elif opt=="-b":
			gene_group_id_column=int(arg.strip())
		elif opt=="-i":
			header=str(arg.strip())
		elif opt=="-n":
			split_line_num=int(arg.strip())

	gene_group_split(target_filename,target_id_column,gene_group_list,gene_group_id_column,header,split_line_num)
if __name__=="__main__":
	main(sys.argv[1:])
