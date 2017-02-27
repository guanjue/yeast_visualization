def name2id(gff_file,id2namefile,outputname):
	data0=open(gff_file,'r')
	data01=[]

	for records in data0:
		data01.append(records.split()[8])

	data1=open(id2namefile,'r')
	data11={}

	for records in data1:
		tmp=records.split()
		data11[tmp[1]]=tmp[0]
	#print(data11)
	result=open(outputname,'w')
	for records in data01:
		if records in data11:
			result.write(data11[records]+'\n')
	result.close()
	data1.close()
	data0.close()


############################################################################
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hg:s:o:")
	except getopt.GetoptError:
		print 'python name2id.py -g <gff_file> -s <id2namefile> -o <outputname>'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python name2id.py -g <gff_file> -s <id2namefile> -o <outputname>'
			sys.exit()
		elif opt=="-g":
			gff_file=str(arg.strip())
		elif opt=="-s":
			id2namefile=str(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())
	name2id(gff_file,id2namefile,outputname)

if __name__=="__main__":
	main(sys.argv[1:])