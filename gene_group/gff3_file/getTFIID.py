def getTFIID(all_id,rp_id,saga_id,tfiid_id):
	data0=open(all_id,'r')
	data01=[]
	for records in data0:
		data01.append(records.split()[0])

	data1=open(rp_id,'r')
	data11={}

	for records in data1:
		tmp=records.split()
		data11[tmp[0]]=' '

	data2=open(saga_id,'r')

	for records in data2:
		tmp=records.split()
		data11[tmp[0]]=' '

	result=open(tfiid_id,'w')
	for records in data01:
		if not(records in data11):
			result.write(records+'\n')
	result.close()
	data2.close()
	data1.close()
	data0.close()


############################################################################
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"ha:r:s:t:")
	except getopt.GetoptError:
		print 'python getTFIID.py -a <all_ids> -r <rp_id> -s <saga_id> -t <tfiid_id>'
		sys.exit(2)
	for opt,arg in opts:
		if opt=="-h":
			print 'python getTFIID.py -a <all_ids> -r <rp_id> -s <saga_id> -t <tfiid_id>'
			sys.exit()
		elif opt=="-a":
			all_id=str(arg.strip())
		elif opt=="-r":
			rp_id=str(arg.strip())
		elif opt=="-s":
			saga_id=str(arg.strip())
		elif opt=="-t":
			tfiid_id=str(arg.strip())
	getTFIID(all_id,rp_id,saga_id,tfiid_id)

if __name__=="__main__":
	main(sys.argv[1:])