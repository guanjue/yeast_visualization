data=open('Xu_2009_ORF_Ts_V64_NC_all_X_Y_tRNA.bed','r')
data0=[]
for rec in data:
	tmp=rec.split()
	if not (tmp[3] in data0):
		data0.append(tmp[3])
	else:
		print(tmp[3])

print(len(data0))
