from __future__ import print_function
import struct
import numpy as np

file1 = open("rct1.mas20.txt")
data1 = []
iline = 0
flagread = 0
flagreadfirst = 0
print("N-Z,N,Z,A,EL,S2n,D_S2n,is_from_experiment,S2p,D_S2p,is_from_experiment,Qa,D_Qa,is_from_experiment,Q2b,D_Q2b,is_from_experiment,Qep,D_Qep,is_from_experiment,Qbn,D_Qbn,is_from_experiment")
while True:
	iline += 1
	line = file1.readline()
	if not line:
		break
	items = line.split()
	if (len(items)==0):
		continue
	if (line[0]=="1"):
		flagread += 1
	if (flagread<2):
		continue
	if (line[0]=="0"):
		flagreadfirst += 1
	if (flagreadfirst==0):
		continue
	if line[0]=="0":
		items = line[1:len(line)].split()

	columnitems = []
	for i in range(0,len(items)):
		columnitems.append(items[i])
		if items[i]=="*":
			columnitems.append("999999")
		#print items[i],
	#print("")
	is_experiment = [1] * 15
	for i in range(len(columnitems)):
		if (columnitems[i]=="*"):
			columnitems[i] = "999999"
		if (columnitems[i][len(columnitems[i])-1]=="#"):
			columnitems[i] = columnitems[i][:len(columnitems[i])-1]
			is_experiment[i]=0
		
	if (len(columnitems)!=15):
		print("Error")
	data1.append({"N-Z": int(columnitems[0])-int(columnitems[2])-int(columnitems[2]),"N": int(columnitems[0])-int(columnitems[2]),"Z": int(columnitems[2]),"A": int(columnitems[0]),"EL": str(columnitems[1]),
		"S2n": float(columnitems[3]), "D_S2n": float(columnitems[4]),"is_ex_S2n": is_experiment[3],
		"S2p": float(columnitems[5]), "D_S2p": float(columnitems[6]),"is_ex_S2p": is_experiment[5],
		"Qa": float(columnitems[7]), "D_Qa": float(columnitems[8]),"is_ex_Qa": is_experiment[7],
		"Q2b": float(columnitems[9]), "D_Q2b": float(columnitems[10]),"is_ex_Q2b": is_experiment[9],
		"Qep": float(columnitems[11]), "D_Qep": float(columnitems[12]),"is_ex_Qep": is_experiment[11],
		"Qbn": float(columnitems[13]), "D_Qbn": float(columnitems[14]),"is_ex_Qbn": is_experiment[13]})


for i in range(len(data1)):
	print("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(data1[i]["N-Z"],data1[i]["N"],data1[i]["Z"],data1[i]["A"],data1[i]["EL"],
		data1[i]["S2n"],data1[i]["D_S2n"],data1[i]["is_ex_S2n"],
		data1[i]["S2p"],data1[i]["D_S2p"],data1[i]["is_ex_S2p"],
		data1[i]["Qa"],data1[i]["D_Qa"],data1[i]["is_ex_Qa"],
		data1[i]["Q2b"],data1[i]["D_Q2b"],data1[i]["is_ex_Q2b"],
		data1[i]["Qep"],data1[i]["D_Qep"],data1[i]["is_ex_Qep"],
		data1[i]["Qbn"],data1[i]["D_Qbn"],data1[i]["is_ex_Qbn"]))

np.save("rct1.mas20.npy",data1)

