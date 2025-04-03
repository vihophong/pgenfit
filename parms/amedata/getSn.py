from __future__ import print_function
import struct
import numpy as np

file1 = open("rct2_1.mas20.txt")
data1 = []
iline = 0
flagread = 0
flagreadfirst = 0
print("N-Z,N,Z,A,EL,Sn,D_Sn,is_from_experiment,Sp,D_Sp,is_from_experiment,Q4b,D_Q4b,is_from_experiment,Qda,D_Qda,is_from_experiment,Qpa,D_Qpa,is_from_experiment,Qna,D_Qna,is_from_experiment")
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
		"Sn": float(columnitems[3]), "D_Sn": float(columnitems[4]),"is_ex_Sn": is_experiment[3],
		"Sp": float(columnitems[5]), "D_Sp": float(columnitems[6]),"is_ex_Sp": is_experiment[5],
		"Q4b": float(columnitems[7]), "D_Q4b": float(columnitems[8]),"is_ex_Q4b": is_experiment[7],
		"Qda": float(columnitems[9]), "D_Qda": float(columnitems[10]),"is_ex_Qda": is_experiment[9],
		"Qpa": float(columnitems[11]), "D_Qpa": float(columnitems[12]),"is_ex_Qpa": is_experiment[11],
		"Qna": float(columnitems[13]), "D_Qna": float(columnitems[14]),"is_ex_Qna": is_experiment[13]})


for i in range(len(data1)):
	print("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(data1[i]["N-Z"],data1[i]["N"],data1[i]["Z"],data1[i]["A"],data1[i]["EL"],
		data1[i]["Sn"],data1[i]["D_Sn"],data1[i]["is_ex_Sn"],
		data1[i]["Sp"],data1[i]["D_Sp"],data1[i]["is_ex_Sp"],
		data1[i]["Q4b"],data1[i]["D_Q4b"],data1[i]["is_ex_Q4b"],
		data1[i]["Qda"],data1[i]["D_Qda"],data1[i]["is_ex_Qda"],
		data1[i]["Qpa"],data1[i]["D_Qpa"],data1[i]["is_ex_Qpa"],
		data1[i]["Qna"],data1[i]["D_Qna"],data1[i]["is_ex_Qna"]))

np.save("rct2_1.mas20.npy",data1)

