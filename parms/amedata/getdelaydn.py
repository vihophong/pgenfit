from __future__ import print_function
import struct
import numpy as np

file1 = open("delayedn_eval.csv")
data1 = []

iline = 0
print("N-Z,N,Z,A,EL,liso,Ex,D_Ex,Bm,D_Bm,Qb1n,D_Qb1n,Qb2n,D_Qb2n,Qb3n,D_Qb3n,T12,D_T12,P1n,D_P1n,P2n,D_P2n,P3n,D_P3n")
while True:
	iline += 1
	line = file1.readline()
	if not line:
		break
	if (iline<3):
		continue
	items = line.split(",")
	for i in range(0,len(items)):
		items[i]=items[i].strip()
		#print(items[i])
	#print("----")
	ele = ""
	for i in range(len(items[0])):
		if (items[0][i].isalpha()):
			ele += items[0][i]
	ele = ele.title()
	data1.append({"N-Z": int(items[2])-int(items[1])-int(items[1]),"N": int(items[2])-int(items[1]),"Z": int(items[1]),"A": int(items[2]),"EL": ele,"liso": int(items[3]),
		 "Ex": float(items[4]),"D_Ex": float(items[5]),
		 "Bm": float(items[6]),"D_Bm": float(items[7]),
		 "Qb1n": float(items[8]),"D_Qb1n": float(items[9]),
		 "Qb2n": float(items[10]),"D_Qb2n": float(items[11]),
		 "Qb3n": float(items[12]),"D_Qb3n": float(items[13]),
		 "T12": float(items[14]),"D_T12": float(items[15]),
		 "P1n": float(items[16]),"D_P1n": float(items[17]),
		 "P2n": float(items[18]),"D_P2n": float(items[19]),
		 "P3n": float(items[20]),"D_P3n": float(items[21])})

for i in range(len(data1)):
	print("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(data1[i]["N-Z"],data1[i]["N"],data1[i]["Z"],data1[i]["A"],data1[i]["EL"],data1[i]["liso"],
		data1[i]["Ex"],data1[i]["D_Ex"],
		data1[i]["Bm"],data1[i]["D_Bm"],
		data1[i]["Qb1n"],data1[i]["D_Qb1n"],
		data1[i]["Qb2n"],data1[i]["D_Qb2n"],
		data1[i]["Qb3n"],data1[i]["D_Qb3n"],
		data1[i]["T12"],data1[i]["D_T12"],
		data1[i]["P1n"],data1[i]["D_P1n"],
		data1[i]["P2n"],data1[i]["D_P2n"],
		data1[i]["P3n"],data1[i]["D_P3n"]))
np.save("delayedn_eval.npy",data1)

