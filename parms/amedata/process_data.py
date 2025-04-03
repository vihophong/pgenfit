from __future__ import print_function
import struct
import numpy as np
import math

datain = np.load("delayedn_eval.npy",allow_pickle='TRUE')

data_masstable = np.load("mass_1.mas20.npy",allow_pickle='TRUE')
data_Qbn = np.load("rct1.mas20.npy",allow_pickle='TRUE')
data_Sn = np.load("rct2_1.mas20.npy",allow_pickle='TRUE')

data1 = []

for i in range(len(datain)):
	if (datain[i]["liso"]!=0):
		data1.append({"N-Z": datain[i]["N-Z"],"N": datain[i]["N"],"Z": datain[i]["Z"],"A": datain[i]["A"],"EL": datain[i]["EL"],"liso": datain[i]["liso"],
		 "Ex": datain[i]["Ex"],"D_Ex": datain[i]["D_Ex"],
		 "Bm": datain[i]["Bm"],"D_Bm": datain[i]["D_Bm"],
		 "Qb1n": datain[i]["Qb1n"] ,"D_Qb1n": datain[i]["D_Qb1n"],
		 "Qb2n": datain[i]["Qb2n"],"D_Qb2n": datain[i]["D_Qb2n"],
		 "Qb3n": datain[i]["Qb3n"],"D_Qb3n": datain[i]["D_Qb3n"],
		 "T12": datain[i]["T12"],"D_T12": datain[i]["D_T12"],
		 "P1n": datain[i]["P1n"],"D_P1n": datain[i]["D_P1n"],
		 "P2n": datain[i]["P2n"],"D_P2n": datain[i]["D_P2n"],
		 "P3n": datain[i]["P3n"],"D_P3n": datain[i]["D_P3n"],
		 "new_Qb": 0,"new_D_Qb": 0,
		 "new_Qb1n": 0,"new_D_Qb1n": 0,
		 "new_Qb2n": 0,"new_D_Qb2n": 0,
		 "new_S1n": 0 ,"new_D_S1n": 0,
		 "new_S2n": 0,"new_D_S2n": 0		 
		 })
	else:
		Qbn_entry = next((item for item in data_Qbn if (item["A"] == datain[i]["A"] and item["Z"] == datain[i]["Z"] and item["EL"] == datain[i]["EL"])), None)
		if (Qbn_entry==None):
			print("Error")
		Qb_entry = next((item for item in data_masstable if (item["A"] == datain[i]["A"] and item["Z"] == datain[i]["Z"] and item["EL"] == datain[i]["EL"])), None)
		if (Qb_entry==None):
			print("Error")
		Sn_daugter_entry = next((item for item in data_Sn if (item["A"] == datain[i]["A"] and item["Z"] == datain[i]["Z"]+1)), None)
		if (Sn_daugter_entry==None):
			print("Error")
		S2n_daugter_entry = next((item for item in data_Qbn if (item["A"] == datain[i]["A"] and item["Z"] == datain[i]["Z"]+1)), None)
		if (S2n_daugter_entry==None):
			print("Error")
		#print(Qbn_entry["Qbn"]-datain[i]["Qb1n"])
		#print(Qbn_entry["D_Qbn"]-datain[i]["D_Qb1n"])
		Qbn = Qb_entry["Qb"]-Sn_daugter_entry["Sn"]
		D_Qbn = math.sqrt(Qb_entry["DQb"]*Qb_entry["DQb"]+Sn_daugter_entry["D_Sn"]*Sn_daugter_entry["D_Sn"])
		Qb2n = Qb_entry["Qb"]-S2n_daugter_entry["S2n"]
		D_Qb2n = math.sqrt(Qb_entry["DQb"]*Qb_entry["DQb"]+S2n_daugter_entry["D_S2n"]*S2n_daugter_entry["D_S2n"])
		#print(Qbn_entry["Qbn"]-Qbn)
		# if (datain[i]["Qb2n"]>0):
		# 	print(Qb2n-datain[i]["Qb2n"])
		# 	print(D_Qb2n-datain[i]["D_Qb2n"])
		# if (datain[i]["Qb2n"]==0):
		# 	print(Qb2n)
		# if (datain[i]["Qb1n"]==0):
		# 	print(Qbn)
		# if (Qbn>1):
		# 	print(datain[i]["A"],datain[i]["EL"],Qbn,datain[i]["Qb1n"])
		# 	print(datain[i]["A"],datain[i]["EL"],Qbn-Qbn_entry["Qbn"])
		# 	print(datain[i]["A"],datain[i]["EL"],D_Qbn-Qbn_entry["D_Qbn"])
		# if (Qb2n>1):
		# 	print(datain[i]["A"],datain[i]["EL"],Qb2n,datain[i]["Qb2n"])
		# 	print(datain[i]["A"],datain[i]["EL"],D_Qb2n,datain[i]["D_Qb2n"])
		if (Qbn<1):
			Qbn = 0
			D_Qbn = 0
		if (Qb2n<1):
			Qb2n = 0
			D_Qb2n = 0
		if (Qb_entry["Qb"]==999999):
			Qb_entry["Qb"]=0
			Qb_entry["DQb"]=0
			Qbn = 0
			D_Qbn = 0
			Qb2n = 0
			D_Qb2n = 0
		data1.append({"N-Z": datain[i]["N-Z"],"N": datain[i]["N"],"Z": datain[i]["Z"],"A": datain[i]["A"],"EL": datain[i]["EL"],"liso": datain[i]["liso"],
		 "Ex": datain[i]["Ex"],"D_Ex": datain[i]["D_Ex"],
		 "Bm": datain[i]["Bm"],"D_Bm": datain[i]["D_Bm"],
		 "Qb1n": datain[i]["Qb1n"],"D_Qb1n": datain[i]["D_Qb1n"],
		 "Qb2n": datain[i]["Qb2n"],"D_Qb2n": datain[i]["D_Qb2n"],
		 "Qb3n": datain[i]["Qb3n"],"D_Qb3n": datain[i]["D_Qb3n"],
		 "T12": datain[i]["T12"],"D_T12": datain[i]["D_T12"],
		 "P1n": datain[i]["P1n"],"D_P1n": datain[i]["D_P1n"],
		 "P2n": datain[i]["P2n"],"D_P2n": datain[i]["D_P2n"],
		 "P3n": datain[i]["P3n"],"D_P3n": datain[i]["D_P3n"],
		 "new_Qb": Qb_entry["Qb"],"new_D_Qb": Qb_entry["DQb"],
		 "new_Qb1n": Qbn_entry["Qbn"],"new_D_Qb1n": Qbn_entry["D_Qbn"],
		 "new_Qb2n": Qb2n,"new_D_Qb2n": D_Qb2n,
		 "new_S1n":Sn_daugter_entry["Sn"],"new_D_S1n":Sn_daugter_entry["D_Sn"],
		 "new_S2n":S2n_daugter_entry["S2n"],"new_D_S2n":S2n_daugter_entry["D_S2n"]
		 })
	

print("N-Z,N,Z,A,EL,liso,Ex,D_Ex,Bm,D_Bm,Qb1n,D_Qb1n,Qb2n,D_Qb2n,Qb3n,D_Qb3n,T12,D_T12,P1n,D_P1n,P2n,D_P2n,P3n,D_P3n,AME2021_Qb,AME2021_D_Qb,AME2021_Qb1n,AME2021_D_Qb1n,AME2021_Qb2n,AME2021_D_Qb2n,AME2021_S1n,AME2021_D_S1n,AME2021_S2n,AME2021_D_S2n")
for i in range(len(data1)):
	print("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(data1[i]["N-Z"],data1[i]["N"],data1[i]["Z"],data1[i]["A"],data1[i]["EL"],data1[i]["liso"],
		data1[i]["Ex"],data1[i]["D_Ex"],
		data1[i]["Bm"],data1[i]["D_Bm"],
		data1[i]["Qb1n"],data1[i]["D_Qb1n"],
		data1[i]["Qb2n"],data1[i]["D_Qb2n"],
		data1[i]["Qb3n"],data1[i]["D_Qb3n"],
		data1[i]["T12"],data1[i]["D_T12"],
		data1[i]["P1n"],data1[i]["D_P1n"],
		data1[i]["P2n"],data1[i]["D_P2n"],
		data1[i]["P3n"],data1[i]["D_P3n"],
		data1[i]["new_Qb"],data1[i]["new_D_Qb"],
		data1[i]["new_Qb1n"],data1[i]["new_D_Qb1n"],
		data1[i]["new_Qb2n"],data1[i]["new_D_Qb2n"],
		data1[i]["new_S1n"],data1[i]["new_D_S1n"],
		data1[i]["new_S2n"],data1[i]["new_D_S2n"]
		))
np.save("delayedn_eval_with_ame20.npy",data1)

