import csv
import numpy as np

fileNamesMeteoCPC=[ \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp002-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp003-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp004-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp005-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp006-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp007-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp008-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp009-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp010-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp011-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp012-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp013-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp014-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp015-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp016-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp017-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp018-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp019-1-MeteoCPC.csv', \
# 	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp020-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp021-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp022-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp023-1-MeteoCPC.csv', \
# 	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp024-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp025-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp026-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp027-1-MeteoCPC.csv',\
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp028-1-MeteoCPC.csv', \
	'../iSKYLAB-data/Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp029-1-MeteoCPC.csv']

metStr=['MeteoCPC-Exp002','MeteoCPC-Exp003','MeteoCPC-Exp004','MeteoCPC-Exp005',\
	'MeteoCPC-Exp006','MeteoCPC-Exp007','MeteoCPC-Exp008','MeteoCPC-Exp009',\
	'MeteoCPC-Exp010','MeteoCPC-Exp011','MeteoCPC-Exp012','MeteoCPC-Exp013',\
	'MeteoCPC-Exp014','MeteoCPC-Exp015','MeteoCPC-Exp016','MeteoCPC-Exp017',\
	'MeteoCPC-Exp018','MeteoCPC-Exp019',\
	#'MeteoCPC-Exp020', 
	'MeteoCPC-Exp021',\
	'MeteoCPC-Exp022','MeteoCPC-Exp023',\
	#'MeteoCPC-Exp024',
	'MeteoCPC-Exp025',\
	'MeteoCPC-Exp026','MeteoCPC-Exp027','MeteoCPC-Exp028','MeteoCPC-Exp029']

def readData(readThis=3,metStr="MeteoCPC-Exp005"):
	time=[];pres=[];tgwm=[];tgws=[];twwm=[];twws=[];
	tdew=[];tfrost=[];cpci=[];cpctb=[];cpctt=[]
	with open(fileNamesMeteoCPC[readThis],'r') as csvfile: 
		reader=csv.DictReader(csvfile) 
		for row in reader: 
			time.append(float(row['Time (Sec)']))
			pres.append(float(row[' Pressure (hPa)']))
			tgwm.append(float(row[' Tgw_mean (degC)']))
			tgws.append(float(row[' Tgw_std (degC)']))
			twwm.append(float(row[' Tww_mean (degC)']))
			twws.append(float(row[' Tww_std (degC)']))
			tdew.append(float(row[' TDew (degC)']))
			tfrost.append(float(row[' TFrost (degC)']))
			cpci.append(float(row[' CPC_Inters (#/cm3)']))
			cpctb.append(float(row[' CPC_TotBot (#/cm3)']))
			cpctt.append(float(row[' CPC_TotTop (#/cm3)']))
			
	data1 = {metStr : \
		{"Time" : np.array(time), "Pressure" : np.array(pres), \
		"Tgw mean" : np.array(tgwm), "Tgw std": np.array(tgws), \
		"Tww mean" : np.array(twwm), "Tww std" : np.array(twws), \
		"TDew" : np.array(tdew), "TFrost" : np.array(tfrost), \
		"CPC_Inters" : np.array(cpci), "CPC_TotBot" : np.array(cpctb), \
		"CPC_TotTop": np.array(cpctt)}}

	return data1

if __name__== "__main__":
	outputModel=True
	readThis=5
	
	if 'data1' in locals():
		pass
	else:
		data1=dict()
	for i in range(len(metStr)):
		data2=readData(readThis=i,metStr=metStr[i])
		data1[metStr[i]]=data2[metStr[i]].copy()
	
	if outputModel:
		# create data for namelist
		num1=len(data1[metStr[readThis]]['Time'])
		print('n_levels_c = ' + str(num1) + ',')
		str1='	time_chamber(1:' + str(num1) + ')   = '
		for i in range(num1):
			str1 = str1 + str(data1[metStr[readThis]]['Time'][i]) + ','
		str1 = str1 + '\n'
		str1=str1+ '	press_chamber(1:' + str(num1) + ')   = '	
		for i in range(num1):
			str1=str1+str(data1[metStr[readThis]]['Pressure'][i]*100.0) + ','
		str1=str1+ '\n'
		str1=str1+ '	temp_chamber(1:' + str(num1) + ')   = '	
		for i in range(num1):
			str1=str1+str(data1[metStr[readThis]]['Tgw mean'][i]+273.15) + ','
		str1=str1+'\n'
		print(str1)
	
	
	
