import csv
import numpy as np

fileNamesMeteoCPC=[ \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp002-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp003-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp004-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp005-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp006-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp007-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp008-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp009-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp010-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp011-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp012-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp013-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp014-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp015-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp016-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp017-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp018-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp019-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp020-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp021-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp022-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp023-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp024-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp025-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp026-1-MeteoCPC.csv', \
	'../../../iSKYLAB01/Datasets/iSKYLAB01-Exp027-1-MeteoCPC.csv']

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
	readThis=4
	metStr='MeteoCPC-Exp006'
	
	data1=readData(readThis=readThis,metStr=metStr)
	
	if outputModel:
		# create data for namelist
		num1=len(data1[metStr]['Time'])
		print('n_levels_c = ' + str(num1) + ',')
		str1='	time_chamber(1:' + str(num1) + ')   = '
		for i in range(num1):
			str1 = str1 + str(data1[metStr]['Time'][i]) + ','
		str1 = str1 + '\n'
		str1=str1+ '	press_chamber(1:' + str(num1) + ')   = '	
		for i in range(num1):
			str1=str1+str(data1[metStr]['Pressure'][i]*100.0) + ','
		str1=str1+ '\n'
		str1=str1+ '	temp_chamber(1:' + str(num1) + ')   = '	
		for i in range(num1):
			str1=str1+str(data1[metStr]['Tgw mean'][i]+273.15) + ','
		str1=str1+'\n'
		print(str1)
	
	
	