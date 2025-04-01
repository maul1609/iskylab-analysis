import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


fileNamesOPC_M=[ \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp002-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp003-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp004-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp005-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp006-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp007-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp008-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp009-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp010-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp011-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp012-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp013-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp014-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp015-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp016-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp017-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp018-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp019-2-OPC-d-MergedW.csv', \
# 	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp020-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp021-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp022-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp023-2-OPC-d-MergedW.csv', \
# 	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp024-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp025-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp026-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp027-2-OPC-d-MergedW.csv',\
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp028-2-OPC-d-MergedW.csv', \
	'../iSKYLAB-data/Datasets/Timeseries-During-Expansion/iSKYLAB01-Exp029-2-OPC-d-MergedW.csv']

opcStr=['MergedOPC-Exp002','MergedOPC-Exp003','MergedOPC-Exp004','MergedOPC-Exp005',\
	'MergedOPC-Exp006','MergedOPC-Exp007','MergedOPC-Exp008','MergedOPC-Exp009',\
	'MergedOPC-Exp010','MergedOPC-Exp011','MergedOPC-Exp012','MergedOPC-Exp013',\
	'MergedOPC-Exp014','MergedOPC-Exp015','MergedOPC-Exp016','MergedOPC-Exp017',\
	'MergedOPC-Exp018','MergedOPC-Exp019',\
	#'MergedOPC-Exp020', 
	'MergedOPC-Exp021',\
	'MergedOPC-Exp022','MergedOPC-Exp023',\
	#'MergedOPC-Exp024',
	'MergedOPC-Exp025',\
	'MergedOPC-Exp026','MergedOPC-Exp027','MergedOPC-Exp028','MergedOPC-Exp029']

def readData(readThis = 3,opcStr="MergedOPC-Exp005"):
	time=[];pres=[];tgwm=[];tgws=[];twwm=[];twws=[];
	tdew=[];tfrost=[];cpci=[];cpctb=[];cpctt=[]
	csvfile = open(fileNamesOPC_M[readThis],'r')
	reader=csv.DictReader(csvfile) 
	
	fieldnames = reader.fieldnames;
	len1=len(fieldnames)
	Dp = [float(val) for val in fieldnames[1:-4]]
	time=[]
	conc=[]
	ntot=[]
	ndrop=[]
	nice=[]
	lwc=[]
	i=0
	for row in reader:
		time.append(float(row[fieldnames[0]]))
		conc.append([float(row[val]) for val in fieldnames[1:-4]])
		ntot.append(float(row[fieldnames[-4]]))
		ndrop.append(float(row[fieldnames[-3]]))
		nice.append(float(row[fieldnames[-2]]))
		lwc.append(float(row[fieldnames[-1]]))
		i=i+1
	csvfile.close()
	conc=np.array(conc)
	Deff=np.sum(conc*np.array(Dp)**3,axis=1) / np.sum(conc*np.array(Dp)**2,axis=1)

	data1=dict()
	data1 = {opcStr : \
		{"Time" : np.array(time), "Conc" : np.array(conc), \
		"ntot" : np.array(ntot), "ndrop": np.array(ndrop), \
		"nice" : np.array(nice), "lwc" : np.array(lwc), \
		'Deff': Deff, \
		'Dp': np.array(Dp)}}
	return data1

if __name__== "__main__":
	plotModel=False
	readThis=4
	
	if 'data1' in locals():
		pass
	else:
		data1=dict()
		
	for i in range(len(opcStr)):
		data2=readData(readThis=i,opcStr=opcStr[i])
		data1[opcStr[i]]=data2[opcStr[i]].copy()
	

	if plotModel:
		plt.pcolormesh(data1['MergedOPC-Exp005']['Time'], \
		data1['MergedOPC-Exp005']['Dp'], \
		data1['MergedOPC-Exp005']['Conc'].T, \
		norm=LogNorm(vmin=1e-3, vmax=data1['MergedOPC-Exp005']['Conc'].max()))  
