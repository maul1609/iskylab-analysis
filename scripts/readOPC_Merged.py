"""Read merged iSKYLAB OPC spectra and bulk cloud properties.

Diameters are supplied in micrometres by the CSV header.  Effective diameter
is calculated as M3/M2 and is NaN for an empty spectrum.
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import iskylab_config as cfg


fileNamesOPC_M=[ \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp002-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp003-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp004-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp005-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp006-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp007-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp008-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp009-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp010-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp011-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp012-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp013-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp014-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp015-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp016-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp017-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp018-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp019-2-OPC-d-MergedW.csv'), \
str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp020-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp021-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp022-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp023-2-OPC-d-MergedW.csv'), \
# 	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp024-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp025-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp026-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp027-2-OPC-d-MergedW.csv'),\
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp028-2-OPC-d-MergedW.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Timeseries-During-Expansion/iSKYLAB01-Exp029-2-OPC-d-MergedW.csv')]

opcStr=['MergedOPC-Exp002','MergedOPC-Exp003','MergedOPC-Exp004','MergedOPC-Exp005',\
	'MergedOPC-Exp006','MergedOPC-Exp007','MergedOPC-Exp008','MergedOPC-Exp009',\
	'MergedOPC-Exp010','MergedOPC-Exp011','MergedOPC-Exp012','MergedOPC-Exp013',\
	'MergedOPC-Exp014','MergedOPC-Exp015','MergedOPC-Exp016','MergedOPC-Exp017',\
	'MergedOPC-Exp018','MergedOPC-Exp019',\
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
	# The merged spectra are dN/dlog10D.  Reconstruct logarithmic bin widths
	# from the diameter centres before calculating number-weighted moments.
	# This avoids the historical hard-coded 0.009839 factor and remains valid if
	# a later OPC processing version changes the diameter grid.
	dp_arr=np.asarray(Dp,dtype=float)
	logc=np.log10(dp_arr)
	loge=np.empty(len(dp_arr)+1,dtype=float)
	loge[1:-1]=0.5*(logc[:-1]+logc[1:])
	loge[0]=logc[0]-0.5*(logc[1]-logc[0])
	loge[-1]=logc[-1]+0.5*(logc[-1]-logc[-2])
	dlog=np.diff(loge)
	dp_edges=10.0**loge
	nbin=conc*dlog[None,:]  # number concentration represented by each OPC bin

	# Cloud-drop bulk moments use the same 2-um lower threshold as the standard
	# iSKYLAB number comparison.  The complete spectrum is still retained.
	drop_mask=dp_arr>2.0
	w=np.where(drop_mask[None,:],nbin,0.0)
	m0=np.nansum(w,axis=1)
	m1=np.nansum(w*dp_arr[None,:],axis=1)
	m2=np.nansum(w*dp_arr[None,:]**2,axis=1)
	m3=np.nansum(w*dp_arr[None,:]**3,axis=1)
	Dmean=np.full_like(m0,np.nan,dtype=float)
	Dvol=np.full_like(m0,np.nan,dtype=float)
	Deff=np.full_like(m0,np.nan,dtype=float)
	rel_disp=np.full_like(m0,np.nan,dtype=float)
	good0=m0>0.0
	Dmean[good0]=m1[good0]/m0[good0]
	Dvol[good0]=(m3[good0]/m0[good0])**(1.0/3.0)
	good2=m2>0.0
	Deff[good2]=m3[good2]/m2[good2]
	variance=np.zeros_like(m0,dtype=float)
	variance[good0]=np.maximum(m2[good0]/m0[good0]-Dmean[good0]**2,0.0)
	good_disp=good0 & (Dmean>0.0)
	rel_disp[good_disp]=np.sqrt(variance[good_disp])/Dmean[good_disp]
	ndrop_psd=np.nansum(w,axis=1)

	data1=dict()
	data1 = {opcStr : \
		{"Time" : np.array(time), "Conc" : np.array(conc), \
		"ntot" : np.array(ntot), "ndrop": np.array(ndrop), \
		"nice" : np.array(nice), "lwc" : np.array(lwc), \
		'Deff': Deff, 'Dmean': Dmean, 'Dvol': Dvol, 'rel_disp': rel_disp, \
		'ndrop_psd': ndrop_psd, 'Dp': np.array(Dp), 'Dp_edges': dp_edges, \
		'dlogD': dlog}}
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
